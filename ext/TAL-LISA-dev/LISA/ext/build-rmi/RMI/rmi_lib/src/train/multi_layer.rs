// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 
use crate::models::*;
use crate::train::{train_model, TrainedRMI};
use log::*;

pub fn train_multi_layer(data: &mut RMITrainingData,
                         model_list: &[String],
                         last_model: String,
                         branch_factor: u64) -> TrainedRMI {
    
    let mut rmi: Vec<Vec<Box<dyn Model>>> = Vec::new();
    let mut data_partitions = vec![data.clone().into_data()];
    let num_rows = data_partitions[0].len();

    let mut current_model_count = 1;
    for (_layer_idx, model_type) in model_list.iter().enumerate() {
        info!("Training {} model layer", model_type);
        // data_partition contains all of our data partitioned into groups
        // based on the previous RMI layer's output
        let next_layer_size = current_model_count * branch_factor;
        let mut next_layer_data =
            vec![Vec::with_capacity(num_rows / next_layer_size as usize); next_layer_size as usize];
        let mut models: Vec<Box<dyn Model>> = Vec::with_capacity(next_layer_size as usize);

        for model_data in data_partitions.into_iter() {
            let mut md_container = RMITrainingData::new(&model_data);

            // not at the last layer -- rescale
            md_container.set_scale(next_layer_size as f64 / num_rows as f64);
            let model = train_model(model_type.as_str(), &md_container);

            // rescale back for next layer
            md_container.set_scale(1.0);

            for (x, y) in md_container.iter_int_int() {
                let model_pred = model.predict_to_int(x.into());
                assert!(model.needs_bounds_check() || model_pred < next_layer_size);
                let target = u64::min(next_layer_size - 1, model_pred) as usize;
                next_layer_data[target].push((x, y));
            }

            models.push(model);
        }

        data_partitions = next_layer_data
            .into_iter()
            .map(ModelData::IntKeyToIntPos)
            .collect();

        current_model_count *= branch_factor;
        rmi.push(models);
    }

    info!("Training last level {} model", last_model);
    let mut last_layer = Vec::new();
    let mut last_layer_max_l1s: Vec<u64> = Vec::new();
    let mut model_avg_error: f64 = 0.0;
    let mut model_avg_l2_error: f64 = 0.0;
    let mut model_avg_log2_error: f64 = 0.0;
    let mut model_max_log2_error: f64 = 0.0;
    let mut model_max_error = 0;
    let mut model_max_error_idx = 0;

    let mut n = 1;
    for (midx, model_data) in data_partitions.into_iter().enumerate() {
        let md_container = RMITrainingData::new(&model_data);
        let last_model = train_model(last_model.as_str(), &md_container);
        let mut max_error = 0;
        
        for (idx, (x, y)) in md_container.iter_int_int().enumerate() {
            let pred = last_model.predict_to_int(x.into());
            let err = u64::max(y, pred) - u64::min(y, pred);

            if let Some(bound) = last_model.error_bound() {
                if err > bound {
                    warn!("Precision issue: model reports max bound of {}, \
                           but an error of {} was observed on input {} at index {}. Prediction: {} Actual: {}",
                          bound, err, x, idx, pred, y);
                }
            }

            max_error = u64::max(max_error, err);
            model_avg_error += ((max_error as f64) - model_avg_error) / (n as f64);
            model_avg_l2_error += ((max_error as f64).powf(2.0) - model_avg_l2_error) / (n as f64);
            let log2_error = ((2 * max_error + 2) as f64).log2();
            model_avg_log2_error += (log2_error - model_avg_log2_error) / (n as f64);
            model_max_log2_error = f64::max(model_max_log2_error, log2_error);
            n += 1;
        }
        if max_error > model_max_error {
            model_max_error = max_error;
            model_max_error_idx = midx;
        }

        last_layer.push(last_model);
        last_layer_max_l1s.push(max_error);
    }
    rmi.push(last_layer);

    return TrainedRMI {
        model_avg_error,
        model_avg_l2_error,
        model_avg_log2_error,
        model_max_error,
        model_max_error_idx,
        model_max_log2_error,
        last_layer_max_l1s,
        rmi,
        models: format!("{},{}", model_list.join(","), last_model),
        branching_factor: branch_factor
    };
}
