// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 

use crate::models::*;
use crate::cache_fix::cache_fix;
use log::*;
use std::time::SystemTime;

mod two_layer;
//mod multi_layer;
mod lower_bound_correction;

pub struct TrainedRMI {
    pub num_rmi_rows: usize,
    pub num_data_rows: usize,
    pub model_avg_error: f64,
    pub model_avg_l2_error: f64,
    pub model_avg_log2_error: f64,
    pub model_max_error: u64,
    pub model_max_error_idx: usize,
    pub model_max_log2_error: f64,
    pub last_layer_max_l1s: Vec<u64>,
    pub rmi: Vec<Vec<Box<dyn Model>>>,
    pub models: String,
    pub branching_factor: u64,
    pub cache_fix: Option<(usize, Vec<(u64, usize)>)>,
    pub build_time: u128
}

fn train_model<T: TrainingKey>(model_type: &str,
                              data: &RMITrainingData<T>) -> Box<dyn Model> {
    let model: Box<dyn Model> = match model_type {
        "linear" => Box::new(LinearModel::new(data)),
        "robust_linear" => Box::new(RobustLinearModel::new(data)),
        "linear_spline" => Box::new(LinearSplineModel::new(data)),
        "cubic" => Box::new(CubicSplineModel::new(data)),
        "loglinear" => Box::new(LogLinearModel::new(data)),
        "normal" => Box::new(NormalModel::new(data)),
        "lognormal" => Box::new(LogNormalModel::new(data)),
        "radix" => Box::new(RadixModel::new(data)),
        "radix8" => Box::new(RadixTable::new(data, 8)),
        "radix18" => Box::new(RadixTable::new(data, 18)),
        "radix22" => Box::new(RadixTable::new(data, 22)),
        "radix26" => Box::new(RadixTable::new(data, 26)),
        "radix28" => Box::new(RadixTable::new(data, 28)),
        "bradix" => Box::new(BalancedRadixModel::new(data)),
        "histogram" => Box::new(EquidepthHistogramModel::new(data)),
        _ => panic!("Unknown model type: {}", model_type),
    };

    return model;
}

fn validate(model_spec: &[String]) {
    let num_layers = model_spec.len();
    let empty_container: RMITrainingData<u64> = RMITrainingData::empty();

    for (idx, model) in model_spec.iter().enumerate() {
        let restriction = train_model(model, &empty_container).restriction();

        match restriction {
            ModelRestriction::None => {}
            ModelRestriction::MustBeTop => {
                assert_eq!(
                    idx, 0,
                    "if used, model type {} must be the root model",
                    model
                );
            }
            ModelRestriction::MustBeBottom => {
                assert_eq!(
                    idx,
                    num_layers - 1,
                    "if used, model type {} must be the bottommost model",
                    model
                );
            }
        }
    }
}

/*fn test_rmi_input(test_key: u64, data: &RMITrainingData, rmi: &TrainedRMI) {
    let correct = data.lower_bound(test_key);
    println!("Predicting {}", test_key);
    let (guess, err) = rmi.test_predict(test_key);
    println!("Model prediction for lookup {}: {} with error {}",
             test_key, guess, err);
    
    println!("({}, {}), {}",
             guess - err,
             guess + err,
             correct);
}*/

pub fn train<T: TrainingKey>(data: &RMITrainingData<T>,
                            model_spec: &str, branch_factor: u64) -> TrainedRMI {

    let start_time = SystemTime::now();
    let (model_list, last_model): (Vec<String>, String) = {
        let mut all_models: Vec<String> = model_spec.split(',').map(String::from).collect();
        validate(&all_models);
        let last = all_models.pop().unwrap();
        (all_models, last)
    };

    if model_list.len() == 1 {
        let mut res = two_layer::train_two_layer(&mut data.soft_copy(), &model_list[0],
                                             &last_model, branch_factor);
        let build_time = SystemTime::now()
            .duration_since(start_time)
            .map(|d| d.as_nanos())
            .unwrap_or(std::u128::MAX);
        res.build_time = build_time;
        
        return res;
    }

    // it is not a simple, two layer rmi
    //return multi_layer::train_multi_layer(data, &model_list, last_model, branch_factor);
    panic!(); // TODO
}

pub fn train_for_size<T: TrainingKey>(data: &RMITrainingData<T>,
                                     max_size: usize) -> TrainedRMI {

    let start_time = SystemTime::now();
    let pareto = crate::find_pareto_efficient_configs(data, 1000);
    // go down the front until we find something small enough

    let config = pareto.into_iter()
        .filter(|x| x.size < max_size as u64)
        .next()
        .expect(format!(
            "Could not find any configurations smaller than {}", max_size).as_str());

    let models = config.models;
    let bf = config.branching_factor;

    info!("Found RMI config {} {} with size {} and average log2 {}",
          models, bf, config.size, config.average_log2_error);
    let mut res = train(data, models.as_str(), bf);
    
    let build_time = SystemTime::now()
            .duration_since(start_time)
            .map(|d| d.as_nanos())
            .unwrap_or(std::u128::MAX);
    res.build_time = build_time;
    return res;
}

pub fn train_bounded(data: &RMITrainingData<u64>,
                     model_spec: &str,
                     branch_factor: u64,
                     line_size: usize) -> TrainedRMI {
    let start_time = SystemTime::now();
    // first, transform our data into error-bounded spline points
    let spline = cache_fix(data, line_size);
    std::mem::drop(data);

    // reindex the spline points so we can build an RMI on top
    let reindexed_splines: Vec<(u64, usize)> = spline.iter()
        .enumerate()
        .map(|(idx, (key, _old_offset))| (*key, idx))
        .collect();
    
    // construct new training data from our spline points
    let mut new_data = RMITrainingData::new(Box::new(reindexed_splines));

    let mut res = crate::train(&mut new_data, model_spec, branch_factor);
    res.cache_fix = Some((line_size, spline));
    res.num_data_rows = data.len();
    
    let build_time = SystemTime::now()
        .duration_since(start_time)
        .map(|d| d.as_nanos())
        .unwrap_or(std::u128::MAX);
    res.build_time = build_time;
    return res;
}
