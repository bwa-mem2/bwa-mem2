// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 
 
use crate::models::*;
use crate::models::utils::radix_index;
use superslice::*;
use log::*;

pub struct EquidepthHistogramModel {
    params: Vec<u64>,
    radix: Vec<u64>
}


fn equidepth_histogram<T: TrainingKey>(data: &RMITrainingData<T>) -> Vec<u64> {
    assert!(data.len() > 0);
    
    let mut splits: Vec<u64> = Vec::new();
    let num_bins = data.get(data.len()-1).1 as usize;
    let items_per_bin = data.len() / num_bins;

    assert!(items_per_bin >= 1, "not enough items for equidepth histogram");
    info!("Equidepth histogram using {} bins", num_bins);
    
    for bin_idx in 0..num_bins {
        let start_idx = bin_idx * items_per_bin;
        let start_val = data.get_key(start_idx).as_uint();
        splits.push(start_val);
    }

    return splits;
}



impl EquidepthHistogramModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> EquidepthHistogramModel {
        if data.len() == 0 {
            return EquidepthHistogramModel { params: Vec::new(), radix: Vec::new() };
        }

        let params = equidepth_histogram(data);
        let radix = radix_index(&params, 20);
        return EquidepthHistogramModel {
            params, radix
        };
    }
}

impl Model for EquidepthHistogramModel {

    fn predict_to_int(&self, inp: &ModelInput) -> u64 {
        let val = inp.as_int();
        
        let val = self.params.upper_bound(&val) - 1;
        return val as u64;

        /*for (idx, &split) in self.params.iter().enumerate() {
            if val <= split { return (idx - 1) as u64; }
        }

        return self.params.len() as u64 - 1;*/
    }

    fn input_type(&self) -> ModelDataType { return ModelDataType::Int; }
    fn output_type(&self) -> ModelDataType { return ModelDataType::Int; }

    fn params(&self) -> Vec<ModelParam> {
        return vec![
            ModelParam::Int(self.params.len() as u64),
            ModelParam::IntArray(self.radix.clone()),
            ModelParam::IntArray(self.params.clone())
        ];
    }
    fn code(&self) -> String {
        return String::from("
inline uint64_t ed_histogram(const uint64_t length,
                             const uint64_t radix[], 
                             const uint64_t pivots[], 
                             uint64_t key) {
    uint64_t key_radix = key >> (64 - 20);
    unsigned int radix_lb = radix[key_radix];
    unsigned int radix_ub = radix[key_radix+1];
    uint64_t li = bs_upper_bound(pivots + radix_lb, radix_ub - radix_lb, key) + radix_lb - 1;
    return li;
}
");
    }

    fn standard_functions(&self) -> HashSet<StdFunctions> {
        let mut to_r = HashSet::new();
        to_r.insert(StdFunctions::BinarySearch);
        return to_r;
    }

    fn function_name(&self) -> String { return String::from("ed_histogram"); }
    fn restriction(&self) -> ModelRestriction { return ModelRestriction::MustBeTop; }
    fn needs_bounds_check(&self) -> bool { return false; }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_ed_hist1() {
        let mut test_data: Vec<(u64, u64)> = Vec::new();

        for i in 0..1000 {
            test_data.push((i*3, i/3));
        }
        
        let md = ModelData::IntKeyToIntPos(test_data);

        let ed_mod = EquidepthHistogramModel::new(&md);

        assert_eq!(ed_mod.predict_to_int((0).into()), 0);
        assert_eq!(ed_mod.predict_to_int((1*3).into()), 0);
        assert_eq!(ed_mod.predict_to_int((4*3).into()), 1);
        assert_eq!(ed_mod.predict_to_int((500*3).into()), 166);
        assert_eq!(ed_mod.predict_to_int((5000*3).into()), 333);
    }

    #[test]
    fn test_empty() {
        EquidepthHistogramModel::new(&ModelData::empty());
    }

}
