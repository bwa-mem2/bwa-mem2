// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 

use crate::models::utils::{common_prefix_size, num_bits};
use crate::models::*;
use log::*;
use std::f64;

#[derive(Debug)]
pub struct BalancedRadixModel {
    params: (u8, u8, u64),
    high: bool,
}

fn chi2<T: TrainingKey>(data: &RMITrainingData<T>,
                       max_bin: u64,
                       model: &BalancedRadixModel) -> f64 {
    // compute the x^2 value of the distribution
    // induced by this model.
    let mut counts = vec![0; max_bin as usize];

    for (x, _y) in data.iter_model_input() {
        counts[model.predict_to_int(&x) as usize] += 1;
    }

    let expected = data.len() as f64 / max_bin as f64;

    return counts
        .into_iter()
        .map(|c| (c as f64 - expected).powf(2.0) / expected)
        .sum();
}

fn bradix<T: TrainingKey>(data: &RMITrainingData<T>, max_output: u64) -> BalancedRadixModel {
    let bits = num_bits(max_output);
    let common_prefix = common_prefix_size(data);
    trace!("Bradix layer common prefix: {}", common_prefix);

    let mut best_result_score = f64::INFINITY;
    let mut best_result = None;
    for test_bits in bits..u8::min(bits + 2, 64) {
        let bits_max = (1 << (test_bits + 1)) - 1;

        let high = BalancedRadixModel {
            params: (common_prefix, test_bits, max_output - 1),
            high: true,
        };
        let high_score = chi2(data, max_output, &high);

        trace!(
            "Bradix high with {} bits had score {}",
            test_bits,
            high_score
        );
        if high_score < best_result_score {
            best_result_score = high_score;
            best_result = Some(high);
        }

        let low = BalancedRadixModel {
            params: (common_prefix, test_bits, max_output - bits_max),
            high: false,
        };
        let low_score = chi2(data, max_output, &low);

        trace!("Bradix low with {} bits had score {}", test_bits, low_score);
        if low_score < best_result_score {
            best_result_score = low_score;
            best_result = Some(low);
        }
    }

    trace!(
        "Best bradix setup: {:?} with score {}",
        best_result,
        best_result_score
    );

    return best_result.unwrap();
}

impl BalancedRadixModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> BalancedRadixModel {
        if data.len() == 0 {
            return BalancedRadixModel {
                params: (0, 0, 0),
                high: true,
            };
        }

        let largest_value = data.iter().map(|(_x, y)| y).max().unwrap();

        return bradix(data, largest_value as u64);
    }
}

impl Model for BalancedRadixModel {
    fn predict_to_int(&self, inp: &ModelInput) -> u64 {
        let (left_shift, num_bits, clamp) = self.params;

        let as_int: u64 = inp.as_int();
        let res = (as_int << left_shift) >> (64 - num_bits);

        if self.high {
            return u64::min(res, clamp);
        } else {
            return if res < clamp { 0 } else { res - clamp };
        }
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Int;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Int;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![
            self.params.0.into(),
            self.params.1.into(),
            self.params.2.into(),
        ];
    }

    fn code(&self) -> String {
        if self.high {
            return String::from(
                "
inline uint64_t bradix_clamp_high(uint64_t prefix_length, 
                                  uint64_t bits, uint64_t clamp, uint64_t inp) {
    uint64_t tmp = (inp << prefix_length) >> (64 - bits);
    return (tmp > clamp ? clamp : tmp);
    
}
",
            );
        } else {
            return String::from(
                "
inline uint64_t bradix_clamp_low(uint64_t prefix_length,
                                 uint64_t bits, uint64_t clamp, uint64_t inp) {
    uint64_t tmp = (inp << prefix_length) >> (64 - bits);
    return (tmp < clamp ? 0 : tmp - clamp);
}
",
            );
        }
    }

    fn function_name(&self) -> String {
        return if self.high {
            String::from("bradix_clamp_high")
        } else {
            String::from("bradix_clamp_low")
        };
    }

    fn needs_bounds_check(&self) -> bool {
        return false;
    }
    fn restriction(&self) -> ModelRestriction {
        return ModelRestriction::MustBeTop;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        BalancedRadixModel::new(&ModelData::empty());
    }

}
