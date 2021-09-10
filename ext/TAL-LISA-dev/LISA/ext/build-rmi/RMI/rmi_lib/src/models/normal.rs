// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 

use crate::models::*;
use std::f64;

fn exp1(inp: f64) -> f64 {
    let mut x = inp;
    x = 1.0 + x / 64.0;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    x *= x;
    return x;
}

fn phi(x: f64) -> f64 {
    return 1.0 / (1.0 + exp1(-1.65451 * x));
}

fn ncdf<T: TrainingKey>(loc_data: &RMITrainingData<T>) -> (f64, f64, f64) {
    let mut scale = -f64::INFINITY;
    let mut mean = 0.0;
    let mut stdev = 0.0;

    let n = loc_data.len() as f64;

    for (inp, y) in loc_data.iter() {
        let x = inp.as_float();
        mean += x / n;
        scale = f64::max(scale, y as f64);
    }
    
    for (inp, _y) in loc_data.iter() {
        let x = inp.as_float();
        stdev += (x - mean).powf(2.0)
    }

    stdev /= n;
    stdev = stdev.sqrt();

    return (mean, stdev, scale);
}

fn lncdf<T: TrainingKey>(loc_data: &RMITrainingData<T>) -> (f64, f64, f64) {
    let mut scale = -f64::INFINITY;
    let mut mean = 0.0;
    let mut stdev = 0.0;

    let n = loc_data.len() as f64;
    for (inp, y) in loc_data.iter() {
        let x = inp.as_float();
        let lnx = if !f64::is_finite(x.ln()) { 0.0 } else { x.ln() };
        mean += lnx / n;
        scale = f64::max(scale, y as f64);
    }

    for (inp, _y) in loc_data.iter() {
        let x = inp.as_float();
        let lnx = if !f64::is_finite(x.ln()) { 0.0 } else { x.ln() };
        stdev += (lnx - mean).powf(2.0)
    }


    stdev /= n;
    stdev = stdev.sqrt();

    return (mean, stdev, scale);
}

pub struct NormalModel {
    params: (f64, f64, f64),
}

impl NormalModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> NormalModel {
        return NormalModel { params: ncdf(data) };
    }
}

impl Model for NormalModel {
    fn predict_to_float(&self, inp: &ModelInput) -> f64 {
        let (mean, stdev, scale) = self.params;
        return phi((inp.as_float() - mean) / stdev) * scale;
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![
            self.params.0.into(),
            self.params.1.into(),
            self.params.2.into(),
        ];
    }

    fn code(&self) -> String {
        return String::from(
            "
inline double ncdf(double mean, double stdev, double scale, double inp) {
    return phi((inp - mean) / stdev) * scale;
}",
        );
    }

    fn function_name(&self) -> String {
        return String::from("ncdf");
    }
    fn standard_functions(&self) -> HashSet<StdFunctions> {
        let mut to_r = HashSet::new();
        to_r.insert(StdFunctions::EXP1);
        to_r.insert(StdFunctions::PHI);
        return to_r;
    }
}

#[cfg(test)]
mod ncdf_tests {
    use super::*;

    #[test]
    fn test_ncdf1() {
        let md = ModelData::IntKeyToIntPos(vec![(1, 1), (2, 3), (3, 5)]);

        let ncdf_mod = NormalModel::new(&md);

        assert_eq!(ncdf_mod.predict_to_int(2.into()), 2);
        assert_eq!(ncdf_mod.predict_to_int(1.into()), 0);
    }

    #[test]
    fn test_empty() {
        NormalModel::new(&ModelData::empty());
    }

}

pub struct LogNormalModel {
    params: (f64, f64, f64),
}

impl LogNormalModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> LogNormalModel {
        return LogNormalModel {
            params: lncdf(data),
        };
    }
}

impl Model for LogNormalModel {
    fn predict_to_float(&self, inp: &ModelInput) -> f64 {
        let (mean, stdev, scale) = self.params;
        let data = inp.as_float();
        return phi((f64::max(data.ln(), 0.0) - mean) / stdev) * scale;
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![
            self.params.0.into(),
            self.params.1.into(),
            self.params.2.into(),
        ];
    }

    fn code(&self) -> String {
        return String::from(
            "
inline double lncdf(double mean, double stdev, double scale, double inp) {
    return phi((fmax(0.0, log(inp)) - mean) / stdev) * scale;
}",
        );
    }

    fn function_name(&self) -> String {
        return String::from("lncdf");
    }
    fn standard_functions(&self) -> HashSet<StdFunctions> {
        let mut to_r = HashSet::new();
        to_r.insert(StdFunctions::EXP1);
        to_r.insert(StdFunctions::PHI);
        return to_r;
    }
}

#[cfg(test)]
mod lncdf_tests {
    use super::*;

    #[test]
    fn test_lncdf1() {
        let md = ModelData::IntKeyToIntPos(vec![(1, 1), (2, 2), (3, 20)]);

        let lncdf_mod = LogNormalModel::new(&md);

        assert_eq!(lncdf_mod.predict_to_int(2.into()), 11);
        assert_eq!(lncdf_mod.predict_to_int(1.into()), 2);
    }

    #[test]
    fn test_empty() {
        LogNormalModel::new(&ModelData::empty());
    }

}
