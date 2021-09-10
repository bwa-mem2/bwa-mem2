// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 
 

use crate::models::*;

fn slr<T: Iterator<Item = (f64, f64)>>(loc_data: T) -> (f64, f64) {

    // compute the covariance of x and y as well as the variance of x in
    // a single pass.

    let mut mean_x = 0.0;
    let mut mean_y = 0.0;
    let mut c = 0.0;
    let mut n: u64 = 0;
    let mut m2 = 0.0;

    let mut data_size = 0;
    for (x, y) in loc_data {
        n += 1;
        let dx = x - mean_x;
        mean_x += dx / (n as f64);
        mean_y += (y - mean_y) / (n as f64);
        c += dx * (y - mean_y);

        let dx2 = x - mean_x;
        m2 += dx * dx2;
        data_size += 1;
    }

    // special case when we have 0 or 1 items
    if data_size == 0 {
        return (0.0, 0.0);
    }

    if data_size == 1 {
        return (mean_y, 0.0);
    }


    let cov = c / ((n - 1) as f64);
    let var = m2 / ((n - 1) as f64);
    assert!(var >= 0.0, "variance of model with {} data items was negative", n);

    if var == 0.0 {
        // variance is zero. pick the mean (only) value.
        return (mean_y, 0.0);
    }

    let beta: f64 = cov / var;
    let alpha = mean_y - beta * mean_x;

    return (alpha, beta);
}

fn loglinear_slr<T: TrainingKey>(data: &RMITrainingData<T>) -> (f64, f64) {
    // log all of the outputs, omit any item that doesn't have a valid log
    let transformed_data: Vec<(f64, f64)> = data
        .iter()
        .map(|(x, y)| (x.as_float(), (y as f64).ln()))
        .filter(|(_, y)| y.is_finite())
        .collect();

    // TODO this currently creates a copy of the data and then calls
    // slr... we can probably do better by moving the log into the slr.
    return slr(transformed_data.into_iter());
}

pub struct LinearModel {
    params: (f64, f64),
}

impl LinearModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> LinearModel {
        let params = slr(data.iter()
                         .map(|(inp, offset)| (inp.as_float(), offset as f64)));
        return LinearModel { params };
    }
}

impl Model for LinearModel {
    fn predict_to_float(&self, inp: &ModelInput) -> f64 {
        let (intercept, slope) = self.params;
        return slope.mul_add(inp.as_float(), intercept);
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![self.params.0.into(), self.params.1.into()];
    }

    fn code(&self) -> String {
        return String::from(
            "
inline double linear(double alpha, double beta, double inp) {
    return std::fma(beta, inp, alpha);
}",
        );
    }

    fn function_name(&self) -> String {
        return String::from("linear");
    }

    fn set_to_constant_model(&mut self, constant: u64) -> bool {
        self.params = (constant as f64, 0.0);
        return true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear1() {
        let md = ModelData::IntKeyToIntPos(vec![(1, 2), (2, 3), (3, 4)]);

        let lin_mod = LinearModel::new(&md);

        assert_eq!(lin_mod.predict_to_int(1.into()), 2);
        assert_eq!(lin_mod.predict_to_int(6.into()), 7);
    }

    #[test]
    fn test_linear_single() {
        let md = ModelData::IntKeyToIntPos(vec![(1, 2)]);

        let lin_mod = LinearModel::new(&md);

        assert_eq!(lin_mod.predict_to_int(1.into()), 2);
    }

    #[test]
    fn test_empty() {
        LinearModel::new(&ModelData::empty());
    }

}

pub struct LogLinearModel {
    params: (f64, f64),
}

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

impl LogLinearModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> LogLinearModel {
        return LogLinearModel {
            params: loglinear_slr(&data),
        };
    }
}

impl Model for LogLinearModel {
    fn predict_to_float(&self, inp: &ModelInput) -> f64 {
        let (alpha, beta) = self.params;
        return exp1(beta.mul_add(inp.as_float(), alpha));
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![self.params.0.into(), self.params.1.into()];
    }

    fn code(&self) -> String {
        return String::from(
            "
inline double loglinear(double alpha, double beta, double inp) {
    return exp1(std::fma(beta, inp, alpha));
}",
        );
    }

    fn function_name(&self) -> String {
        return String::from("loglinear");
    }
    fn standard_functions(&self) -> HashSet<StdFunctions> {
        let mut to_r = HashSet::new();
        to_r.insert(StdFunctions::EXP1);
        return to_r;
    }
}

#[cfg(test)]
mod loglin_tests {
    use super::*;

    #[test]
    fn test_loglinear1() {
        let md = ModelData::IntKeyToIntPos(vec![(2, 2), (3, 4), (4, 16)]);

        let loglin_mod = LogLinearModel::new(&md);

        assert_eq!(loglin_mod.predict_to_int(2.into()), 1);
        assert_eq!(loglin_mod.predict_to_int(4.into()), 13);
    }

    #[test]
    fn test_empty() {
        LogLinearModel::new(&ModelData::empty());
    }
}


pub struct RobustLinearModel {
    params: (f64, f64),
}


impl RobustLinearModel {
    pub fn new<T: TrainingKey>(data: &RMITrainingData<T>) -> RobustLinearModel {
        let total_items = data.len();
        if data.len() == 0 {
            return RobustLinearModel {
                params: (0.0, 0.0)
            };
        }
        
        let bnd = usize::max(1, ((total_items as f64) * 0.0001) as usize);
        assert!(bnd*2+1 < data.len());
        
        let iter = data.iter()
            .skip(bnd)
            .take(data.len() - 2*bnd);

        let robust_params = slr(iter
                                .map(|(inp, offset)| (inp.as_float(), offset as f64)));
        
        return RobustLinearModel {
            params: robust_params
        };
    }
}

impl Model for RobustLinearModel {
    fn predict_to_float(&self, inp: &ModelInput) -> f64 {
        let (alpha, beta) = self.params;
        return beta.mul_add(inp.as_float(), alpha);
    }

    fn input_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }
    fn output_type(&self) -> ModelDataType {
        return ModelDataType::Float;
    }

    fn params(&self) -> Vec<ModelParam> {
        return vec![self.params.0.into(), self.params.1.into()];
    }

    fn code(&self) -> String {
        return String::from(
            "
inline double linear(double alpha, double beta, double inp) {
    return std::fma(beta, inp, alpha);
}",
        );
    }
    
    fn function_name(&self) -> String {
        return String::from("linear");
    }

    fn set_to_constant_model(&mut self, constant: u64) -> bool {
        self.params = (constant as f64, 0.0);
        return true;
    }
}
