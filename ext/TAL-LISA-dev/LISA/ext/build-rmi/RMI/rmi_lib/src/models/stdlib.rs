// < begin copyright > 
// Copyright Ryan Marcus 2020
// 
// See root directory of this project for license terms.
// 
// < end copyright > 
 

#[derive(Debug, PartialEq, Eq, Hash)]
pub enum StdFunctions {
    EXP1,
    PHI,
    BinarySearch,
}

impl StdFunctions {
    pub fn decl(&self) -> &'static str {
        match self {
            StdFunctions::EXP1 => "inline double exp1(double x);",
            StdFunctions::PHI => "inline double phi(double x);",
            StdFunctions::BinarySearch => {
                "uint64_t bs_lower_bound(const uint64_t a[], uint64_t n, uint64_t x);"
            }
        }
    }

    pub fn code(&self) -> &'static str {
        match self {
            StdFunctions::EXP1 => {
                "
inline double exp1(double x) {
  x = 1.0 + x / 64.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}
"
            }
            StdFunctions::PHI => {
                "
inline double phi(double x) {
  return 1.0 / (1.0 + exp1(- 1.65451 * x));
}
"
            }
            StdFunctions::BinarySearch => {
                "
uint64_t bs_upper_bound(const uint64_t a[], uint64_t n, uint64_t x) {
    int l = 0;
    int h = n; // Not n - 1
    while (l < h) {
        int mid = (l + h) / 2;
        if (x >= a[mid]) {
            l = mid + 1;
        } else {
            h = mid;
        }
    }
    return l;
}

"
            }
        }
    }
}
