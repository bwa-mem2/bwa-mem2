use log::*;
use crate::models::TrainingKey;
use crate::RMITrainingData;

#[derive(Debug)]
pub struct Spline {
    from_x: u64,
    from_y: usize,
    to_x: u64,
    to_y: usize
}

impl Spline {
    fn from(pt1: (u64, usize), pt2: (u64, usize)) -> Spline {
        assert!(pt1.0 <= pt2.0,
                "Cannot construct spline from {:?} to {:?}", pt1, pt2);
        assert!(pt1.1 <= pt2.1,
                "Cannot construct spline from {:?} to {:?}", pt1, pt2);
        return Spline { from_x: pt1.0, from_y: pt1.1,
                        to_x: pt2.0, to_y: pt2.1 };
    }

    fn with_new_dest(&self, dest: (u64, usize)) -> Spline {
        assert!(dest.0 >= self.from_x,
                "When source x is {}, cannot set dest x to {}",
                self.from_x, dest.0);
        assert!(dest.1 >= self.from_y);
        return Spline { from_x: self.from_x, from_y: self.from_y,
                        to_x: dest.0, to_y: dest.1 };
    }
    
    fn end(&self) -> (u64, usize) {
        return (self.to_x, self.to_y);
    }

    fn predict(&self, inp: u64) -> usize {
        let v0 = self.from_y as f64;
        let v1 = self.to_y as f64;
        let t = ((inp - self.from_x) as f64) / (self.to_x - self.from_x) as f64;

        return (1.0 - t).mul_add(v0, t * v1) as usize;
    }
}

struct SplineFit {
    spline: Option<Spline>,
    curr_pts: Vec<(u64, usize)>,
    line_size: usize
}

impl SplineFit {

    pub fn new(line_size: usize) -> SplineFit {
        return SplineFit {
            spline: None,
            curr_pts: Vec::new(),
            line_size
        };
    }

    pub fn add_point(&mut self, point: (u64, usize)) -> Option<(u64, usize)> {
        if self.spline.is_none() {
            self.spline = Some(Spline::from(point, point));
            return Some(point);
        }
        
        // check to see if the current spline can include this point
        let last_spline = self.spline.as_ref().unwrap();
        let proposed_spline = last_spline.with_new_dest(point);

        self.curr_pts.push(last_spline.end());
        if self.check_spline(&proposed_spline) {
            // accept this proposal, it works.
            self.spline = Some(proposed_spline);
            return None;
        } else {
            // reject this proposal, start a new spline.
            let prev_pt = last_spline.end();
            assert!(point.0 > prev_pt.0,
                    "new point: {:?} prev point: {:?}",
                    point, prev_pt);

            self.spline = Some(Spline::from(prev_pt, point));
            self.curr_pts.clear();
            self.curr_pts.push(point);
            return Some(prev_pt);
        }
    }


    pub fn finish(self) -> Option<(u64, usize)> {
        return self.spline.map(|s| s.end());
    }

    
    fn check_spline(&self, spline: &Spline) -> bool {
        return self.curr_pts.iter()
            .all(|pt| {
                let predicted_line: usize = spline.predict(pt.0) / self.line_size;
                let correct_line: usize = pt.1 / self.line_size;
                predicted_line == correct_line
            });
    }
}

pub fn cache_fix(data: &RMITrainingData<u64>, line_size: usize) -> Vec<(u64, usize)> {
    assert!(data.len() > line_size,
            "Cannot apply a cachefix with fewer items than the line size");
    info!("Fitting cachefix spline to {} datapoints", data.len());
    
    let mut fit = SplineFit::new(line_size);
    let mut spline = Vec::new();

    // Potential speedup here by carefully building a spline over the first
    // and last element of each cache line. Requires careful handling of duplicates,
    // especially when they cross cache lines.
    let mut last_key = 0;
    for (key, offset) in data.iter_unique() {
        assert!(key.minus_epsilon() >= last_key,
                "key: {:?} last key: {:?}, key - e: {:?}",
                key, last_key, key.minus_epsilon());
        
        if key.minus_epsilon() != last_key {
            match fit.add_point((key.minus_epsilon(), offset)) {
                None => {},
                Some(p) => spline.push(p)
            };
        }
        
        match fit.add_point((key, offset)) {
            None => {},
            Some(p) => spline.push(p)
        };
        
        last_key = key;
    }


    match fit.finish() {
        None => {},
        Some(p) => spline.push(p)
    };
    
    info!("Bounded spline compressed data to {}% of original ({} points, constructed from {} points).",
          ((spline.len() as f64 / data.len() as f64)*100.0).round(),
          spline.len(), data.len());

    
    return spline;
}
