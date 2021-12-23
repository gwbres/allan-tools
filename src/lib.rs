//! Allan-tools 
//! 
//! This lib is a portage / equivalent package to
//! Allantools (python) library to compute Allan & related
//! statistics over some data.
//!
//! url: <https://github.com/gwbres/allan-tools>  
//! url: <https://github.com/aewallin/allantools>

mod tau;
mod utilities;

/// describes all known computations
#[derive(Clone, Copy)]
pub enum Deviation {
    Standard, // `std` deviation
    Allan,    // `allan` deviation
    Modified, // `modified` deviation
    Hadamard, // `hadamard` deviation
    Time,     // `time` deviation
    Total,    // `total` deviation
    Theo,     // `theo` deviation
}

#[cfg(test)]
mod tests {
    use super::*;
    //#[test]
    /*
     * tests lib against white noise
     */
    /*fn white_noise() {
        let mut rng = rand::thread_rng();
        let between = Range::new(0.0, 1.0);
        
        for _ in 0..10_000 {
            let v = between.ind_sample(&mut rng);
        }
    }*/
}
