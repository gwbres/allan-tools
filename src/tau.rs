/// Lists all `TauAxis` known to the generator
#[derive(Clone, Copy)]
pub enum TauAxis {
    Log2,  /* log2 axis */
    Log10, /* log10 axis */
    All,   /* complete axis */
}

impl Default for TauAxis {
    fn default() -> TauAxis {
        TauAxis::Log10
    }
}

/// Creates `tau` axis [1: `tau_m`]
/// `tau_m` < 2^32 for TauAxis::All
pub fn tau_generator (axis: TauAxis, tau_m: f64) -> Vec<f64> {
    match axis {
        TauAxis::Log2  => log2_axis(tau_m),
        TauAxis::Log10 => log10_axis(tau_m),
        TauAxis::All   => (1..tau_m as u32).map(f64::from).collect(),
    }
}

fn log2_axis (tau_m: f64) -> Vec<f64> {
    let mut tau = 1.0_f64;
    let mut ret: Vec<f64> = Vec::new();
    while tau <= tau_m {
        ret.push(tau);
        tau *= 2.0_f64;
    }
    ret
}

fn log10_axis (tau_m: f64) -> Vec<f64> {
    let mut tau = 1.0_f64;
    let mut ret: Vec<f64> = Vec::new();
    while tau <= tau_m {
        ret.push(tau);
        tau *= 10.0_f64;
    }
    ret
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /*
     * Tests Log2 tau axis tool
     */
    fn test_log2_tau_generator() {
        let taus = tau_generator(TauAxis::Log2, 2.0_f64.powf(16.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 2.0_f64.powf(i as f64))
        }
    }

    #[test]
    /*
     * Tests Log10 tau axis tool
     */
    fn test_log10_tau_generator() {
        let taus = tau_generator(TauAxis::Log10, 10.0_f64.powf(10.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 10.0_f64.powf(i as f64))
        }
    }

    #[test]
    /*
     * Tests tau axis tool
     */
    fn test_full_tau_axis() {
        let taus = tau_generator(TauAxis::All, 1_000_000.0_f64);
        for i in 0..taus.len() {
            assert_eq!(taus[i], i as f64 +1.0)
        }
    }

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
