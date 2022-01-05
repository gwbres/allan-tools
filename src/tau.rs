use thiserror::Error;

/// Lists all `TauAxis` known to the generator
#[derive(Clone, Copy)]
pub enum TauAxis {
    Octave, // octave axis
    Decade, // decade axis
    All,    // full linear space
}

/// `TauAxis` related errors 
#[derive(Error, Debug)]
pub enum Error {
    #[error("encountered non valid `tau` < 0 value")]
    NegativeTauValue, 
    #[error("encountered non valid `tau` == 0")]
    NullTauValue, 
    #[error("`tau` axis should only comprise increasing values")]
    InvalidTauShape,
}

impl Default for TauAxis {
    /// Builds a default `TauAxis` 
    fn default() -> TauAxis {
        TauAxis::Octave
    }
}

/// Returns Ok() if given tau axis passes standard sanity checks
pub fn tau_sanity_checks (taus: &Vec<f64>) -> Result<(), Error> {
    for i in 0..taus.len() {
        if taus[i] < 0.0_f64 {
            return Err(Error::NegativeTauValue)
        }

        if taus[i] == 0.0_f64 {
            return Err(Error::NullTauValue)
        }

        if i > 0 {
            if taus[i] <= taus[i-1] {
                return Err(Error::InvalidTauShape)
            }
        }
    }
    Ok(())
}

/// Generate log(base) `TauAxis`
/// ranging from [tau_0: tau_m]
fn log_n_tau_generator (tau_0: f64, tau_m: f64, base: f64) -> Vec<f64> {
    let mut tau = tau_0; 
    let mut ret: Vec<f64> = Vec::new();
    while tau <= tau_m {
        ret.push(tau);
        tau *= base
    }
    ret
}

/// Generates `Log2` axis ranging from [tau_0: tau_m]
fn log2_tau_generator (tau_0: f64, tau_m: f64) -> Vec<f64> { log_n_tau_generator(tau_0, tau_m, 2.0_f64) }
/// Generates `Log10` axis ranging from [tau_0: tau_m]
fn log10_tau_generator (tau_0: f64, tau_m: f64) -> Vec<f64> { log_n_tau_generator(tau_0, tau_m, 10.0_f64) }

/// Crates `tau` axis [`tau_0`: `tau_m`]    
/// `tau_0` is samling rate is standard use and adev calculations,   
/// `tau_m` < 2^32 for TauAxis::All
pub fn tau_generator (axis: TauAxis, tau_0: f64, tau_m: f64) -> Vec<f64> {
    match axis {
        TauAxis::Octave => log2_tau_generator(tau_0, tau_m),
        TauAxis::Decade => log10_tau_generator(tau_0, tau_m),
        TauAxis::All    => (tau_0 as u32..tau_m as u32)
                        .map(f64::from)
                            .collect(),
    }
}

#[cfg(test)]
    mod tests {
    use super::*;

    #[test]
    /// Tests `Tau` generator
    fn test_tau_generator() {
        // Octave axis
        let taus = tau_generator(TauAxis::Octave, 1.0_f64, 2.0_f64.powf(16.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 2.0_f64.powf(i as f64))
        }
        // Decade axis
        let taus = tau_generator(TauAxis::Decade, 1.0_f64, 10.0_f64.powf(16.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 10.0_f64.powf(i as f64))
        }
        // Full axis
        let taus = tau_generator(TauAxis::All, 1.0_f64, 1_000_000.0_f64);
        for i in 0..taus.len() {
            assert_eq!(taus[i], i as f64 +1.0)
        }
    }
}
