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
pub enum TauAxisError {
    #[error("encountered non valid `tau` < 0 value")]
    NegativeTauValue, 
}

impl Default for TauAxis {
    /// Builds a default `TauAxis` 
    fn default() -> TauAxis {
        TauAxis::Octave
    }
}

/// Returns Ok() if given tau axis passes standard sanity checks
pub fn tau_sanity_checks (taus: Vec<f64>) -> Result<(), TauAxisError> {
    for i in 0..taus.len() {
        if taus[i] < 0.0_f64 {
            return Err(TauAxisError::NegativeTauValue)
        }
    }
    Ok(())
}

/// Crates `tau` axis [1: `tau_m`]
/// `tau_m` < 2^32 for TauAxis::All
pub fn tau_generator (axis: TauAxis, tau_m: f64) -> Vec<f64> {
    match axis {
        TauAxis::Octave => log2_tau_generator(tau_m),
        TauAxis::Decade => log10_tau_generator(tau_m),
        TauAxis::All    => (1..tau_m as u32)
                        .map(f64::from)
                            .collect(),
    }
}

/// Generate log(base) `TauAxis`
fn log_n_tau_generator (tau_m: f64, base: f64) -> Vec<f64> {
    let mut tau = 1.0_f64;
    let mut ret: Vec<f64> = Vec::new();
    while tau <= tau_m {
        ret.push(tau);
        tau *= base
    }
    ret
}

/// Generates `Log2` axis
fn log2_tau_generator (tau_m: f64) -> Vec<f64> { log_n_tau_generator(tau_m, 2.0_f64) }
/// Generates `Log10` axis
fn log10_tau_generator (tau_m: f64) -> Vec<f64> { log_n_tau_generator(tau_m, 10.0_f64) }

#[cfg(test)]
    mod tests {
    use super::*;

    #[test]
    /// Tests `Tau` generator
    fn test_tau_generator() {
        // Octave axis
        let taus = tau_generator(TauAxis::Octave, 2.0_f64.powf(16.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 2.0_f64.powf(i as f64))
        }
        // Decade axis
        let taus = tau_generator(TauAxis::Decade, 10.0_f64.powf(16.0));
        for i in 0..taus.len() {
            assert_eq!(taus[i], 10.0_f64.powf(i as f64))
        }
        // Full axis
        let taus = tau_generator(TauAxis::All, 1_000_000.0_f64);
        for i in 0..taus.len() {
            assert_eq!(taus[i], i as f64 +1.0)
        }
    }
}
