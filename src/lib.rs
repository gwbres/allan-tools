//! Allan-tools 
//! 
//! This lib is a portage / equivalent package to
//! Allantools (python) library to compute Allan & related
//! statistics over some data.
//!
//! This lib only implements basic (biased) 
//! error bar estimates at the moment.
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

/// Computes Allan variance of given input data 
/// for given tau values.  
/// data: input vector   
/// tau: input `tau` axis
/// sampling_rate: acquisition rate (Hz)   
/// is_fractionnal: true if input vector is made of fractionnal (n.a) data
/// returns: (adev, err) : deviations and error bars
pub fn avar (data: Vec<f64>, taus: Vec<f64>, sampling_rate: f64, is_fractionnal: bool) -> Result<Vec<f64>, tau::TauAxisError> {
    match tau::tau_sanity_checks(&taus) {
        Ok(_) => {},
        Err(e) => return Err(e),
    }
    
    let data = match is_fractionnal {
        false => data.clone(),
        true => utilities::fractionnal_integral(data, sampling_rate),
    };

    let tau_0 = 1.0_f64 / sampling_rate; 
    let mut devs: Vec<f64> = Vec::with_capacity(taus.len());

    for i in 0..taus.len() {
        devs[i] = calc_avar(&data, taus[i], tau_0); 
    }

    Ok(devs)
}

/// Computes Allan variance @ given tau
/// on input data sampled every tau_0 (s) [sampling period]
fn calc_avar (data: &Vec<f64>, tau: f64, tau_0: f64) -> f64 {
    let n = (tau / tau_0) as usize;
    let m = (data.len() -1)/n; 
    let mut sum = 0.0_f64;

    for i in 0..m-2 {
        sum += (data[n*i+2*n] - 2.0_f64*data[n*i+n] + data[n*i]).powf(2.0_f64)  
    }

    let norm = 2.0_f64 
        * (n as f64).powf(2.0_f64) 
            * tau_0.powf(2.0_f64) 
                * ((m-1) as f64);

    sum / norm
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn test_avar() {
        //let taus = tau::tau_generator(tau::TauAxis::Octave, 256.0_f64); 
        //let normal = rand::distributions::Normal::new(0.5, 3.0); 
        //let mut rng = rand::thread_rng();
        //let mut data: Vec<f64> = Vec::with_capacity(1024);
        //for i in 0..1024 {
        //    data[i] = rand::thread_rng().sample(normal);
        //}
        //let avar = avar(data, taus, 1.0, false);
    }
}
