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

pub mod tau;
pub mod noise;
pub mod utils;

use thiserror::Error;

/// describes error related to deviation computations
#[derive(Error, Debug)]
pub enum Error {
    #[error("`tau` axis error")]
    TauAxisEror(#[from] tau::Error), 
    #[error("non feasible deviation - missing some more samples")]
    NotEnoughSamplesError, 
}

#[derive(Clone, Copy)]
/// describes all known computations
pub enum Calculation {
    Allan,    // `allan` deviation/variance
    Modified, // `modified allan` deviation/variance
    Time,     // `time` deviation/variance
}

/// Computes desired deviation of given input data 
/// for given tau values.  
/// data: input vector   
/// taus: desired `tau` offsets (s)   
/// sampling_rate: acquisition rate (Hz)   
/// is_fractionnal: true if input vector is made of fractionnal (n.a) data
/// overlapping: true if using overlapping interval (increase confidence / errbar narrows down faster)
/// returns: (adev, err) : deviations & error bars for each feasible `tau`
pub fn deviation (data: Vec<f64>, taus: &Vec<f64>, calc: Calculation, is_fractionnal: bool, overlapping: bool) 
        -> Result<(Vec<f64>,Vec<f64>), Error> 
{
    tau::tau_sanity_checks(&taus)?;
    
    let data = match is_fractionnal {
        true => utils::fractionnal_integral(data, 1.0_f64),
        false => data.clone(),
    };

    let mut devs: Vec<f64> = Vec::new();
    let mut errs: Vec<f64> = Vec::new();

    for i in 0..taus.len() {
        match calc {
            Calculation::Allan => {
                if let Ok((dev,err)) = calc_allan(&data, taus[i], false, overlapping) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Calculation::Modified => {
                if let Ok((dev,err)) = calc_modified(&data, taus[i], false) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Calculation::Time => {
                if let Ok((dev,err)) = calc_time(&data, taus[i], false) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
        }
    }
    Ok((devs, errs))
}

/// Computes desired variance at selected `taus` offset (s). 
pub fn variance (data: Vec<f64>, taus: &Vec<f64>, calc: Calculation, is_fractionnal: bool, overlapping: bool) 
        -> Result<(Vec<f64>,Vec<f64>), Error> 
{
    tau::tau_sanity_checks(&taus)?;
    
    let data = match is_fractionnal {
        true => utils::fractionnal_integral(data, 1.0_f64),
        false => data.clone(),
    };

    let mut devs: Vec<f64> = Vec::new();
    let mut errs: Vec<f64> = Vec::new();

    for i in 0..taus.len() {
        match calc {
            Calculation::Allan => {
                if let Ok((dev,err)) = calc_allan(&data, taus[i], true, overlapping) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Calculation::Modified => {
                if let Ok((dev,err)) = calc_modified(&data, taus[i], true) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Calculation::Time => {
                if let Ok((dev,err)) = calc_time(&data, taus[i], true) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
        }
    }
    Ok((devs, errs))
}
/// Computes Allan deviation/variance 
/// @ given tau on input data.   
/// overlapping: true for overlapped deviation
fn calc_allan (data: &Vec<f64>, tau: f64, is_var: bool, overlapping: bool) -> Result<(f64,f64), Error> {
    let tau_u: usize = tau as usize;
    let stride: usize = match overlapping {
        true => 1,
        false => tau as usize 
    };

    if tau_u > (data.len()-1) / 2 {
        return Err(Error::NotEnoughSamplesError)
    }

    let mut i: usize = 0;
    let mut n = 0.0_f64;
    let mut sum = 0.0_f64;

    while i < data.len() -2*tau_u {
        sum += (data[i] - 2.0_f64*data[i+tau_u] + data[i+2*tau_u]).powf(2.0_f64);
        n += 1.0_f64;
        i += stride
    }
    
    let mut dev = (sum/2.0/n);
    if !is_var {
        dev = dev.powf(0.5_f64)
    }
    dev /= tau; // * rate
    
    Ok((dev, dev/(n.powf(0.5_f64))))
}

/// Computes modified Allan deviation/variance 
/// @ given tau on input data.   
/// Mdev is always computed in overlapping fashion
fn calc_modified (data: &Vec<f64>, tau: f64, is_var: bool) -> Result<(f64,f64), Error> {
    let tau_u: usize = tau as usize;
    if tau_u > (data.len()-1) / 2 {
        return Err(Error::NotEnoughSamplesError)
    }

    let mut i: usize = 0;
    let mut n = 0.0_f64;
    let (mut v, mut sum) = (0.0_f64, 0.0_f64);

    while (i < data.len() -2*tau_u) && (i < tau_u) {
        v += data[i] - 2.0_f64*data[i+tau_u] + data[i+2*tau_u];
        i += 1
    }
    sum += v.powf(2.0_f64);
    n += 1.0_f64;

    i = 0;
    while i < data.len() -3*tau_u {
        sum += (data[i] - 3.0_f64*data[i+tau_u] + 3.0_f64*data[i+2*tau_u] - data[i+3*tau_u]).powf(2.0_f64);
        n += 1.0_f64;
        i += 1 
    }
    let mut dev = (sum/2.0/n);
    if !is_var {
        dev = dev.powf(0.5_f64)
    }
    dev /= tau; // * rate

    Ok((dev, dev/(n.powf(0.5_f64))))
}

/// Computes `time` deviation / variance
/// at desired `tau` offset (s)
fn calc_time (data: &Vec<f64>, tau: f64, is_var: bool) -> Result<(f64,f64), Error> {
    let (mdev,err) = calc_modified(data, tau, is_var)?;
    //Ok((mdev * tau / 3.0), err)
    Ok((mdev * tau / (3.0_f64).powf(0.5_f64), err))
}

/// Computes desired statistics in `Three Cornerned Hat` fashion.   
/// data_ab: A against B data   
/// data_bc: B against C data   
/// data_ca: C against A data   
/// taus: desired tau offsets
/// sample_rate: sampling rate (Hz)   
/// is_fractionnal: true if measured data are fractionnal errors   
/// overlapping: true if computing in overlapped fashion    
/// deviation: which deviation to compute    
/// returns  ((dev_a,err_a),(dev_b,err_b),(dev_c,err_c))   
/// where dev_a: deviation clock(a) and related error bar for all feasible tau offsets,
///       same thing for clock(b) and (c) 
fn three_cornered_hat(data_ab: Vec<f64>, data_bc: Vec<f64>, data_ca: Vec<f64>,
        taus: Vec<f64>, sample_rate: f64, is_fractionnal: bool, 
            overlapping: bool, calc: Calculation) 
                -> Result<((Vec<f64>,Vec<f64>), (Vec<f64>,Vec<f64>), (Vec<f64>,Vec<f64>)), Error>
{
    let (var_ab, err_ab) = variance(data_ab, &taus, calc, is_fractionnal, overlapping)?;
    let (var_bc, err_bc) = variance(data_bc, &taus, calc, is_fractionnal, overlapping)?;
    let (var_ca, err_ca) = variance(data_ca, &taus, calc, is_fractionnal, overlapping)?;
    let (mut a, mut b, mut c): (Vec<f64>,Vec<f64>,Vec<f64>) = 
        (Vec::with_capacity(var_ab.len()),Vec::with_capacity(var_ab.len()),Vec::with_capacity(var_ab.len()));
    for i in 0..var_ab.len() {
        a.push((0.5 * (var_ab[i] - var_bc[i] + var_ca[i])).powf(0.5_f64));
        b.push((0.5 * (var_bc[i] - var_ca[i] + var_ab[i])).powf(0.5_f64));
        c.push((0.5 * (var_ca[i] - var_ab[i] + var_bc[i])).powf(0.5_f64));
    }
    Ok(((a, err_ab),(b, err_bc),(c, err_ca)))
}

/// Structure optimized for `real time` / `rolling` computation,   
/// refer to dedicated documentation
pub struct RealTime {
    tau0: u64,
    buffer: Vec<f64>,
    taus: Vec<f64>,
    devs: Vec<f64>,
}

impl RealTime {

    /// Builds a new `real time` core
    pub fn new (tau_0: u64) -> RealTime {
        RealTime {
            tau0: tau_0,
            buffer: Vec::new(),
            taus: Vec::new(),
            devs: Vec::new(),
        }
    }

    /// Pushes 1 symbol into the core
    pub fn push (&mut self, sample: f64) { 
        self.buffer.push(sample); 
        if let Some((tau,dev)) = self.process() {
            self.taus.push(tau);
            self.devs.push(dev)
        }
    }
    
    /// Pushes n symbols into the core
    pub fn push_n (&mut self, samples: &Vec<f64>) { 
        for i in 0..samples.len() {
            self.buffer.push(samples[i]) 
        }
        if let Some((tau,dev)) = self.process() {
            self.taus.push(tau);
            self.devs.push(dev)
        }
    }

    /// Returns (taus, devs, errs) triplet:   
    /// taus: currently evaluated time offsets (s)   
    /// devs: related adev estimates (n.a)   
    /// errs: error bar for each adev estimate
    pub fn get (&self) -> (&Vec<f64>,&Vec<f64>,Vec<f64>) {
        let errs: Vec<f64> = Vec::with_capacity(self.devs.len());
        (&self.taus, &self.devs, errs)
    }

    /// Processes internal data & evaluates
    /// new deviation, if possible
    fn process (&self) -> Option<(f64, f64)> {
        None
    }

}

#[cfg(test)]
pub mod plotutils;

mod tests {
    use super::*;
    #[test]
    fn test_adev_whitepm() {
        let wn = noise::white_noise(-10.0, 1.0, 128);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 128.0); 
        let (dev, err)= deviation(wn.clone(), &taus, Calculation::Allan, false, false)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "ADEV", 
            vec!["ADEV (White PM)"], 
            "tests/adev-white-pm.png"
        );
    }
    #[test]
    fn test_adev_whitefm() {
        let wn = noise::white_noise(-10.0, 1.0, 1000000);
        let taus = tau::tau_generator(tau::TauAxis::Decade, 10000.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Allan, true, false)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "ADEV", 
            vec!["ADEV (White FM)"], 
            "tests/adev-white-fm.png"
        );
    }
    #[test]
    fn test_adev_pinkpm() {
        let wn = noise::pink_noise(-10.0, 1.0, 1000000);
        let taus = tau::tau_generator(tau::TauAxis::All, 100000.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Allan, false, false)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "ADEV", 
            vec!["ADEV (Pink PM)"], 
            "tests/adev-pink-pm.png"
        );
    }
    #[test]
    fn test_adev_pinkfm() {
        let wn = noise::pink_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Allan, true, false)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "ADEV", 
            vec!["ADEV (Pink FM)"], 
            "tests/adev-pink-fm.png"
        );
    }
    #[test]
    fn test_oadev_whitepm() {
        let wn = noise::white_noise(-10.0, 1.0, 128);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 128.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Allan, false, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "oADEV", 
            vec!["oADEV (White PM)"], 
            "tests/oadev-white-pm.png"
        );
    }
    #[test]
    fn test_oadev_whitefm() {
        let wn = noise::white_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Allan, true, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "oADEV", 
            vec!["oADEV (White FM)"], 
            "tests/oadev-white-fm.png"
        );
    }
/*
    #[test]
    fn test_oadev_pinkpm() {
        let wn = noise::pink_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let adev = deviation(wn.clone(), &taus, Calculation::Allan, false, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &adev.0)], 
            "oADEV", 
            vec!["oADEV (Pink PM)"], 
            "tests/oadev-pink-pm.png"
        );
    }
    #[test]
    fn test_oadev_pinkfm() {
        let wn = noise::pink_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let adev = deviation(wn.clone(), &taus, Calculation::Allan, true, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &adev.0)], 
            "oADEV", 
            vec!["oADEV (Pink FM)"], 
            "tests/oadev-pink-fm.png"
        );
    }
*/
    #[test]
    fn test_mdev_whitepm() {
        let wn = noise::white_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Modified, false, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "MDEV", 
            vec!["MDEV (White PM)"], 
            "tests/mdev-white-pm.png"
        );
    }
    #[test]
    fn test_mdev_pinkpm() {
        let wn = noise::pink_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Modified, false, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "MDEV", 
            vec!["MDEV (Pink PM)"], 
            "tests/mdev-pink-pm.png"
        );
    }
    #[test]
    fn test_mdev_whitefm() {
        let wn = noise::white_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Modified, true, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "MDEV", 
            vec!["MDEV (White FM)"], 
            "tests/mdev-white-fm.png"
        );
    }
    #[test]
    fn test_mdev_pinkfm() {
        let wn = noise::pink_noise(-10.0, 1.0, 10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let (dev, err) = deviation(wn.clone(), &taus, Calculation::Modified, true, true)
            .unwrap();
        plotutils::plot2d(
            vec![(&taus, &dev, &err)], 
            "MDEV", 
            vec!["MDEV (Pink FM)"], 
            "tests/mdev-pink-fm.png"
        );
    }
}
