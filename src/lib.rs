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
pub enum Deviation {
    Allan,    // `allan` deviation/variance
    Modified, // `modified allan` deviation/variance
    Time,     // `time` deviation/variance
}

/// Computes desired deviation of given input data 
/// for desired tau values.  
/// data: input vector   
/// taus: desired `tau` offsets (s)   
/// sampling_rate: acquisition rate (Hz)   
/// is_fractional: true if input vector is made of fractional (n.a) data
/// overlapping: true if using overlapping interval (increase confidence / errbar narrows down faster)
/// returns: (dev, err) : deviation & statistical error bars for each
/// feasible `tau`
pub fn deviation (data: &Vec<f64>, taus: &Vec<f64>, calc: Deviation, sampling_rate: f64, is_fractional: bool, overlapping: bool) 
        -> Result<(Vec<f64>,Vec<f64>), Error> 
{
    tau::tau_sanity_checks(&taus)?;
    let data = match is_fractional {
        true => utils::fractional_integral(data, 1.0_f64),
        false => data.clone(),
    };

    let mut devs: Vec<f64> = Vec::new();
    let mut errs: Vec<f64> = Vec::new();

    for i in 0..taus.len() {
        match calc {
            Deviation::Allan => {
                if let Ok((dev,err)) = calc_adev(&data, taus[i], sampling_rate, overlapping) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Deviation::Modified => {
                if let Ok((dev,err)) = calc_mdev(&data, taus[i], sampling_rate) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Deviation::Time => {
                if let Ok((dev,err)) = calc_tdev(&data, taus[i], sampling_rate) {
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

/// Computes desired variance of given input data 
/// for desired tau values.  
/// data: input vector   
/// taus: desired `tau` offsets (s)   
/// sampling_rate: acquisition rate (Hz)   
/// is_fractional: true if input vector is made of fractional (n.a) data
/// overlapping: true if using overlapping interval (increase confidence / errbar narrows down faster)
/// returns: (var, err) : variance & statistical error bars for each
/// feasible `tau`
pub fn variance (data: &Vec<f64>, taus: &Vec<f64>, dev: Deviation, sampling_rate: f64, is_fractional: bool, overlapping: bool) 
        -> Result<(Vec<f64>,Vec<f64>), Error> 
{
    let (mut var, err) = deviation(&data, &taus, dev, sampling_rate, is_fractional, overlapping)?;
    for i in 0..var.len() {
        var[i] *= var[i]
    }
    Ok((var, err))
}
/// Computes Allan deviation
/// @ given tau on input data.   
/// overlapping: true for overlapped deviation
fn calc_adev (data: &Vec<f64>, tau: f64, sampling_rate: f64, overlapping: bool) -> Result<(f64,f64), Error> {
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
        sum += (data[i+2*tau_u] - 2.0_f64*data[i+tau_u] + data[i]).powf(2.0_f64);
        n += 1.0_f64;
        i += stride
    }
    
    let mut dev = sum /2.0;
    dev = (dev / n).powf(0.5_f64) / tau / sampling_rate; 
    Ok((dev, dev/(n.powf(0.5_f64))))
}

/// Computes modified Allan deviation
/// @ given tau on input data.   
/// Mdev is always computed in overlapping fashion
fn calc_mdev (data: &Vec<f64>, tau: f64, sampling_rate: f64) -> Result<(f64,f64), Error> {
    let tau_u: usize = tau as usize;
    if tau_u > (data.len()-1) / 2 {
        return Err(Error::NotEnoughSamplesError)
    }

    let mut i: usize = 0;
    let mut n = 0.0_f64;
    let (mut v, mut sum) = (0.0_f64, 0.0_f64);

    while (i < data.len() -2*tau_u) && (i < tau_u) {
        v += data[i+2*tau_u] - 2.0_f64*data[i+tau_u] + data[i];
        i += 1
    }
    sum += v.powf(2.0_f64);
    n += 1.0_f64;

    i = 0;
    while i < data.len() -3*tau_u {
        v += data[i+3*tau_u] - 3.0_f64*data[i+2*tau_u] + 3.0_f64*data[i+tau_u] - data[i];
        sum += v.powf(2.0_f64);
        n += 1.0_f64;
        i += 1 
    }
    let mut dev = sum /2.0 /tau /tau;
    dev = (dev / n).powf(0.5_f64) / tau / sampling_rate;
    Ok((dev, dev/(n.powf(0.5_f64))))
}

/// Computes `time` deviation at desired `tau` offset (s)
fn calc_tdev (data: &Vec<f64>, tau: f64, sampling_rate: f64) -> Result<(f64,f64), Error> {
    let (mdev,mderr) = calc_mdev(data, tau, sampling_rate)?;
    Ok((
        mdev * tau / (3.0_f64).powf(0.5_f64),
        mderr // mderr / ns.powf(0.5_f64)
    ))
}

/// Computes desired statistics in `Three Cornerned Hat` fashion.   
/// data_ab: A against B data   
/// data_bc: B against C data   
/// data_ca: C against A data   
/// taus: desired tau offsets
/// sample_rate: sampling rate (Hz)   
/// is_fractional: true if measured data are fractional errors   
/// overlapping: true if computing in overlapped fashion    
/// deviation: which deviation to compute    
/// returns  ((dev_a,err_a),(dev_b,err_b),(dev_c,err_c))   
/// where dev_a: deviation clock(a) and related error bar for all feasible tau offsets,
///       same thing for clock(b) and (c) 
pub fn three_cornered_hat(data_ab: &Vec<f64>, data_bc: &Vec<f64>, data_ca: &Vec<f64>,
        taus: &Vec<f64>, sample_rate: f64, is_fractional: bool, 
            overlapping: bool, calc: Deviation) 
                -> Result<((Vec<f64>,Vec<f64>), (Vec<f64>,Vec<f64>), (Vec<f64>,Vec<f64>)), Error>
{
    let (var_ab, err_ab) = variance(&data_ab, &taus, calc, sample_rate, is_fractional, overlapping)?;
    let (var_bc, err_bc) = variance(&data_bc, &taus, calc, sample_rate, is_fractional, overlapping)?;
    let (var_ca, err_ca) = variance(&data_ca, &taus, calc, sample_rate, is_fractional, overlapping)?;
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
#[derive(Debug)]
pub struct RealTime {
    tau0: u64,
    buffer: Vec<f64>,
    taus: Vec<usize>,
    devs: Vec<f64>,
}

impl RealTime {
    /// Builds a new `real time` core.   
    /// tau_0: sampling period (s)
    pub fn new (tau_0: u64) -> RealTime {
        RealTime {
            tau0: tau_0,
            buffer: Vec::new(),
            taus: vec![tau_0 as usize],
            devs: Vec::new(),
        }
    }

    /// Pushes 1 symbol into the core
    pub fn push (&mut self, sample: f64) { 
        self.buffer.push(sample); 
        if self.new_tau() {
            self.taus.push(
                self.taus[self.taus.len()-1] * 2);
        }
    }
    
    /// Pushes n symbols into the core
    pub fn push_n (&mut self, samples: &Vec<f64>) { 
        for i in 0..samples.len() {
            self.buffer.push(samples[i]) 
        }
        if self.new_tau() {
            self.taus.push(
                self.taus[self.taus.len()-1] * 2);
        }
    }

    /// Returns (taus, devs, errs) triplet:   
    /// taus: currently evaluated time offsets (s)   
    /// devs: related adev estimates (n.a)   
    /// errs: error bar for each adev estimate
    pub fn get (&self) -> (&Vec<usize>,&Vec<f64>,Vec<f64>) {
        let errs: Vec<f64> = Vec::with_capacity(self.devs.len());
        (&self.taus, &self.devs, errs)
    }

    fn new_tau (&self) -> bool {
        self.buffer.len() >= 2*self.taus[self.taus.len()-1]+1
    }

    /// Processes internal data & evaluates
    /// new deviation, if possible
    fn process (&self) -> Option<(f64, f64)> {
        None
    }

    fn calc_estimator (&self, tau: f64) -> f64 {
        let tau_u = tau as usize;
        let mut sum = 0.0_f64;
        for i in 0..self.buffer.len()-2*tau_u {
            sum += (self.buffer[i] - 2.0_f64*self.buffer[i+tau_u] - self.buffer[i+2*tau_u]).powf(2.0_f64)
        }
        sum
    }

    /// Resets internal core
    pub fn reset (&mut self) {
        self.buffer.clear();
        self.taus.clear();
        self.devs.clear()
    }

}

#[cfg(test)]
pub mod plotutils;
mod tests {
    use super::*;
	use std::str::FromStr;
    #[test]
    fn test_deviation() {
        let N: usize = 10000;
        let noises: Vec<&str> = vec![
            "whitepm",
            "flickpm",
            "whitefm",
            "pinkfm",
        ];
        let axes: Vec<tau::TauAxis> = vec![
            tau::TauAxis::Octave,
            tau::TauAxis::Decade,
            tau::TauAxis::All,
        ];
        let calcs: Vec<Deviation> = vec![
            Deviation::Allan,
            Deviation::Modified,
            Deviation::Time,
        ];
        // test against pure noise
        for noise in noises {
            let mut input: Vec<f64>;
            if noise.eq("white") {
                input = noise::white_noise(-10.0,1.0, N)
            } else {
                input = noise::pink_noise(-10.0,1.0, N)
            };
            let is_fract = noise.contains("fm");

            for ax in &axes {
                let taus = tau::tau_generator(*ax, 1.0, 1000.0); 
                for overlapping in vec![false, true] {
                    for calc in &calcs {
                        let (dev, err) = deviation(
                            &input,
                            &taus,
                            *calc,
                            1.0_f64,
                            is_fract,
                            overlapping)
                                .unwrap();
                        let mut fp = String::from("tests/");
                        fp.push_str(noise);
                        fp.push_str("-");
                        if overlapping {
                            fp.push_str("o")
                        }
                        match calc {
                            Deviation::Allan => fp.push_str("adev"),
                            Deviation::Modified => fp.push_str("mdev"),
                            Deviation::Time => fp.push_str("tdev"),
                        }
                        fp.push_str(".png");
                        plotutils::plot1d_err(
                            vec![(&taus, &dev, &err)],
                                "test deviation",
                                vec![&fp],
                                &fp,
                        );
                    }
                }
            }
        }
    }
    /*
    #[test]
    fn test_against_models() {
       let names: Vec<&str> = vec!["adev","mdev","tdev"];
        let (mut xm_adev, mut ym_adev): (Vec<f64>,Vec<f64>) = (Vec::new(),Vec::new());
        let (mut xm_mdev, mut ym_mdev): (Vec<f64>,Vec<f64>) = (Vec::new(),Vec::new());
        let (mut xm_tdev, mut ym_tdev): (Vec<f64>,Vec<f64>) = (Vec::new(),Vec::new());
        // read models
        for i in 0..names.len() {
            let content = std::fs::read_to_string(
                std::path::PathBuf::from(
                    env!("CARGO_MANIFEST_DIR").to_owned()
                    +"/tests/holdover-" +names[i] +".csv"))
                .unwrap();
            let mut lines = content.lines();
            let mut line = lines.next()
                .unwrap();
            loop {
                let c = line.find(",").unwrap();
                let x = f64::from_str(line.split_at(c).0.trim()).unwrap();
                let y = f64::from_str(line.split_at(c+1).1.trim()).unwrap();

                if names[i].eq("adev") {
                    xm_adev.push(x);
                    ym_adev.push(y)
                } else if names[i].eq("mdev") {
                    xm_mdev.push(x);
                    ym_mdev.push(y);
                } else {
                    xm_tdev.push(x);
                    ym_tdev.push(x)
                }

                if let Some(l) = lines.next() {
                    line = l;
                } else {
                    break
                }
            }
        }
        // read raw
        let mut raw_data: Vec<f64> = Vec::new();
        let mut frac_raw_data: Vec<f64> = Vec::new();
        let mut x_adev: Vec<f64> = Vec::new();
        let mut y_adev: Vec<f64> = Vec::new();
        let names: Vec<&str> = vec!["holdover"];
        for i in 0..names.len() {
            let content = std::fs::read_to_string(
                std::path::PathBuf::from(
                    env!("CARGO_MANIFEST_DIR").to_owned()
                    +"/tests/" +names[i] +"-phase.csv"))
                .unwrap();
            let mut lines = content.lines();
            let mut line = lines.next()
                .unwrap();
            loop {
                let raw = f64::from_str(line.trim()).unwrap();
                raw_data.push(raw);
                if let Some(l) = lines.next() {
                    line = l;
                } else {
                    break
                }
            }
        }
        let taus = tau::tau_generator(tau::TauAxis::Decade, 1.0_f64, 1000000.0_f64);
        let (adev, err) = deviation(&raw_data, &taus, Deviation::Allan, 0.1_f64, false, true).unwrap();
        let (mdev, err) = deviation(&raw_data, &taus, Deviation::Modified, 0.1_f64, false, true).unwrap();

        let taus = tau::tau_generator(tau::TauAxis::Decade, 0.1_f64, 100000.0_f64);
        plotutils::plotmodel_tb(
            // models
            &xm_adev, 
            &ym_adev,
            &xm_mdev,
            &ym_mdev,
            // adev from phase
            &taus, &adev, &mdev,
            )
    }*/
    #[test]
    fn test_three_cornered_hat() {
        let pm_pink  = utils::diff(&noise::pink_noise(-10.0,1.0,10000),None);
        let fm_white = noise::white_noise(-10.0,1.0,10000);
        let fm_pink = noise::pink_noise(-10.0,1.0,10000);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1.0_f64, 10000.0);

        let ((dev_a, err_a),(dev_b,err_b),(dev_c,err_c)) =
            three_cornered_hat(&pm_pink, &fm_white, &fm_pink,
                &taus, 1.0, false, true, Deviation::Allan).unwrap();

        let (adev_ab, err_ab) = deviation(&pm_pink , &taus, Deviation::Allan, 1.0_f64, false, true).unwrap();
        let (adev_bc, err_bc) = deviation(&fm_white, &taus, Deviation::Allan, 1.0_f64, false, true).unwrap();
        let (adev_ca, err_ca) = deviation(&fm_pink , &taus, Deviation::Allan, 1.0_f64, false, true).unwrap();

        plotutils::plot3corner(
            &taus,
            (&adev_ab, &err_ab),
            (&dev_a,   &err_a),
            (&adev_bc, &err_bc),
            (&dev_b,   &err_b),
            (&adev_ca, &err_ca),
            (&dev_c,   &err_c),
        );
    }

    #[test]
    fn test_realtime_core() {
        let mut rt = RealTime::new(1);
        let noise = noise::white_noise(1.0,1.0,100);
        for i in 0..noise.len() {
            rt.push(noise[i]);
            println!("{:#?}", rt)
        }
    }
}
