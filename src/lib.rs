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
mod noise;
mod utils;

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
    Allan,    // `allan` deviation
    ModifiedAdev, // `modified allan` deviation
    TimeDeviation, // `tdev`
}

/// Computes desired deviation of given input data 
/// for given tau values.  
/// data: input vector   
/// taus: desired `tau` offsets (s)   
/// sampling_rate: acquisition rate (Hz)   
/// is_fractionnal: true if input vector is made of fractionnal (n.a) data
/// overlapping: true if using overlapping interval (increase confidence / errbar narrows down faster)
/// returns: (adev, err) : deviations & error bars for each feasible `tau`
pub fn deviation (data: Vec<f64>, taus: &Vec<f64>, deviation: Deviation, is_fractionnal: bool, overlapping: bool) 
        -> Result<(Vec<f64>,Vec<f64>), Error> 
{
    tau::tau_sanity_checks(&taus)?;
    
    let data = match is_fractionnal {
        false => data.clone(),
        true => utils::fractionnal_integral(data, 1.0_f64),
    };

    let mut i: u32 = 0;
    let mut devs: Vec<f64> = Vec::new();
    let mut errs: Vec<f64> = Vec::new();

    for i in 0..taus.len() {
        match deviation {
            Deviation::Allan => {
                if let Ok((dev,err)) = calc_adev(&data, taus[i], overlapping) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Deviation::ModifiedAdev => {
                if let Ok((dev,err)) = calc_mdev(&data, taus[i]) {
                    devs.push(dev);
                    errs.push(err)
                } else {
                    break
                }
            },
            Deviation::TimeDeviation => {
                if let Ok((dev,err)) = calc_tdev(&data, taus[i]) {
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

/// Computes Allan deviation @ given tau on input data.   
/// overlapping: true for overlapped deviation
fn calc_adev (data: &Vec<f64>, tau: f64, overlapping: bool) -> Result<(f64,f64), Error> {
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
        sum += (data[i + 2*tau_u] - 2.0_f64*data[i+tau_u] + data[i]).powf(2.0_f64);
        n += 1.0_f64;
        i += stride
    }
    
    sum /= 2.0_f64;
    let dev = sum.powf(0.5_f64 / n) / tau;
    Ok((dev, dev/n.powf(0.5_f64)))
}

/// Computes modified Allan deviation @ given tau on input data.   
/// Mdev is always computed in overlapping fashion
fn calc_mdev (data: &Vec<f64>, tau: f64) -> Result<(f64,f64), Error> {
    let tau_u: usize = tau as usize;
    if tau_u > (data.len()-1) / 2 {
        return Err(Error::NotEnoughSamplesError)
    }

    let mut i: usize = 0;
    let mut n = 0.0_f64;
    let (mut v, mut sum) = (0.0_f64, 0.0_f64);

    while (i < data.len() -2*tau_u) && (i < tau_u) {
        v += data[i + 2*tau_u] - 2.0_f64*data[i+tau_u] + data[i];
        i += 1
    }
    sum += v.powf(2.0_f64);
    n += 1.0_f64;

    i = 0;
    while i < data.len() -3*tau_u {
        sum += (data[i + 3*tau_u] - 3.0_f64*data[i+2*tau_u] + 3.0_f64*data[i+tau_u] - data[i]).powf(2.0_f64);
        n += 1.0_f64;
        i += 1 
    }
    sum /= 2.0_f64 * tau.powf(2.0_f64);

    let dev = sum.powf(0.5_f64 / n) / tau;
    Ok((dev, dev/n.powf(0.5_f64)))
}

/// Computes `tdev` time deviation
/// at desired `tau` offset (s)
fn calc_tdev (data: &Vec<f64>, tau: f64) -> Result<(f64,f64), Error> {
    let (mdev,err) = calc_mdev(data, tau)?;
    Ok((mdev * tau / (3.0_f64).powf(0.5_f64), err))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gnuplot::{Figure, Caption, Color, PointSymbol, PointSize, AxesCommon};

    #[test]
    fn test_adev() {
        let mut data = utils::random(100000); 
        println!("DATA \n {:#?}", data);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 1024.0); 
        let dev = deviation(data.clone(), &taus, Deviation::Allan, true, false)
            .unwrap();
        
        let x: Vec<u32> = (0..1024).collect(); 
        let mut fg = Figure::new();
        fg.set_title("ADEV");
        fg.axes2d()
            .set_x_grid(true)
            .set_x_log(Some(10.0_f64))
            .set_y_grid(true)
            .set_y_log(Some(10.0_f64))
            .lines(&taus, &dev.0, &[
                Caption("ADEV"), 
                Color("blue"), 
                PointSize(10.0),
                PointSymbol('x'),
            ]);
        
        fg.show();
        let dev = deviation(data.clone(), &taus, Deviation::ModifiedAdev, true, false)
            .unwrap();
        
        fg.axes2d()
            .set_x_grid(true)
            .set_x_log(Some(10.0_f64))
            .set_y_grid(true)
            .set_y_log(Some(10.0_f64))
            .lines(&taus, &dev.0, &[
                Caption("MDEV"), 
                Color("red"), 
                PointSize(10.0),
                PointSymbol('x'),
            ]);
        
        fg.show();
    }
/*    
    #[test]
    fn test_mdev() {
        let mut data = utils::random(1024); 
        println!("DATA \n {:#?}", data);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 128.0); 
        let dev = deviation(data, taus, Deviation::ModifiedAdev, false, false);
        println!("MDEV \n {:#?}", dev)
    }
    
    #[test]
    fn test_tdev() {
        let mut data = utils::random(1024); 
        println!("DATA \n {:#?}", data);
        let taus = tau::tau_generator(tau::TauAxis::Octave, 128.0); 
        let dev = deviation(data, taus, Deviation::TimeDeviation, false, false);
        println!("TDEV \n {:#?}", dev)
    }
*/
}
