//! tools / utilities to generate noise distributions

use crate::utils;

/// Generates `white` noise distribution of desired `size`
/// and desired Power Spectral Density [dBc/Hz]
pub fn white_noise (psd: f64, sample_rate: f64, size: usize) -> Vec<f64> {
    let rand = utils::random(size);
    let psd = 10.0_f64.powf(psd/20.0);
    utils::normalize(rand, (2.0_f64/psd/sample_rate).powf(0.5_f64))
}

pub fn pink_noise (psd: f64, sample_rate: f64, size: usize) -> Vec<f64> {
    let mut ret: Vec<f64> = Vec::with_capacity(size);
    let white = white_noise(psd, sample_rate, size);
    let (mut b0,mut b1,mut b2,mut b3,
        mut b4,mut b5,mut b6) = 
            (0.0_f64,0.0_f64,0.0_f64,
                0.0_f64,0.0_f64,0.0_f64,0.0_f64);
    for i in 0..size {
        b0 = 0.99886 * b0 + white[i] * 0.0555179;
        b1 = 0.99332 * b1 + white[i] * 0.0750759;
        b2 = 0.96900 * b2 + white[i] * 0.1538520;
        b3 = 0.86650 * b3 + white[i] * 0.3104856;
        b4 = 0.55000 * b4 + white[i] * 0.5329522;
        b5 = -0.7616 * b5 - white[i] * 0.0168980;
        let pink = b0 + b1 + b2 + b3 + b4 + b5 + b6 + white[i] * 0.5362;
        b6 = white[i] * 0.115926;
        ret.push(pink)
    }
    ret
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;
    use crate::plotutils;

    #[test]
    fn test_white_noise_generator() {
        let samples = white_noise(-10.0, 1.0, 1000);
        assert_eq!(utils::nist_power_law_identifier(&samples), -1/2);
        plotutils::plot1d(samples, "", "White noise", "tests/white-noise.png");
    }
    #[test]
    fn test_pink_noise_generator() {
        let samples = pink_noise(-10.0, 1.0, 1000);
        assert_eq!(utils::nist_power_law_identifier(&samples), 0);
        plotutils::plot1d(samples, "", "Pink noise", "tests/pink-noise.png");
    }
}
