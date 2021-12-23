//! tools / utilities to generate noise distributions

use crate::utils;
//use rand::prelude::*;

/// Generates `white` noise distribution of desired `size`
/// and desired Power Spectral Density [dBc/Hz]
pub fn white_noise (psd: f64, sample_rate: f64, size: usize) -> Vec<f64> {
    let mut rand = utils::random(size);
    //utils::normalize(rand, (psd*sample_rate/2.0_f64).powf(0.5_f64))
    rand
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_white_noise() {
        let noise = white_noise(-10.0, 10.0E3, 16);
        println!("{:#?}", noise)
    }
}
