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
    let rand = utils::random(size);
    let psd = 10.0_f64.powf(psd/20.0);
    utils::normalize(rand, (2.0_f64/psd/sample_rate).powf(0.5_f64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;
    //use gnuplot::{Figure, Caption, Color, PointSymbol, PointSize, AxesCommon};

    #[test]
    fn test_white_noise() {
        let samples = white_noise(-10.0, 1.0, 1000);
        let m = statistical::mean(&samples);
        let mu = statistical::variance(&samples, Some(m)); 
        println!("m: {} mu: {}", m, mu);
        assert_eq!(
            vec![0], 
            utils::nist_power_law_identifier(&samples, Some(samples.len()))
        );
        
        let x: Vec<i32> = (0..samples.len() as i32).collect();
/*
        let mut fg = Figure::new();
        fg.axes2d()
            .set_x_grid(true)
            .set_y_grid(true)
            .lines(&x, &samples,
                &[Caption("White noise"),
                    Color("blue")
                 ]);
        fg.show();
*/
    }
    #[test]
    fn test_pink_noise() {
        let samples = pink_noise(-10.0, 1.0, 1000);
        let m = statistical::mean(&samples);
        let mu = statistical::variance(&samples, Some(m)); 
        println!("m: {} mu: {}", m, mu);
        assert_eq!(
            vec![0],
            utils::nist_power_law_identifier(&samples, Some(samples.len()))
        );

        let x: Vec<i32> = (0..samples.len() as i32).collect();
/*
        let mut fg = Figure::new();
        fg.axes2d()
            .set_x_grid(true)
            .set_y_grid(true)
            .lines(&x, &samples,
                &[Caption("Pink noise"),
                    Color("blue")
                 ]);
        fg.show();
*/    
    }
}
