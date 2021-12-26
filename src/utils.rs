//! tools / utilities to manipulate
//! datasets & vectors

use rand::prelude::*;
use rand_distr::StandardNormal;

/// numpy::cumsum direct equivalent 
pub fn cumsum (data: Vec<f64>, normalization: Option<f64>) -> Vec<f64> {
    let mut ret: Vec<f64> = Vec::with_capacity(data.len());
    ret.push(data[0]);
    for i in 1..data.len() {
        ret.push((ret[i-1] + data[i]) * normalization.unwrap_or(1.0_f64))
    }
    ret
}

/// numpy::diff direct equivalent
pub fn diff (data: Vec<f64>, normalization: Option<f64>) -> Vec<f64> {
    let mut ret: Vec<f64> = Vec::with_capacity(data.len()-1);
    for i in 1..data.len() {
        ret.push((data[i] - data[i-1]) * normalization.unwrap_or(1.0_f64))
    }
    ret    
}

/// Generate `size` random symbols  0 < x <= 1.0f
pub fn random (size: usize) -> Vec<f64> {
    let mut ret: Vec<f64> = Vec::with_capacity(size);
    for i in 0..size {
        ret.push(rand::thread_rng().sample(StandardNormal))
    }
    ret
}

/// normalizes data[] by 1/norm.   
pub fn normalize (data: Vec<f64>, norm: f64) -> Vec<f64> {
    let mut data = data.clone();
    for i in 0..data.len() {
        data[i] /= norm
    }
    data
}

/// Macro to convert frequency data (Hz) to fractionnal frequency (n.a) 
pub fn to_fractionnal_frequency (frequency: Vec<f64>, f_0: f64) -> Vec<f64> {  normalize(frequency, f_0) }

/// Integrates fractionnal data (n.a).
/// data: raw fractionnal data (n.a)   
/// sample_rate: sampling rate (Hz) during fract acquisition   
/// returns: integrated data ((s) if input is fract. frequency)  
pub fn fractionnal_integral (data: Vec<f64>, sample_rate: f64) -> Vec<f64> {
    let dt = 1.0_f64 / sample_rate;
    //let mean = statistical::mean(&data);
    // Substract mean value before cumsum
    // in order to avoir precision issues when we have
    // small frequency fluctuations on a large average frequency
    //let mut data = data.clone();
    //for i in 0..data.len() {
    //    data[i] -= mean
    //}
    cumsum(data, Some(dt))
}

/// Macro to convert fractionnal frequency data (n.a) to phase time (s) 
pub fn fractional_freq_to_phase_time (frequency: Vec<f64>, f_0: f64) -> Vec<f64> { fractionnal_integral(frequency, f_0) }

/// Computes derivative,
/// converts data to fractionnal data   
/// data: integrated data   
/// returns: fractionnal data
pub fn derivative (data: Vec<f64>, sample_rate: f64) -> Vec<f64> {
    diff(data, Some(sample_rate))
}

/// Utility function to convert
/// phase data (s) to phase (rad)   
/// phase: phase data vector    
/// f_0: norminal frequency
pub fn phase_to_radians (phase: Vec<f64>, f_0: f64) -> Vec<f64> {
    let mut ret: Vec<f64> = Vec::with_capacity(phase.len());
    for i in 0..phase.len() {
        ret.push(2.0_f64 * std::f64::consts::PI * f_0 * phase[i])
    }
    ret
}

/// Identifies power law contained in given serie.   
/// data: input data serie    
/// min_dist: minimal distance between two samples in the serie,   
///       before considering identification.
///       If none is passed: min_dist = 10    
/// returns: Vec<f64> of identified exponents, fractionnal
///          power laws are not supported @ the moment
pub fn nist_power_law_identifier (data: &Vec<f64>, min_dist: Option<usize>) -> Vec<i32> {
    let min_dist: usize = match min_dist {
        Some(d) => d,
        _ => 10
    };
    let S = data.len() / min_dist;
    let mut ret: Vec<i32> = Vec::with_capacity(S);
    for i in 0..S {
        let s = data.len() / S;
        let p = &data[i*s..(i+1)*s];
        println!("slicing size {}", p.len());
        let m = statistical::mean(&p);
        let d = 0.0_f64;
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for j in 0..p.len()-1 {
            num += (p[j] - m) * (p[j+1] - m);
            den += (p[j] - m).powf(2.0_f64);
        }
        let r = num / den;
        ret.push((-2.0*(r/(r+1.0))).round() as i32)
    }
    ret
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::noise::*;

    #[test]
    fn test_cumsum() {
        let input: Vec<f64> = vec![1.0_f64,1.0_f64,1.0_f64,1.0_f64];
        let output = cumsum(input, None);
        assert_eq!(output, vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64]);
    }
    
    #[test]
    fn test_fractionnal_integral() {
        let input: Vec<f64> = vec![1.0_f64,1.0_f64,1.0_f64,1.0_f64];
        let output = fractionnal_integral(input, 1.0_f64);
        assert_eq!(output, vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64]);
    }
    
    #[test]
    fn test_diff() {
        let input: Vec<f64> = vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64];
        let output = diff(input, None);
        assert_eq!(output, vec![1.0_f64,1.0_f64,1.0_f64]);
    }
    
    #[test]
    fn test_derivative() {
        let input: Vec<f64> = vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64];
        let output = derivative(input, 1.0_f64);
        assert_eq!(output, vec![1.0_f64,1.0_f64,1.0_f64]);
    }
    
    #[test]
    fn test_normalization() {
        let input: Vec<f64> = vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64];
        let output = normalize(input.clone(), 0.5_f64);
        assert_eq!(output, vec![2.0_f64,4.0_f64,6.0_f64,8.0_f64]);
    }

    #[test]
    fn test_power_law_identifier() {
        // pre validated noise processes 
        let data : Vec<f64> = vec!
[ 2.07764065e-01 , 3.68862105e-01 ,-3.76029189e-01 ,-2.67755620e-01, 
  6.88324037e-02 , 4.92092027e-01 ,-3.63910303e-01 , 3.74569714e-02,
 -8.84877497e-01 , 1.45030793e-01 , 1.34257424e-01 ,-3.56311376e-01,
  4.79076274e-01 , 7.98373045e-02 ,-7.86561576e-01 , 5.22930841e-01,
  5.33880123e-02 , 1.39962411e-01 , 5.07671897e-01 , 4.62672192e-01,
 -1.61125249e-01 ,-3.70219440e-02 , 1.86826848e-01 , 5.28533577e-01,
  1.16174161e+00 , 3.83249578e-01 ,-9.82527105e-01 ,-2.02268846e-01,
  6.82462721e-01 ,-1.20303678e-01 , 5.47292470e-02 , 5.31605915e-01,
 -5.29275283e-01 ,-4.77554967e-01 ,-2.00716794e-01 ,-3.79065308e-01,
 -4.33380608e-01 ,-5.42013347e-01 ,-1.38853863e-01 ,-6.99814917e-02,
 -1.18668459e+00 ,-4.72269956e-01 , 5.10866729e-01 , 9.55441552e-01,
  1.14343396e+00 , 3.05873649e-03 ,-7.84706374e-01 , 1.45764688e+00,
  1.36230710e+00 ,-3.83855214e-01 , 7.45302180e-01 , 3.07589727e-01,
 -1.49520765e+00 ,-1.60477434e+00 , 3.58264421e-02 ,-2.38223386e-01,
  2.10949298e-01 ,-6.62465088e-01 ,-1.94086696e-01 ,-9.30881468e-02,
  1.04371773e+00 ,-2.56467326e-01 ,-5.49109182e-01 , 2.71383778e-01,
 -2.40870857e-01 ,-2.01569668e-01 , 5.60955654e-01 ,-1.45863553e-02,
  7.16700996e-01 , 2.91320665e-01 , 6.81573918e-01 ,-1.06102044e+00,
  7.50204249e-01 , 2.37768696e-01 , 2.28015949e-01 , 2.06907025e-01,
 -3.68067595e-01 ,-3.89504649e-01 ,-2.47430651e-01 , 1.51238221e+00,
 -1.11786691e+00 ,-9.76931020e-01 , 9.77843564e-01 ,-9.58333797e-01,
  6.97095860e-01 , 6.32559971e-01 ,-6.76156378e-01 , 4.56158781e-01,
  3.37404000e-01 , 8.79057400e-01 , 7.49421791e-01 , 1.89524008e+00,
  4.64763523e-01 , 7.23756803e-01 , 3.06029033e-02 ,-9.78870271e-02,
 -2.82990828e-01 , 2.74069455e-01 , 9.22666572e-01 , 5.05482633e-01,
 -1.42015628e-01 ,-4.31839635e-01 , 7.73124027e-02 ,-7.56018230e-01,
 -2.40737989e-01 ,-1.02493063e-01 , 5.06388261e-01 ,-8.65145075e-01,
 -7.86580941e-02 ,-1.37397075e-01 , 5.04396892e-02 , 6.36973381e-01,
 -1.26294636e-01 ,-1.26664782e-01 ,-1.63002278e-01 , 7.64732118e-01,
  1.12082039e+00 ,-9.78517654e-01 , 2.28354542e-01 , 5.47501089e-01,
 -4.29473690e-01 ,-4.55977563e-01 ,-4.75583890e-01 ,-5.43303180e-02,
  6.57129954e-01 ,-1.52803895e-01 , 1.15782447e+00 ,-5.95232376e-01,
  8.04612608e-01 , 5.03855699e-01 , 5.26305054e-01 , 1.29038392e+00,
  8.13280500e-01 , 3.99568045e-01 ,-1.27310566e+00 ,-1.64579022e-01,
 -8.15051625e-01 ,-3.78045141e-01 ,-7.10126788e-01 ,-1.16874827e-01,
  7.06080568e-01 , 1.08007091e+00 , 2.28018776e-01 , 8.41845158e-01,
  3.86545435e-01 ,-5.63738163e-01 , 2.04115587e-02 ,-6.68489867e-01,
 -1.29312972e+00 , 8.18136371e-01 ,-1.19950723e+00 ,-1.53889093e-01,
 -1.93405375e+00 ,-5.20695613e-02 , 1.02492079e+00 ,-5.56127388e-01,
  3.24470374e-01 , 1.02689086e-01 , 8.89406488e-01 , 9.43971554e-02,
 -1.21814516e+00 ,-4.95967503e-01 ,-2.20985667e-01 ,-1.28830930e+00,
 -7.48579011e-01 ,-8.56991889e-01 ,-5.33679468e-01 , 4.95836952e-01,
  3.93049864e-01 , 7.56126529e-01 ,-1.96282807e-01 , 6.39855177e-02,
 -8.68553104e-01 , 5.40117775e-01 ,-5.44116958e-01 , 1.30773537e+00,
  1.45960388e+00 ,-1.48937475e-01 , 1.20105726e+00 ,-4.49560168e-01,
  2.00993891e-01 ,-3.79137022e-01 ,-1.16663580e+00 , 3.34511428e-01,
  4.54562168e-01 ,-1.07848098e+00 ,-3.85319167e-01 ,-6.30513668e-04,
  8.77628295e-01 ,-9.19868468e-01 ,-1.92678017e-01 , 7.43819441e-02,
  6.13739619e-01 ,-7.65085374e-01 ,-6.16063684e-01 , 1.12692180e+00,
  2.79225025e-01 ,-3.26114715e-01 , 1.01596729e+00 , 2.12728640e-01,
 -6.62698206e-01 ,-6.31746808e-01 ,-2.38516102e-01 , 3.02165723e-01,
  9.73644197e-01 , 1.12314914e-01 , 1.12048659e-01 ,-5.42811220e-02,
  9.63611722e-01 , 1.13937042e+00 , 1.21273554e+00 , 5.16357249e-01,
  2.59049527e-01 ,-1.48139997e+00 , 7.76485238e-01 , 8.62111698e-01,
  1.74430097e-01 , 6.75369813e-01 ,-1.39570298e-01 , 5.92042399e-01,
 -1.40800936e-01 , 1.03909467e+00 ,-8.23156164e-01 ,-5.09867347e-02,
 -9.61884770e-02 ,-5.80447849e-02 , 1.94808623e-01 , 1.14750476e+00,
  5.43358735e-02 ,-7.24678116e-01 ,-3.86074387e-01 ,-8.90714193e-01,
 -5.25760359e-01 ,-1.36721751e-01 , 5.91974464e-01 , 1.26329202e+00,
  5.28378972e-01 ,-7.97153538e-01 ,-5.41641081e-01 ,-4.53963689e-02,
  8.05788850e-01 , 5.76202110e-01 , 5.35887240e-01 ,-2.78022932e-01,
  2.98414976e-01 ,-1.48909782e+00 ,-1.39607468e-01 , 8.22942071e-01,
  8.43177791e-01 , 5.31933871e-02 , 2.91092717e-01 , 1.35944283e-01,
  2.66120336e-01 , 1.32896610e+00 , 5.51378198e-01 ,-1.88890741e-01];
        
		println!("POWER LAW WHITE {:#?}", nist_power_law_identifier(&data, Some(256)));
        println!("POWER LAW WHITE {:#?}", nist_power_law_identifier(&data, Some(128)));
        println!("POWER LAW WHITE {:#?}", nist_power_law_identifier(&data, Some(64)));

        let data : Vec<f64> = vec!
[ -7.66907221 , -4.13157021 , -4.42014171 , -3.98676303 , -5.65894193,  
  -7.30877706 , -6.96322893 , -4.87095624 , -6.58496316 , -6.47917415,
  -4.24285331 , -7.66441971 , -7.37749678 , -7.84193056 , -8.5759014
  -7.95451593 , -5.59906751 , -4.96521878 , -0.86901549 , -4.63224933,
  -1.1717595  ,  3.13342943 , -1.51773286 , -4.94326461 , -4.71275918,
  -3.68144799 , -3.83556567 , -2.19927632 , -5.12111462 , -8.28671313,
  -5.76004122 , -6.16818954 , -6.43084291 , -8.62011866 ,-10.43754325,
  -4.30817021 , -1.51097428 , -3.45799746 , -1.49714891 , -1.40429162,
  -1.45385766 , -1.92291151 , -2.65517297 , -4.36250499 , -5.06881069,
  -6.4675506  , -5.50916498 , -7.12342172 , -8.25830009 , -6.67526559,
  -4.71158571 , -6.62886454 , -5.02420455 , -3.27752954 , -4.66764717,
  -2.41565159 , -2.33325624 , -5.76433504 , -6.86541043 , -5.97552522,
  -4.46966194 , -4.490252   , -2.98448281 , -4.61504106 , -4.62600544,
  -5.14593469 , -5.24956138 , -6.12430194 , -7.96027992 , -4.10870057,
  -6.99710943 , -5.48648609 ,-10.66548876 , -9.49144581 , -8.54355603,
 -10.9367372  , -8.67967832 , -9.55962306 , -6.39722543 , -4.20258217,
  -4.43620008 , -7.14913031 ,-11.52870122 , -8.22233485 , -6.71806289,
  -6.7439879  , -5.20779632 , -7.97599884 , -2.90937778 , -4.6323461
  -3.38233104 , -4.1435255  , -4.56123344 , -4.59221843 , -6.63034727,
  -6.21808899 , -2.79249502 , -5.79180348 , -5.66253061 , -3.94063431,
  -2.59074194 , -4.12988808 , -3.20473602 , -5.51362572 , -7.75108616,
 -10.49820615 , -9.55922904 , -8.21928612 , -7.34711528 , -6.99684734,
  -7.61224524 , -7.87308362 , -9.72043006 , -8.37656567 , -9.98300553,
  -7.97092973 , -3.47108355 , -7.26171772 , -7.19775465 , -6.31449415,
  -5.70492654 , -7.72574286 , -8.37939271 , -7.53216658 , -3.59931843,
  -5.95299711 , -3.17798086 , -1.48560548 , -1.40122677 ,  1.09024672,
  -2.36229643 , -0.98016845 , -5.45637119 , -7.69669249 , -4.88246366,
  -5.4615448  , -6.20038134 , -1.66688459 , -6.3926793  , -7.00976557,
  -4.36633149 , -5.60700611 , -6.99129889 , -8.86050606 , -7.22980502,
  -5.48686632 , -3.39837897 , -3.80123924 , -3.61700435 , -5.40144424,
  -4.91759555 , -5.99965745 , -1.97730712 , -3.13956914 , -6.34291137,
  -3.07218343 , -4.56659077 , -3.0619463  , -0.40300699 , -1.38919313,
  -0.11799494 , -3.16166023 , -1.41730191 ,  1.33500625 ,  3.44681479,
  -1.0341273  , -1.05760556 ,  2.43371746 ,  2.5239576  , -1.52297202,
  -0.39334404 ,  1.07530775 ,  0.82953937 , -2.97009851 , -0.92303325,
   3.2848313  , -1.11401293 , -3.16757839 , -4.08763292 , -4.41510957,
  -5.91405759 , -4.68991145 , -3.6824841  , -1.55904736 , -1.00577421,
  -5.85846513 , -5.74879901 , -6.67812233 , -5.18754264 , -5.74660752,
  -4.89414559 , -0.68738116 , -1.67968166 , -3.55772947 , -2.02504749,
  -0.99115951 ,  0.6077482  , -6.36102187 , -5.88986044 , -9.60114611,
 -10.71276645 , -6.96220994 , -6.35319122 , -7.09801324 , -7.09098893,
  -6.56315725 , -9.1072789  , -9.33228824 , -6.66611163 , -6.82078187,
  -4.83453468 , -3.02105813 , -4.11149087 , -5.20411102 , -5.46345643,
  -9.14043123 , -8.31806714 , -5.90415889 , -7.77977125 , -8.61263795,
  -3.54679506 , -3.97316797 , -2.75924883 , -6.68191755 , -4.8386098
  -3.6066944  , -1.70557263 , -3.17033357 , -5.14457096 , -2.3826246
  -5.08151303 , -2.80588967 , -7.2395737  , -5.1642679  , -6.17095197,
  -7.64152626 , -3.14424759 , -7.19871501 , -4.68509663 , -6.82089387,
  -6.77008846 , -6.6572087  , -7.83152801 , -5.93162472 , -5.44302367,
  -6.74944052 , -5.91822956 , -4.06507675 , -5.67555957 , -5.02364584,
  -0.95419866 , -6.7999319  , -5.36216373 , -4.74022267 , -4.14945157, 
  -6.55556785];
		println!("POWER LAW PINK {:#?}", nist_power_law_identifier(&data, Some(256)));
        println!("POWER LAW PINK {:#?}", nist_power_law_identifier(&data, Some(128)));
        println!("POWER LAW PINK {:#?}", nist_power_law_identifier(&data, Some(64)));

    }
}
