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

/// Integrates fractionnal data (n.a).
/// data: raw fractionnal data (n.a)   
/// sample_rate: sampling rate (Hz) during fract acquisition   
/// returns: integrated data ((s) if input is fract. frequency)  
pub fn fractionnal_integral (data: Vec<f64>, sample_rate: f64) -> Vec<f64> {
    let dt = 1.0_f64 / sample_rate;
    let mean = statistical::mean(&data);
    // Substract mean value before cumsum
    // in order to avoir precision issues when we have
    // small frequency fluctuations on a large average frequency
    let mut data = data.clone();
    for i in 0..data.len() {
        data[i] -= mean
    }
    cumsum(data, Some(dt))
}

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
        ret.push(2.0_f64*std::f64::consts::PI * f_0 * phase[i])
    }
    ret
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_cumsum() {
        let input: Vec<f64> = vec![1.0_f64,1.0_f64,1.0_f64,1.0_f64];
        let output = cumsum(input, None);
        println!("{:#?}", output)
    }
    
    #[test]
    fn test_fractionnal_integral() {
        let input: Vec<f64> = vec![1.0_f64,1.0_f64,1.0_f64,1.0_f64];
        let output = fractionnal_integral(input, 1.0_f64);
        println!("{:#?}", output)
    }
    
    #[test]
    fn test_diff() {
        let input: Vec<f64> = vec![1.0_f64,2.0_f64,3.0_f64,4.0_f64];
        let output = diff(input, None);
        println!("{:#?}", output)
    }
    
    #[test]
    fn test_derivative() {
        let input: Vec<f64> = vec![1.0_f64,1.0_f64,1.0_f64,1.0_f64];
        let output = derivative(input, 1.0_f64);
        println!("{:#?}", output)
    }
}
