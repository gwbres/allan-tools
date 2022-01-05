# allan-tools

[![Rust](https://github.com/gwbres/allan-tools/actions/workflows/rust.yml/badge.svg)](https://github.com/gwbres/allan-tools/actions/workflows/rust.yml)
[![crates.io](https://docs.rs/allan-tools/badge.svg)](https://docs.rs/allan-tools/badge.svg)

[![crates.io](https://img.shields.io/crates/v/allan-tools.svg)](https://crates.io/crates/allan-tools)
[![crates.io](https://img.shields.io/crates/d/allan-tools.svg)](https://crates.io/crates/allan-tools)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Allantools (python lib) portage to Rust

This library allows easy computations of Allan deviation & similar statistics.   
These statistical methods are mostly used in system stability studies.

### Variances / Deviations

Compute Allan deviation over raw data:

```rust
  use allantools::*;
  let taus = tau::generator(tau::TauAxis::Octave, 128);
  let (dev, errs) = deviation(&data, &taus, Calculation::Allan, false, false).unwrap();
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/model.png" alt="alt text" width="500"/>

This lib against stable32 on a well known signal.

### Overlapping

Improve statiscal confidence by using _overlapped_ formulae 

```rust
  let data: Vec<f64> = some_data();
  let taus = tau::generator(tau::TauAxis::Octave, 128);
  let overlapping = true;
  let (var, errs) = deviation(&data, &taus, Calculation::Allan, false, overlapping).unwrap();
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/oadev-white-pm.png" width="500"/>

### Fractionnal data
`is fractional` can be used to compute statistics over fractional
(n.a) data:

```rust
  let data: Vec<f64> = some_data();
  let taus = tau::generator(tau::TauAxis::Octave, 10000);
  let is_fractional = true;
  let ( adev, errs) = deviation(&data, &taus, Calculation::Allan, is_fractional, false).unwrap();
  let (oadev, errs) = deviation(&data, &taus, Calculation::Allan, is_fractional, true).unwrap();
```

### Tau axis generator

The user can pass any &#964; serie to all computation methods.   

This lib integrates a &#964; axis generator too, which is a convenient
method to quickly pass a standard axis to a computation method.
Several axis are known:  

* TauAxis::Octave is the most efficient
* TauAxis::Decade is the standard and is efficient
* TauAxis::All requires more computation

```rust
  let taus = tau::generator(tau::TauAxis::Decade, 10000); //log10
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/adev-white-fm.png" alt="alt text" width="500"/>

Using TauAxis::All requires more computation but gives a total
time granularity 

```rust
  let taus = tau::generator(tau::TauAxis::All, 10000);
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/adev-pink-pm.png" alt="alt text" width="500"/>

### Tau offset and error management
 
This library computes the requested statistics for all &#964; values, as long as 
$#964;(n) can be evaluated.   
If &#964; (n) cannot be evaluated, computation stops and returns all
previously evaluated offsets.

If not a single &#964; value is feasible, the lib returns Error::NotEnoughSamplesError

The user must pass a valid &#964; serie, otherwise:

* TauAxis::NullTauValue: is returned when &#964; = 0 (non sense) is requested
* TauAxis::NegativeTauValue: is return when &#964; < 0 (non physical) is requested
* TauAxis::InvalidTauShape: shape is not an increasing (not necessarily steady) shape

### Data & Noise generators

Some data generators were integrated or develpped for testing purposes:

* White noise generator

```rust
  let x = allantools::noise::white_noise(
    -140.0_f64, // dBc/Hz
    1.0_f64, // (Hz)
    10000); // 10k samples
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/white-noise.png" alt="alt text" width="200"/>

* Pink noise generator

```rust
  let x = allantools::noise::pink_noise(
    -140.0_f64, // dBc @ 1Hz
    1.0_f64, // (Hz)
    1024); // 1k samples
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/pink-noise.png" alt="alt text" width="200"/>

|  Noise |          White PM         |        Flicker PM        |   White FM   |  Flicker FM |
|:------:|:-------------------------:|:------------------------:|:------------:|:-----------:|
|  adev  |            -3/2           |            -1            |     -1/2     |      0      |
|  mdev  |             -1            |            -1            |     -1/2     |      0      |
| method | utils::diff(noise::white) | utils::diff(noise::pink) | noise::white | noise::pink |

### Power Law Identification

#### NIST LAG1D autocorrelation
[NIST Power Law identification method[[46]]](https://www.nist.gov/publications/handbook-frequency-stability-analysis)   

This macro works well on homogeneous data series (single noise process),
or series where one noise process is very dominant.

```rust
  let r = allantools::nist_lag1d_autocorr(&some_data);
```

#### Bias1 + R(n) identification method
TODO

### Three Cornered Hat

Three cornered hat fashion statistics, to estimate
a/b/c from a against b, b against c and c against a measurements.

```rust
   let a_against_b = some_measurements("a", "b");
   let b_against_c = some_measurements("b", "c");
   let c_against_a = some_measurements("c", "a");
   
   let taus = tau::tau_generator(tau::TauAxis::Octave, 10000.0);

   let ((dev_a, err_a),(dev_b,err_b),(dev_c,err_c)) =
      three_cornered_hat(&a_against_b, &b_against_c, &c_against_a,
         &taus, 1.0, false, true, Calculation::Allan).unwrap();
```

<img src="https://github.com/gwbres/allan-tools/blob/main/tests/3corner.png" alt="alt text" width="450"/>

### Tools & utilities

__cumsum__ : (python::numpy like) returns cummulative sum of a serie
```rust
   let data: Vec<f64> = some_data();
   allantools::utilities::cumsum(data, None);
   allantools::utilities::cumsum(data, Some(10E6_f64)); // opt. normalization
```

__diff__ : (python::numpy like) returns 1st order derivative of a serie
```rust
   let data: Vec<f64> = some_data();
   allantools::utilities::diff(data, None);
   allantools::utilities::diff(data, Some(10E6_f64)); // opt. normalization
```

__random__ : generates a pseudo random sequence 0 < x <= 1.0
```rust
   let data = allantools::utilities::random(1024); // 1k symbols 
   println!("{:#?}", data);
```

__normalize__ : normalizes a sequence to 1/norm :
```rust
   let data: Vec<f64> = somedata(); 
   let normalized = allantools::utilities::normalize(
       data, 
       2.0_f64 * std::f64::consts::PI); // 1/(2pi)
```

__to\_fractional\_frequency__ : converts a raw data serie
to fractional data.   
```rust
   let data: Vec<f64> = somedata(); // sampled @ 10kHz
   let fract = allantools::utilities::to_fractional_frequency(data, 10E3); // :)
```

__fractional_integral__ : converts a serie of fractional measurements
to integrated measurements (like fractional frequency (n.a) to phase time (s)).
```rust
   let data: Vec<f64> = somedata(); // (n.a) 
   let fract = allantools::utilities::fractional_integral(data, 1.0); // sampled @ 1Hz :)
```

__fractional\_freq\_to\_phase\_time__ : macro wrapper of previous function

__phase\_to\_radians__ : converts phase time (s) to phase radians (rad)
```rust
   let data: Vec<f64> = somedata(); // (s)
   let data_rad = allantools::utilities::phase_to_radians(data);
```
