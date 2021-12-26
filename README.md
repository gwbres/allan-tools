# allan-tools

[![crates.io](https://img.shields.io/crates/v/allan-tools.svg)](https://crates.io/crates/allan-tools)
[![Rust](https://github.com/gwbres/allan-tools/actions/workflows/rust.yml/badge.svg)](https://github.com/gwbres/allan-tools/actions/workflows/rust.yml)
[![crates.io](https://img.shields.io/crates/d/allan-tools.svg)](https://crates.io/crates/allan-tools)

Allantools (python lib) portage to Rust

This library allows easy computation of 
Allan & related statistics.   
These statistics are mostly used in system stability
studies.

### Variances / Deviations

Compute Allan Deviation over a raw data serie

```rust
  use allantools::*;
  let taus = tau::generator(tau::TauAxis::Decade, 10000);
  let (adev, errs) = deviation(&data, taus, Deviation::Allan, false, false);
```

Improve statiscal confidence by using _overlapped_ formulas 

```rust
  let taus = tau::generator(tau::TauAxis::Decade, 10000);
  let (adev, errs) = deviation(&data, taus, Deviation::Allan, false, true);
```

Compute Allan Deviation over a serie of fractionnal error

```rust
  let taus = tau::generator(tau::TauAxis::Decade, 10000);
  let ( adev, errs) = deviation(&data, taus, Deviation::Allan, true, false);
  let (oadev, errs) = deviation(&data, taus, Deviation::Allan, true, true);
```

### Integrated data & noise generators

Some data generators were integrated or develpped for testing purposes:

* White noise generator produces scaled normal distribution

```rust
  let psd = -140; // [dBcHz]
  let fs = 10.0E6; // [Hz]
  let x = allantools::noise::white_noise(psd, fs, 10000); // 10k samples
```

Some data generators were integrated or develpped for testing purposes:

* Pink noise generator produces a -10dB/dec shape when raw data is considered,
or a -5dB/dec shape if we're considering fractionnal data

```rust
  let psd = -140; // [dBcHz]
  let fs = 10.0E6; // [Hz]
  let a0_1hz = -10; // [dB] = level @ 1Hz
  let x = allantools::noise::pink_noise(a0_1hz, psd, fs, 1024); // 1k samples
```

### Tools / utilities

TODO
