# nanofft

nanofft is a dependency-free library for computing fast fourier transforms
that targets embedded systems. It provides an in-place implementation of
the Radix-2 FFT algorithm. All computations are performed directly on the
input buffer and require no additional allocations. This makes nanofft
suitable for `no_std` environments.

The main target of nanofft are microcontrollers that do not support floating
point operations, especially the rp2040. To maintain an acceptable dynamic
range, a global power-of-2 multiplier (exponent) is maintained. While computing
the fft, values are periodically bit-shifted along with incrementing or
decrementing the exponent.

## Example

The following example demonstrates computing a 16-point RFFT on a set of
samples generated from a sine signal:

```rust
// generate 16 samples of a sine wave at frequency 3
let sample_count = 16;
let signal_freq = 3.;
let sample_interval = 1. / sample_count as f32;
let mut samples: Vec<_> = (0..sample_count)
    .map(|i| {
        let phase = i as f32 * 2. * core::f32::consts::PI * signal_freq * sample_interval;
        let (im, re) = phase.sin_cos();
        (
            (i16::MAX as f32 * re) as i16,
            (i16::MAX as f32 * im) as i16,
        )
    })
    .collect();

// compute the DFT of the samples
let log2 = nanofft::i16::fft_pairs_dyn(&mut samples);

// compute magnitude and convert type
let mult = 2_f32.powi(log2 as _);
let amplitudes: Vec<_> = samples.iter()
    .map(|c| {
        let re = c.0 as f32 * mult;
        let im = c.1 as f32 * mult;
        (re.powi(2) + im.powi(2)).sqrt()
    })
    .collect();

// the spectrum has a spike at index `signal_freq`
for f in &amplitudes {
    print!("{:.2}  ", f);
}
println!();
```

### Supported FFT Sizes

nanofft only supports FFT point-sizes that are powers of two, a limitation of
the Radix-2 algorithm.

### Precison

The following table presents the RMS relative error for different data types
and FFT sizes. The test used random data as input and the comparisons were
made against rustfft using f64. Note that those results were obtained
with the `wide_trig_lut` feature enabled.

|points|   f32   |   f64   |   i16   |   i32   |
|-----:|:-------:|:-------:|:-------:|:-------:|
|     4| 3.511e-8|  0.000e0| 1.884e-2| 1.872e-2|
|     8| 3.586e-8|5.577e-11| 8.368e-3| 8.236e-3|
|    16| 3.303e-8|5.331e-11| 4.047e-3| 3.916e-3|
|    32| 2.986e-8|5.139e-11| 8.543e-4| 7.301e-4|
|    64| 2.596e-8|4.972e-11| 3.628e-4| 2.477e-4|
|   128| 2.211e-8|4.500e-11| 5.731e-4| 4.656e-4|
|   256| 1.775e-8|3.799e-11| 9.621e-5| 1.596e-9|
|   512| 1.395e-8|2.989e-11| 8.145e-5| 7.910e-6|
|  1024| 1.087e-8|2.411e-11| 6.430e-5| 2.539e-6|
|  2048| 8.510e-9|1.859e-11| 7.099e-5| 1.594e-5|
|  4096| 6.779e-9|1.474e-11| 4.506e-5|7.092e-10|
|  8192| 5.065e-9|1.151e-11| 3.480e-5|5.559e-10|
| 16384| 3.736e-9|8.596e-12| 2.688e-5|4.209e-10|
| 32768| 2.889e-9|6.640e-12| 2.089e-5|3.278e-10|

### Performance

All benchmarks were performed on a raspberry pi pico board running at 125 MHz.
The table below shows time needed per FFT calculation.

|points| rp2040 |
|-----:|:------:|
|     4| 2.61 ms|
|     8| 3.70 ms|
|    16| 5.02 ms|
|    32| 6.32 ms|
|    64| 7.93 ms|
|   128| 9.83 ms|
|   256|12.37 ms|
|   512|15.64 ms|
|  1024|20.21 ms|
|  2048|27.08 ms|
|  4096|38.53 ms|
|  8192|61.48 ms|
| 16384|130.2 ms|

## License

This project is licensed under the MIT license ([LICENSE](LICENSE) or
http://opensource.org/licenses/MIT).

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in nanofft by you, shall be licensed as above, without any
additional terms or conditions.

