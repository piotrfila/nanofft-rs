# nanofft

nanofft is a library for computing fast fourier transforms that targets
embedded systems. It provides an in-place implementation of the Radix-2 FFT
algorithm. All computations are performed directly on the input buffer and
require no additional allocations. This makes nanofft suitable for `no_std`
environments.

The main target of nanofft are microcontrollers that do not support floating
point operations. To do this while maintaining an acceptable dynamic range,
a global power-of-2 multiplier (exponent) is maintained. While computing
the fft, values are periodically bit-shifted along with incrementing
or decrementing the exponent.

### Supported FFT Sizes

nanofft only supports FFT point-sizes that are powers of two, a limitation of
the Radix-2 algorithm.

### Precison

The following table presents the RMS relative error for different data types
and FFT sizes. The test used random data as input and the comparisons were
made against rustfft using f64. Note that those results were obtained
with the `wide_trig_lut` feature enabled.

|points|   f32   |   f64   |   i32   |   i16   |
|-----:|:-------:|:-------:|:-------:|:-------:|
|    64| 2.235e-8|3.847e-11| 6.964e-9| 4.448e-4|
|   128| 2.131e-8|4.960e-11| 2.258e-9| 1.377e-4|
|   256| 1.351e-8|2.940e-11| 1.927e-9| 1.213e-4|
|   512| 1.397e-8|3.094e-11| 1.561e-9| 9.269e-5|
|  1024| 1.046e-8|2.243e-11| 1.070e-9| 6.819e-5|
|  2048| 6.842e-9|1.736e-11| 1.041e-9| 6.873e-5|
|  4096| 5.679e-9|1.175e-11| 1.035e-9| 6.619e-5|
|  8192| 5.028e-9|1.051e-11|6.627e-10| 4.185e-5|
| 16384| 4.683e-9|9.411e-12|7.535e-10| 4.190e-5|
| 32768| 5.521e-9|7.424e-12|4.147e-10| 4.240e-5|

## License

This project is licensed under the MIT license ([LICENSE](LICENSE) or
http://opensource.org/licenses/MIT).

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in microfft by you, shall be licensed as above, without any
additional terms or conditions.

