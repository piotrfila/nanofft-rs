extern crate rustfft;
extern crate rand;
extern crate nanofft;

use rustfft::{ FftPlanner, num_complex::Complex };
use rand::{ Rng, thread_rng };
use nanofft::{ Arithmetic, fft_arrays };

fn test<T: Arithmetic, const N: usize>(data: &mut [Complex<f64>; N]) {
    let mut data_re = data.map(|x| T::from_f64(x.re));
    let mut data_im = data.map(|x| T::from_f64(x.im));
    let range = fft_arrays::<T, N>(&mut data_re, &mut data_im);
    for i in 0..N {
        data[i].re = data_re[i].into_f64(range);
        data[i].im = data_im[i].into_f64(range);
    }
}

fn test_size<const N: usize>(planner: &mut FftPlanner<f64>) -> Vec<f64> {
    let mut rng = thread_rng();
    let mut rand_complex = || Complex {
        re: rng.gen::<f64>() * 2. - 1.,
        im: rng.gen::<f64>() * 2. - 1.,
    };
    
    let data = [(); N].map(|_| rand_complex());
    let mut baseline = data;
    let rustfft = planner.plan_fft_forward(N);
    rustfft.process(&mut baseline);

    [
        test::<f32, N>,
        test::<f64, N>,
        test::<i32, N>,
        test::<i16, N>,
    ]
    .into_iter()
    .map(|f| {
        let mut data_copy = data;
        f(&mut data_copy);
        let mut diff = 0_f64;

        for (c1, c2) in data_copy.iter().zip(baseline.iter()) {
            let error = (c1.re - c2.re).powi(2) + (c1.im - c2.im).powi(2);
            let avg_magnitude = 0.5 * (c1.re.powi(2) + c1.im.powi(2) + c2.re.powi(2) + c2.im.powi(2));
            diff += error / avg_magnitude;
        }
        diff.sqrt() / (N as f64)
    })
    .collect::<Vec<_>>()
}

fn main() {
    let mut planner = FftPlanner::new();
    let results = [
        test_size::<    4>(&mut planner),
        test_size::<    8>(&mut planner),
        test_size::<   16>(&mut planner),
        test_size::<   32>(&mut planner),
        test_size::<   64>(&mut planner),
        test_size::<  128>(&mut planner),
        test_size::<  256>(&mut planner),
        test_size::<  512>(&mut planner),
        test_size::< 1024>(&mut planner),
        test_size::< 2048>(&mut planner),
        test_size::< 4096>(&mut planner),
        test_size::< 8192>(&mut planner),
        test_size::<16384>(&mut planner),
        test_size::<32768>(&mut planner),
    ];
    println!("|points|   f32   |   f64   |   i32   |   i16   |");
    println!("|-----:|:-------:|:-------:|:-------:|:-------:|");
    for (i, errors) in results.iter().enumerate() {
        print!("|{:>6}|", 2_usize.pow(i as u32 + 2));
        for e in errors {
            print!("{:9.3e}|", e);
        }
        println!();
    }
}