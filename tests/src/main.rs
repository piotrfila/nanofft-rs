#![feature(const_mut_refs)]
#![feature(const_swap)]

extern crate rustfft;
extern crate rand;
extern crate nanofft;

use rustfft::FftPlanner;

trait Convert where Self: Sized {
    type RangeInfo;

    fn from_f64(x: f64) -> Self;
    fn into_f64(self, s: Self::RangeInfo) -> f64;
}

macro_rules! convert_float {
    ($($t:ty)*) => { $(

    impl Convert for $t {
        type RangeInfo = ();
    
        fn from_f64(x: f64) -> Self { x as Self }
        fn into_f64(self, _: Self::RangeInfo) -> f64 { self as f64 }
    }
    
    )* };
}

macro_rules! convert_int {
    ($($t:ty)*) => { $(

    impl Convert for $t {
        type RangeInfo = i16;
    
        fn from_f64(x: f64) -> Self { (Self::MAX as f64 * x) as Self }
        fn into_f64(self, i: Self::RangeInfo) -> f64 {
            let scale = f64::from_bits(((i as i32 + 2 - f64::MIN_EXP) as u64) << (f64::MANTISSA_DIGITS - 1));
            self as f64 * scale
        }
    }
    
    )* };
}

convert_float!(f32 f64);
convert_int!(i8 i16 i32 i64);

macro_rules! mktest {
    ($($namespace:ident)*; $($size:literal)*) => {
        $( mod $namespace {
            use super::Convert;
            use rustfft::{ FftPlanner, num_complex::Complex };
            use rand::{ Rng, thread_rng };

            fn test_pairs<const N: usize>(data: &mut [Complex<f64>; N]) {
                let mut data_re = data.map(|x| <_ as Convert>::from_f64(x.re));
                let mut data_im = data.map(|x| <_ as Convert>::from_f64(x.im));
                let range = nanofft::$namespace::fft_arrays(&mut data_re, &mut data_im);
                for (dst, src) in data.iter_mut().zip(data_re.iter().zip(data_im.iter())) {
                    dst.re = src.0.into_f64(range);
                    dst.im = src.1.into_f64(range);
                }
            }
            
            fn test_arrays<const N: usize>(data: &mut [Complex<f64>; N]) {
                let mut data_t = data.map(|x| (<_ as Convert>::from_f64(x.re), <_ as Convert>::from_f64(x.im)));
                let range = nanofft::$namespace::fft_pairs(&mut data_t);
                for (dst, src) in data.iter_mut().zip(data_t.iter()) {
                    dst.re = src.0.into_f64(range);
                    dst.im = src.1.into_f64(range);
                }
            }

            fn test_arrays_dyn<const N: usize>(data: &mut [Complex<f64>; N]) {
                let mut data_t = data.map(|x| (<_ as Convert>::from_f64(x.re), <_ as Convert>::from_f64(x.im)));
                let range = nanofft::$namespace::fft_pairs_dyn(&mut data_t);
                for (dst, src) in data.iter_mut().zip(data_t.iter()) {
                    dst.re = src.0.into_f64(range);
                    dst.im = src.1.into_f64(range);
                }
            }

            pub fn test_size<const N: usize>(planner: &mut FftPlanner<f64>) -> f64 {
                let mut rng = thread_rng();
                let mut rand_complex = || Complex {
                    re: rng.gen::<f64>() * 2. - 1.,
                    im: rng.gen::<f64>() * 2. - 1.,
                };
                
                let data = [(); N].map(|_| rand_complex());
                let mut baseline = data;
                let rustfft = planner.plan_fft_forward(N);
                rustfft.process(&mut baseline);
            
                let results = [
                    test_pairs::<N>,
                    test_arrays::<N>,
                    test_arrays_dyn,
                ]
                .map(|f| {
                    let mut data_copy = data;
                    f(&mut data_copy);
                    let mut diff = 0_f64;
            
                    for (c1, c2) in data_copy.iter().zip(baseline.iter()) {
                        let error = (c1.re - c2.re).powi(2) + (c1.im - c2.im).powi(2);
                        let magnitude = c2.re.powi(2) + c2.im.powi(2);
                        diff += error / magnitude;
                    }
                    diff.sqrt() / (N as f64)
                });
                for x in &results[1..] {
                    if *x != results[0] {
                        println!("bruh");
                    }
                }
                results[0]
            }
        } )* 

        fn run_tests(repeats: usize) {
            let mut planner = FftPlanner::new();

            fn run_tests_for_size<const N: usize>(planner: &mut FftPlanner<f64>, repeats: usize) {
                print!("|{:>6}|", N);
                $({
                    let mut e = 0.;
                    for _ in 0..repeats {
                        e += $namespace::test_size::<N>(planner);
                    }
                    print!("{:9.3e}|", e / (repeats as f64));
                })*
                println!();
            }
            $(
                run_tests_for_size::<$size>(&mut planner, repeats);
            )*
            
        }
    };

}

mktest!(f32 f64 i16 i32; 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768);

fn main() {
    println!("|points|   f32   |   f64   |   i16   |   i32   |");
    println!("|-----:|:-------:|:-------:|:-------:|:-------:|");
    // let mut results = vec![Vec::new(); fns.len()];
    #[cfg(debug_assertions)]
    run_tests(1);
    #[cfg(not(debug_assertions))]
    run_tests(1024);
}