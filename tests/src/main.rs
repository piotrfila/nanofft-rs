extern crate rustfft;
extern crate rand;
extern crate nanofft;

use rustfft::{ FftPlanner, num_complex::Complex };
use rand::{ Rng, thread_rng };

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
    ($t:ty, arrays, $f:expr) => {
        |data: &mut [Complex<f64>; N]| {
            let mut data_re = data.map(|x| <$t as Convert>::from_f64(x.re));
            let mut data_im = data.map(|x| <$t as Convert>::from_f64(x.im));
            let range = $f(&mut data_re, &mut data_im);
            for (dst, src) in data.iter_mut().zip(data_re.iter().zip(data_im.iter())) {
                dst.re = src.0.into_f64(range);
                dst.im = src.1.into_f64(range);
            }
        }
    };

    ($t:ty, pairs, $f:expr) => {
        |data: &mut [Complex<f64>; N]| {
            let mut data_t = data.map(|x| (<$t as Convert>::from_f64(x.re), <$t as Convert>::from_f64(x.im)));
            let range = $f(&mut data_t);
            for (dst, src) in data.iter_mut().zip(data_t.iter()) {
                dst.re = src.0.into_f64(range);
                dst.im = src.1.into_f64(range);
            }
        }
    };
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

    #[cfg(debug_assertions)]
    let fns = [
        mktest!(f32, pairs, nanofft::f32::fft_pairs_dyn),
        mktest!(f64, pairs, nanofft::f64::fft_pairs_dyn),
        mktest!(i16, pairs, nanofft::i16::fft_pairs_dyn),
        mktest!(i32, pairs, nanofft::i32::fft_pairs_dyn),
    ];
    #[cfg(not(debug_assertions))]
    let fns = [
        mktest!(f32, pairs, nanofft::f32::fft_pairs_dyn),
        mktest!(f64, pairs, nanofft::f64::fft_pairs_dyn),
        mktest!(i16, pairs, nanofft::i16::fft_pairs_dyn),
        mktest!(i32, pairs, nanofft::i32::fft_pairs_dyn),
    ];
    
    fns.into_iter()
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
    })
    .collect::<Vec<_>>()
}

fn main() {
    let fns = [
        test_size::<    4>,
        test_size::<    8>,
        test_size::<   16>,
        test_size::<   32>,
        test_size::<   64>,
        test_size::<  128>,
        test_size::<  256>,
        test_size::<  512>,
        test_size::< 1024>,
        test_size::< 2048>,
        test_size::< 4096>,
        test_size::< 8192>,
        test_size::<16384>,
        test_size::<32768>,
    ];

    let mut planner = FftPlanner::new();
    println!("|points|   f32   |   f64   |   i16   |   i32   |");
    println!("|-----:|:-------:|:-------:|:-------:|:-------:|");
    let mut results = vec![Vec::new(); fns.len()];
    #[cfg(debug_assertions)]
    let repeats = 1;
    #[cfg(not(debug_assertions))]
    let repeats = 1024;
    for _ in 0..repeats {
        for (r, f) in results.iter_mut().zip(fns.iter()) {
            let ret = f(&mut planner);
            for (i, e) in ret.iter().enumerate() {
                if r.len() <= i { r.push(0.) }
                r[i] += e;
            }
        }
    }
    for (i, errors) in results.iter().enumerate() {
        print!("|{:>6}|", 2_usize.pow(i as u32 + 2));
        for e in errors {
            print!("{:9.3e}|", e / (repeats as f64));
        }
        println!();
    }
}