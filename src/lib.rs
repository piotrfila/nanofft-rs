mod tables;

use core::ops::{ Add, AddAssign, Sub };
use tables::*;

pub trait Arithmetic where Self: Sized + Copy + Add<Output = Self> + Sub<Output = Self> + AddAssign {
    type ScaleInfo: Copy;

    fn from_f64(x: f64) -> Self;
    fn into_f64(self) -> f64;

    fn one() -> Self { Self::from_f64(1.0) }
    fn zero() -> Self { Self::from_f64(0.0) }

    fn scale_init() -> Self::ScaleInfo;
    fn mul(a: Self, b: Self, scale: Self::ScaleInfo) -> Self;

    fn sin_cos(angle: f32) -> (Self, Self);
}

impl Arithmetic for f32 {
    type ScaleInfo = ();

    fn from_f64(x: f64) -> Self { x as Self }
    fn into_f64(self) -> f64 { self as f64 }

    fn scale_init() -> Self::ScaleInfo { () }
    fn mul(a: Self, b: Self, _: Self::ScaleInfo) -> Self { a * b }

    fn sin_cos(angle: f32) -> (Self, Self) {
        let len = TRIG_TABLE.len() - 1;
        let (a, b) = if angle < 0.5 {
            let idx = (len as f32 * angle * 2.) as usize;
            (TRIG_TABLE[idx], -TRIG_TABLE[len - idx])
        } else {
            let idx = (len as f32 * (angle - 0.5) * 2.) as usize;
            (TRIG_TABLE[len - idx], TRIG_TABLE[idx])
        };
        (
            a as Self / (TrigTableT::MAX as Self),
            b as Self / (TrigTableT::MAX as Self),
        )
        // ((-angle * std::f32::consts::PI) as Self).sin_cos()
    }
}


impl Arithmetic for f64 {
    type ScaleInfo = ();

    fn from_f64(x: f64) -> Self { x as Self}
    fn into_f64(self) -> f64 { self as f64 }

    fn scale_init() -> Self::ScaleInfo { () }
    fn mul(a: Self, b: Self, _: Self::ScaleInfo) -> Self { a * b }

    fn sin_cos(angle: f32) -> (Self, Self) {
        let len = TRIG_TABLE.len() - 1;
        let (a, b) = if angle < 0.5 {
            let idx = (len as f32 * angle * 2.) as usize;
            (TRIG_TABLE[idx], -TRIG_TABLE[len - idx])
        } else {
            let idx = (len as f32 * (angle - 0.5) * 2.) as usize;
            (TRIG_TABLE[len - idx], TRIG_TABLE[idx])
        };
        (
            a as Self / (TrigTableT::MAX as Self),
            b as Self / (TrigTableT::MAX as Self),
        )
        // ((-angle * std::f32::consts::PI) as Self).sin_cos()
    }
}

// Use of different types for each array is intentional.
// It allows to also use this function to rearrange a single array
// by passing `&mut [(); N]` as the second argument.
// (and hoping the compiler optimizes out operations on ()s)
// Based on microfft's implementation
fn bit_reverse_reorder<A, B, const N: usize>(re: &mut [A; N], im: &mut [B; N]) {
    debug_assert!(N.is_power_of_two());

    let shift = core::mem::size_of::<usize>() as u32 * 8 - N.trailing_zeros();
    for i in 0..N {
        let rev = i.reverse_bits();
        let j = rev >> shift;
        if j > i {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
}

// Based on [this implementation](https://lloydrochester.com/post/c/example-fft/#test-cases-for-the-fft)
fn compute_arrays<T: Arithmetic, const N: usize>(re: &mut [T; N], im: &mut [T; N]) {
    let scale = T::scale_init();
    let mut step = 1;
    while step < N {
        let jump = step << 1;
        let step_d = step as f32;
        let mut twiddle_re = T::one();
        let mut twiddle_im = T::zero();
        for group in 0..step {
            let mut pair = group;
            while pair < N {
                let matchh = pair + step;
                let product_re = T::mul(twiddle_re, re[matchh], scale) - T::mul(twiddle_im, im[matchh], scale);
                let product_im = T::mul(twiddle_re, im[matchh], scale) + T::mul(twiddle_im, re[matchh], scale);
                re[matchh] = re[pair] - product_re;
                im[matchh] = im[pair] - product_im;
                re[pair] += product_re;
                im[pair] += product_im;

                pair+=jump
            }
            
            // we need the factors below for the next iteration
            // if we don't iterate then don't compute
            if group+1 == step { continue }

            let angle = (group as f32 + 1.) / step_d;
            (twiddle_im, twiddle_re) = T::sin_cos(angle);
        }
        step <<= 1;
    }
}

pub fn fft_arrays<T: Arithmetic, const N: usize>(data_re: &mut [T; N], data_im: &mut [T; N]) {
    bit_reverse_reorder(data_re, data_im);
    compute_arrays(data_re, data_im);
}