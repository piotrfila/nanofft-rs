#![no_std]

mod arithmetic;
mod tables;

pub use arithmetic::Arithmetic;

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
fn compute_arrays<T: Arithmetic, const N: usize>(re: &mut [T; N], im: &mut [T; N]) -> T::RangeInfo {
    let mut range = T::range_init();
    let mut step = 1;
    while step < N {
        let jump = step << 1;
        let step_d = step as f32;
        let mut twiddle_re = T::one();
        let mut twiddle_im = T::zero();
        let mut scale = T::scale_init();
        for re in re.iter() {
            re.scale_update(&mut scale);
        }
        for im in im.iter() {
            im.scale_update(&mut scale);
        }
        T::range_update(&mut scale, &mut range);
        for group in 0..step {
            let mut pair = group;
            while pair < N {
                let matchh = pair + step;
                let product_re = T::mul(twiddle_re, re[matchh], scale) - T::mul(twiddle_im, im[matchh], scale);
                let product_im = T::mul(twiddle_re, im[matchh], scale) + T::mul(twiddle_im, re[matchh], scale);
                re[matchh] = T::mul(re[pair], T::from_f64(1.), scale) - product_re;
                im[matchh] = T::mul(im[pair], T::from_f64(1.), scale) - product_im;
                re[pair] = T::mul(re[pair], T::from_f64(1.), scale) + product_re;
                im[pair] = T::mul(im[pair], T::from_f64(1.), scale) + product_im;
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
    range
}

pub fn fft_arrays<T: Arithmetic, const N: usize>(data_re: &mut [T; N], data_im: &mut [T; N]) -> T::RangeInfo {
    bit_reverse_reorder(data_re, data_im);
    compute_arrays(data_re, data_im)
}