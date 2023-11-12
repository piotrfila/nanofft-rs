#![no_std]
#![cfg_attr(feature = "const", feature(const_mut_refs))]
#![cfg_attr(feature = "const", feature(const_swap))]
mod tables;

use crate::tables::*;

pub type Angle = u32;

#[cfg(feature = "narrow_index_type")]
type Index = u16;
#[cfg(not(feature = "narrow_index_type"))]
type Index = u32;


macro_rules! maybe_const {
    (pub $($x:tt)*) => {
        #[cfg(feature = "const")]
        pub const $($x)*

        #[cfg(not(feature = "const"))]
        pub $($x)*
    };
    ($($x:tt)*) => {
        #[cfg(feature = "const")]
        const $($x)*

        #[cfg(not(feature = "const"))]
        $($x)*
    };
}

macro_rules! generic_fn_variant {
    (
        fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body
        fn $($const_generic)* $(-> $ret)? $body
    };
    (
        const fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        const fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body
        const fn $($const_generic)* $(-> $ret)? $body
    };
    (
        const? fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        maybe_const! { fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body }
        maybe_const! { fn $($const_generic)* $(-> $ret)? $body }
    };
    (
        pub fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        pub fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body
        pub fn $($const_generic)* $(-> $ret)? $body
    };
    (
        pub const fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        pub const fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body
        pub const fn $($const_generic)* $(-> $ret)? $body
    };
    (
        pub const? fn $dynamic:ident$(<$($generic:ident),*>)?($($arg:ident: $arg_type:ty)*) $(-> $ret:ty)? $body:block
        $($const_generic:tt)*
    ) => {
        maybe_const! { pub fn $dynamic$(<$($generic),*>)?($($arg: $arg_type)*) $(-> $ret)? $body }
        maybe_const! { pub fn $($const_generic)* $(-> $ret)? $body }
    };
}


// Use of different types for each array is intentional.
// It allows to also use this function to rearrange a single array
// by passing `&mut [(); N]` as the second argument.
// (and hoping the compiler optimizes out operations on ()s)
// Based on microfft's implementation
generic_fn_variant! { 
    pub const? fn bit_reverse_reorder_dyn<T>(data: &mut [T]) {
        debug_assert!(data.len().is_power_of_two());
        debug_assert!(Index::MAX as usize + 1 >= data.len());

        let shift = Index::BITS - data.len().trailing_zeros();
        let mut it = 0;
        while it < (data.len() as u32) {
            let i = it as Index;
            let rev = i.reverse_bits();
            let j = rev >> shift;
            if j > i {
                data.swap(i as usize, j as usize);
            }
            it += 1;
        }
    }
    bit_reverse_reorder<T, const N: usize>(data: &mut [T; N])
}

generic_fn_variant!{
    pub const? fn interleave_dyn<T>(arr: &mut [(T, T)]) {
        use core::mem::swap;

        let mut log2 = arr.len().trailing_zeros() - 2;
        for i in 0.. {
            for j in 0..(1 << i) {
                let j = j * 4;
                arr[(j + 1) << log2..(j + 3) << log2].rotate_left(1 << log2);
            }
            if log2 == 0 { break; }
            log2 -= 1;
        }
        for j in 0..arr.len() / 2 {
            let j = j * 2;
            let (a, b) = arr.split_at_mut(j + 1);
            swap(&mut a[j].1, &mut b[0].0)
        }
    }
    interleave<T, const N: usize>(arr: &mut [(T, T); N])
}

// Angle represens an angle in range [0, pi)
// other angles are not used in this fft implementation
const fn sin_cos(angle: Angle) -> (TrigTableType, TrigTableType) {
    let shift = Angle::BITS - TRIG_TABLE_BITS;

    let angle_signed = 0_i32.wrapping_add_unsigned(angle);

    if angle_signed < 0 {
        if angle_signed << 1 < 0 {
            let idx = (!(angle << 2) >> shift) as usize;
            (TRIG_TABLE[idx].0, TRIG_TABLE[idx].1)
        }
        else {
            let idx = ((angle << 2) >> shift) as usize;
            (TRIG_TABLE[idx].1, TRIG_TABLE[idx].0)
        }
    }
    else {
        if angle_signed << 1 < 0 {
            let idx = (!(angle << 2) >> shift) as usize;
            (TRIG_TABLE[idx].1, -TRIG_TABLE[idx].0)
        }
        else {
            let idx = ((angle << 2) >> shift) as usize;
            (TRIG_TABLE[idx].0, -TRIG_TABLE[idx].1)
        }
    }
}


macro_rules! fft_impl {
    (
        float; $t:ty; $len:expr;
        ($($arg:ident: $arg_type:ty),*);
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        fn $($signature:tt)*
    ) => {
        fft_impl!(
            $len,
            loop_init: let (mut twiddle_re, mut twiddle_im) = (1., 0.),
            multiply: {
                let product_re = twiddle_re * $y_re - twiddle_im * $y_im;
                let product_im = twiddle_re * $y_im + twiddle_im * $y_re;
                $y_re = $x_re - product_re;
                $y_im = $x_im - product_im;
                $x_re += product_re;
                $x_im += product_im;
            },
            next_twiddle: |angle| {
                let (sin, cos) = crate::sin_cos(angle);
                twiddle_re = cos as $t / (crate::TrigTableType::MAX as $t);
                twiddle_im = sin as $t / (crate::TrigTableType::MAX as $t);
            },
            ($($arg: $arg_type),*);
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            fn $($signature)*
        );
    };

    (
        int; $t:ty; $wide:ty; $len:expr;
        ($($arg:ident: $arg_type:ty),*) $(-> $ret:ident: $ret_type:ty)?;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($signature:tt)*
    ) => {
        fft_impl!(
            $len,
            loop_init: let (mut twiddle_re, mut twiddle_im, scale) = {
                let mut scale = 0;

                let mut $x = 0;
                while $x < $len {
                    let combined = $x_re as $wide | (($x_im as $wide) << (0 as $t).count_zeros());
                    scale |= combined ^ (combined << 1);
                    $x += 1;
                };
                let scale = ((scale >> (1 as $wide).count_zeros()) | (scale >> (1 as $t).count_zeros())) & 1;

                $($ret += (scale as $ret_type))?;

                let (a, b) = (crate::TrigTableType::MAX as $wide, (!(1 as $t).reverse_bits()) as $wide);
                let one = if a < b { a } else { b };
                (one, 0, scale)
            },
            multiply: {
                let (a, b) = ((1 as crate::TrigTableType).count_zeros(), (1 as $t).count_zeros());
                let shift = if a < b { a } else { b };
                if scale != 0 {
                    $y_re >>= 1;
                    $y_im >>= 1;
                    $x_re >>= 1;
                    $x_im >>= 1;
                    let product_re = ((($y_re as $wide * twiddle_re) - ($y_im as $wide * twiddle_im)) >> shift) as $t;
                    let product_im = ((($y_im as $wide * twiddle_re) + ($y_re as $wide * twiddle_im)) >> shift) as $t;
                    $y_re = $x_re - product_re;
                    $y_im = $x_im - product_im;
                    $x_re += product_re;
                    $x_im += product_im;
                }
                else {
                    let product_re = ((($y_re as $wide * twiddle_re) - ($y_im as $wide * twiddle_im)) >> shift) as $t;
                    let product_im = ((($y_im as $wide * twiddle_re) + ($y_re as $wide * twiddle_im)) >> shift) as $t;
                    $y_re = $x_re - product_re;
                    $y_im = $x_im - product_im;
                    $x_re += product_re;
                    $x_im += product_im;
                };
            },
            next_twiddle: |angle| {
                let (sin, cos) = crate::sin_cos(angle);
                let shift = crate::TrigTableType::BITS.saturating_sub((0 as $t).count_zeros());
                twiddle_re = (cos >> shift) as $wide;
                twiddle_im = (sin >> shift) as $wide;
            },
            ($($arg: $arg_type),*) $(-> $ret: $ret_type)?;
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            $($signature)*
        );
    };

    (
        $len:expr,
        loop_init: $loop_init:stmt,
        multiply: $mul:block,
        next_twiddle: |$angle:ident| $next_twiddle:block,
        ($($arg:ident: $arg_type:ty),*) $(-> $ret:ident: $ret_type:ty)?;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($signature:tt)*
    ) => {

    $($signature)* ($($arg: $arg_type),* $(, mut $ret: $ret_type)?) $(-> $ret_type)? {
        let mut step_log2 = 0;
        let mut step = 1;
        while {
            step < $len
        } {
            let jump = step << 1;
            $loop_init
            let mut group = 0;
            loop {
                let mut $x = group;
                while $x < $len {
                    let $y = $x + step;
                    $mul
                    $x += jump;
                }
                
                // we need the factors below for the next iteration
                // if we don't iterate then don't compute
                group += 1;
                if group == step { break }

                let $angle = (group as crate::Angle) << (crate::Angle::BITS - step_log2);

                $next_twiddle
            }
            step_log2 += 1;
            step = 1 << step_log2;
        }
        $($ret)?
    }

    };
}

macro_rules! type_impl {
    ($kind:tt; $div2:expr; $($ret:ident = $ret_init:literal: $ret_type:ty)?; $mod:ident, $t:ty, $($wide:ty)? $(,$qualifier:tt)?) => { pub mod $mod {
    fft_impl!(
        $kind; $t; $($wide;)? N;
        (data: &mut [($t, $t); N]) $(-> $ret: $ret_type)?;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
        $($qualifier)? fn compute_pairs<const N: usize>
    );

    fft_impl!(
        $kind; $t; $($wide;)? data.len();
        (data: &mut [($t, $t)]) $(-> $ret: $ret_type)?;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
        $($qualifier)? fn compute_pairs_dyn
    );

    fft_impl!(
        $kind; $t; $($wide;)? N;
        (re: &mut [$t; N], im: &mut [$t; N]) $(-> $ret: $ret_type)?;
        a; re[a]; im[a];
        b; re[b]; im[b];
        $($qualifier)? fn compute_arrays<const N: usize>
    );

    generic_fn_variant!{
        pub const? fn rfft_twice_postprocess_dyn(arr: &mut [($t, $t)]) {
            use core::mem::swap;

            let mid = arr.len() / 2;
            let (a, b) = arr.split_at_mut(mid);
            swap(&mut a[0].1, &mut b[0].0);
            for (x, y) in a[1..].iter_mut().zip(b.iter_mut().rev()) {
                let z_re = $div2(x.0);
                let z_im = $div2(x.1);
                let w_re = $div2(y.0);
                let w_im = $div2(y.1);
                x.0 = z_re + w_re;
                x.1 = z_im - w_im;
                y.0 = w_im + z_im;
                y.1 = w_re - z_re;
            }

            for i in 1..b.len() / 2 {
                b.swap(i, b.len() - i);
            }
        }
        rfft_twice_postprocess<const N: usize>(arr: &mut [($t, $t); N])
    }

    pub $($qualifier)? fn rfft_pairs_twice<const N: usize>(data: &mut [($t, $t); N]) $(-> $ret_type)? {
        debug_assert!(data.len().is_power_of_two());
        super::interleave(data);
        let ret = fft_pairs(data);
        rfft_twice_postprocess(data);
        ret
    }

    pub $($qualifier)? fn rfft_pairs_twice_dyn(data: &mut [($t, $t)]) $(-> $ret_type)? {
        debug_assert!(data.len().is_power_of_two());
        super::interleave_dyn(data);
        let ret = fft_pairs_dyn(data);
        rfft_twice_postprocess_dyn(data);
        ret
    }

    pub $($qualifier)? fn fft_pairs<const N: usize>(data: &mut [($t, $t); N]) $(-> $ret_type)? {
        debug_assert!(data.len().is_power_of_two());
        super::bit_reverse_reorder(data);
        compute_pairs(data $(, $ret_init)?)
    }

    pub $($qualifier)? fn fft_pairs_dyn(data: &mut [($t, $t)]) $(-> $ret_type)? {
        debug_assert!(data.len().is_power_of_two());
        super::bit_reverse_reorder_dyn(data);
        compute_pairs_dyn(data $(, $ret_init)?)
    }

    pub $($qualifier)? fn fft_arrays<const N: usize>(data_re: &mut [$t; N], data_im: &mut [$t; N]) $(-> $ret_type)? {
        debug_assert!(N.is_power_of_two());
        super::bit_reverse_reorder(data_re);
        super::bit_reverse_reorder(data_im);
        compute_arrays(data_re, data_im $(, $ret_init)?)
    }

    } };
}

type_impl!(float; |x| x * 0.5_f32;; f32, f32,);
type_impl!(float; |x| x * 0.5_f64;; f64, f64,);
#[cfg(feature = "const")]
type_impl!(int; |x| x >> 1; lsb_mult_log2 = -15: i16; i16, i16, i32, const);
#[cfg(feature = "const")]
type_impl!(int; |x| x >> 1; lsb_mult_log2 = -31: i16; i32, i32, i64, const);
#[cfg(not(feature = "const"))]
type_impl!(int; |x| x >> 1; lsb_mult_log2 = -15: i16; i16, i16, i32);
#[cfg(not(feature = "const"))]
type_impl!(int; |x| x >> 1; lsb_mult_log2 = -31: i16; i32, i32, i64);
