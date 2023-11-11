#![no_std]
mod tables;

use crate::tables::*;

pub type Angle = u32;

#[cfg(feature = "narrow_indyex_type")]
type Index = u16;
#[cfg(not(feature = "narrow_index_type"))]
type Index = u32;


// Use of different types for each array is intentional.
// It allows to also use this function to rearrange a single array
// by passing `&mut [(); N]` as the second argument.
// (and hoping the compiler optimizes out operations on ()s)
// Based on microfft's implementation
fn bit_reverse_reorder<A, B, const N: usize>(re: &mut [A; N], im: &mut [B; N]) {
    debug_assert!(N.is_power_of_two());
    debug_assert!(Index::MAX as usize + 1 >= N);

    let shift = Index::BITS - N.trailing_zeros();
    for i in 0..(N as u32) {
        let i = i as Index;
        let rev = i.reverse_bits();
        let j = rev >> shift;
        if j > i {
            re.swap(i as usize, j as usize);
            im.swap(i as usize, j as usize);
        }
    }
}

fn bit_reverse_reorder_dyn<T>(data: &mut [T]) {
    debug_assert!(data.len().is_power_of_two());
    debug_assert!(Index::MAX as usize + 1 >= data.len());

    let shift = Index::BITS - data.len().trailing_zeros();
    for i in 0..(data.len() as u32) {
        let i = i as Index;
        let rev = i.reverse_bits();
        let j = rev >> shift;
        if j > i {
            data.swap(i as usize, j as usize);
        }
    }
}


// Angle represens an angle in range [0, pi)
// other angles are not used in this fft implementation
fn sin_cos(angle: Angle) -> (TrigTableType, TrigTableType) {
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
        $fname:ident($($arg:ident: $arg_type:ty),*);
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
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
            $fname($($arg: $arg_type),*);
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            $($generics)*
        );
    };

    (
        int; $t:ty; $wide:ty; $len:expr;
        $fname:ident($($arg:ident: $arg_type:ty),*) $(-> $ret:ident: $ret_type:ty)?;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
    ) => {
        fft_impl!(
            $len,
            loop_init: let (mut twiddle_re, mut twiddle_im, scale) = {
                let mut scale = 0;

                for $x in 0..$len {
                    let combined = $x_re as $wide | (($x_im as $wide) << (0 as $t).count_zeros());
                    scale |= combined ^ (combined << 1);
                };
                let scale = ((scale >> (1 as $wide).count_zeros()) | (scale >> (1 as $t).count_zeros())) & 1;

                $($ret += (scale as $ret_type))?;

                let one = (crate::TrigTableType::MAX as $wide).min((!(1 as $t).reverse_bits()) as $wide);
                (one, 0, scale)
            },
            multiply: {
                let shift = (1 as crate::TrigTableType).count_zeros().min((1 as $t).count_zeros());
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
            $fname($($arg: $arg_type),*) $(-> $ret: $ret_type)?;
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            $($generics)*
        );
    };

    (
        $len:expr,
        loop_init: $loop_init:stmt,
        multiply: $mul:block,
        next_twiddle: |$angle:ident| $next_twiddle:block,
        $fname:ident($($arg:ident: $arg_type:ty),*) $(-> $ret:ident: $ret_type:ty)?;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
    ) => {

    fn $fname $($generics)* ($($arg: $arg_type),* $(, mut $ret: $ret_type)?) $(-> $ret_type)? {
        let mut step_log2 = 0;
        let mut step = 1;
        while {
            step < $len
        } {
            let jump = step << 1;
            $loop_init
            for group in 0..step {
                let mut $x = group;
                while $x < $len {
                    let $y = $x + step;
                    $mul
                    $x += jump;
                }
                
                // we need the factors below for the next iteration
                // if we don't iterate then don't compute
                if group + 1 == step { continue }

                let $angle = (group as crate::Angle + 1) << (crate::Angle::BITS - step_log2);

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
    ($kind:tt; $($ret:ident = $ret_init:literal: $ret_type:ty)?; $mod:ident, $t:ty $(, $wide:ty)?) => { pub mod $mod {
    fft_impl!(
        $kind; $t; $($wide;)? N;
        compute_pairs(data: &mut [($t, $t); N]) $(-> $ret: $ret_type)?;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
        <const N: usize>
    );

    fft_impl!(
        $kind; $t; $($wide;)? data.len();
        compute_pairs_dyn(data: &mut [($t, $t)]) $(-> $ret: $ret_type)?;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
    );

    fft_impl!(
        $kind; $t; $($wide;)? N;
        compute_arrays(re: &mut [$t; N], im: &mut [$t; N]) $(-> $ret: $ret_type)?;
        a; re[a]; im[a];
        b; re[b]; im[b];
        <const N: usize>
    );

    pub fn fft_pairs<const N: usize>(data: &mut [($t, $t); N]) $(-> $ret_type)? {
        super::bit_reverse_reorder(data, &mut [(); N]);
        compute_pairs(data $(, $ret_init)?)
    }

    pub fn fft_pairs_dyn(data: &mut [($t, $t)]) $(-> $ret_type)? {
        super::bit_reverse_reorder_dyn(data);
        compute_pairs_dyn(data $(, $ret_init)?)
    }

    pub fn fft_arrays<const N: usize>(data_re: &mut [$t; N], data_im: &mut [$t; N]) $(-> $ret_type)? {
        super::bit_reverse_reorder(data_re, data_im);
        compute_arrays(data_re, data_im $(, $ret_init)?)
    }

    } };
}

type_impl!(float;; f32, f32);
type_impl!(float;; f64, f64);
type_impl!(int; lsb_mult_log2 = -15: i16; i16, i16, i32);
type_impl!(int; lsb_mult_log2 = -31: i16; i32, i32, i64);
