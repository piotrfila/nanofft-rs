#![no_std]
mod arithmetic;
mod tables;

pub use arithmetic::Arithmetic;

#[cfg(feature = "narrow_index_type")]
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
    debug_assert!(Index::MAX as usize >= N);

    let shift = Index::BITS - N.trailing_zeros();
    for i in 0..(N as Index) {
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
    debug_assert!(Index::MAX as usize >= data.len());

    let shift = Index::BITS - data.len().trailing_zeros();
    for i in 0..(data.len() as Index) {
        let rev = i.reverse_bits();
        let j = rev >> shift;
        if j > i {
            data.swap(i as usize, j as usize);
        }
    }
}

macro_rules! fft_impl {
    (
        float; $t:ty; $len:expr;
        $fname:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
    ) => {
        fft_impl!(
            0. as $t; 1. as $t; (); (); $len;
            |x, y, _| x * y;
            |ang| <$t as Arithmetic>::sin_cos(ang);
            |_, _| ();
            |_, _| ();
            |x, _| x;
            $fname($($arg: $arg_type),*) -> ();
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            $($generics)*
        );
    };

    (
        int; $t:ty; $wide:ty; $len:expr;
        $fname:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
    ) => {
        fft_impl!(
            0 as $t; (!(1 as $t).reverse_bits()); 0 as $t; -((1 as $t).count_zeros() as i16); $len;
            |x, y, s| ((x as $wide * (y as $wide)) >> (s as u32)) as $t;
            |ang| <$t as Arithmetic>::sin_cos(ang);
            |x: $t, s| s | x.abs();
            |s: &mut $t, r: &mut i16| {
                let margin = s.leading_zeros();
                *s = (1 + (0 as $t).count_zeros() - margin) as _;
                *r += 2 - (margin as i16);
            };
            |x, s| x >> ((s as u32 + 1) - (0 as $t).count_zeros());
            $fname($($arg: $arg_type),*) -> i16;
            $x; $x_re; $x_im;
            $y; $y_re; $y_im;
            $($generics)*
        );
    };

    (
        $zero:expr; $one:expr; $scale_init:expr; $range_init:expr; $len:expr;
        $mul:expr;
        $sin_cos:expr;
        $scale_update:expr;
        $range_update:expr;
        $scale:expr;
        $fname:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty;
        $x:ident; $x_re:expr; $x_im:expr;
        $y:ident; $y_re:expr; $y_im:expr;
        $($generics:tt)*
    ) => {

    fn $fname $($generics)* ($($arg: $arg_type),*) -> $ret {
        let mut range = $range_init;
        let mut step = 1;
        while step < $len {
            let jump = step << 1;
            let mut twiddle_re = $one;
            let mut twiddle_im = $zero;
            let mut scale = $scale_init;
            for $x in 0..$len {
                scale = $scale_update($x_re, scale);
                scale = $scale_update($x_im, scale);
            }
            $range_update(&mut scale, &mut range);
            for group in 0..step {
                let mut $x = group;
                while $x < $len {
                    let $y = $x + step;
                    let product_re = $mul(twiddle_re, $y_re, scale) - $mul(twiddle_im, $y_im, scale);
                    let product_im = $mul(twiddle_re, $y_im, scale) + $mul(twiddle_im, $y_re, scale);
                    $y_re = $scale($x_re, scale) - product_re;
                    $y_im = $scale($x_im, scale) - product_im;
                    $x_re = $scale($x_re, scale) + product_re;
                    $x_im = $scale($x_im, scale) + product_im;
                    $x += jump
                }
                
                // we need the factors below for the next iteration
                // if we don't iterate then don't compute
                if group + 1 == step { continue }
    
                let angle = (group + 1) * (65536 / step);
                (twiddle_im, twiddle_re) = $sin_cos(angle as _);
            }
            step <<= 1;
        }
        range
    }

    };
}

macro_rules! type_impl {
    ($kind:tt; $($mod:ident, $t:ty $(, $wide:ty)?);*) => { $( pub mod $mod {

    use super::Arithmetic;

    fft_impl!(
        $kind; $t; $($wide;)? N;
        compute_pairs(data: &mut [($t, $t); N]) -> <$t as Arithmetic>::RangeInfo;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
        <const N: usize>
    );

    fft_impl!(
        $kind; $t; $($wide;)? data.len();
        compute_pairs_dyn(data: &mut [($t, $t)]) -> <$t as Arithmetic>::RangeInfo;
        a; data[a].0; data[a].1;
        b; data[b].0; data[b].1;
    );

    fft_impl!(
        $kind; $t; $($wide;)? N;
        compute_arrays(re: &mut [$t; N], im: &mut [$t; N]) -> <$t as Arithmetic>::RangeInfo;
        a; re[a]; im[a];
        b; re[b]; im[b];
        <const N: usize>
    );

    pub fn fft_pairs<const N: usize>(data: &mut [($t, $t); N]) -> <$t as Arithmetic>::RangeInfo {
        super::bit_reverse_reorder(data, &mut [(); N]);
        compute_pairs(data)
    }

    pub fn fft_pairs_dyn(data: &mut [($t, $t)]) -> <$t as Arithmetic>::RangeInfo {
        super::bit_reverse_reorder_dyn(data);
        compute_pairs_dyn(data)
    }

    pub fn fft_arrays<const N: usize>(data_re: &mut [$t; N], data_im: &mut [$t; N]) -> <$t as Arithmetic>::RangeInfo {
        super::bit_reverse_reorder(data_re, data_im);
        compute_arrays(data_re, data_im)
    }

    } )* };
}

type_impl!(float; f32, f32; f64, f64);
type_impl!(int; i16, i16, i32; i32, i32, i64);
