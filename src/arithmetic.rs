
use core::ops::{ Add, AddAssign, Sub };
use crate::tables::*;

pub trait Arithmetic where Self: Sized + Copy + Add<Output = Self> + Sub<Output = Self> + AddAssign {
    type ScaleInfo: Copy;
    type RangeInfo: Copy;

    fn from_f64(x: f64) -> Self;
    fn into_f64(self, s: Self::RangeInfo) -> f64;

    fn one() -> Self { Self::from_f64(1.0) }
    fn zero() -> Self { Self::from_f64(0.0) }

    fn scale_init() -> Self::ScaleInfo;
    fn scale_update(self, _: &mut Self::ScaleInfo) { }
    fn range_init() -> Self::RangeInfo;
    fn range_update(_: &mut Self::ScaleInfo, _: &mut Self::RangeInfo) { }
    fn mul(a: Self, b: Self, scale: Self::ScaleInfo) -> Self;

    fn sin_cos(angle: f32) -> (Self, Self);
}

macro_rules! arithmetic_float {
    ($($t:ty)*) => { $(

    impl Arithmetic for $t {
        type ScaleInfo = ();
        type RangeInfo = ();
    
        fn from_f64(x: f64) -> Self { x as Self }
        fn into_f64(self, _: Self::RangeInfo) -> f64 { self as f64 }
    
        fn scale_init() -> Self::ScaleInfo { () }
        fn range_init() -> Self::RangeInfo { () }
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
        }
    }
    
    )* };
}

macro_rules! arithmetic_int {
    ($(($t:ty, $wide:ty))*) => { $(

    impl Arithmetic for $t {
        type ScaleInfo = $t;
        type RangeInfo = i16;
    
        fn from_f64(x: f64) -> Self { (Self::MAX as f64 * x) as Self }
        fn into_f64(self, i: Self::RangeInfo) -> f64 { 
            let exponent = (i as i32 - f64::MIN_EXP + 2) as u64 & ((1 << (u64::BITS - f64::MANTISSA_DIGITS - 1)) - 1);
            let scale = f64::from_bits(exponent << f64::MANTISSA_DIGITS);
            scale * self as f64
        }
    
        fn scale_init() -> Self::ScaleInfo { 0 }
        fn scale_update(self, i: &mut Self::ScaleInfo) { *i |= self.abs() }
        fn range_init() -> Self::RangeInfo { 1 - (Self::BITS as Self::RangeInfo) }
        fn range_update(s: &mut Self::ScaleInfo, r: &mut Self::RangeInfo) {
            let margin = s.leading_zeros();
            *s = (1 + Self::BITS - margin) as _;
            *r += 2 - (margin as Self::RangeInfo);
        }
    
        fn mul(a: Self, b: Self, s: Self::ScaleInfo) -> Self { ((a as $wide * (b as $wide)) >> (s as u32)) as Self }
    
        fn sin_cos(angle: f32) -> (Self, Self) {
            let len = TRIG_TABLE.len() - 1;
            let (a, b) = if angle < 0.5 {
                let idx = (len as f32 * angle * 2.) as usize;
                (TRIG_TABLE[idx], -TRIG_TABLE[len - idx])
            } else {
                let idx = (len as f32 * (angle - 0.5) * 2.) as usize;
                (TRIG_TABLE[len - idx], TRIG_TABLE[idx])
            };
            #[cfg(not(feature = "wide_trig_lut"))]
            let (a, b) = (
                (a as Self) << (Self::BITS - TrigTableT::BITS),
                (b as Self) << (Self::BITS - TrigTableT::BITS),
            );
            #[cfg(feature = "wide_trig_lut")]
            let (a, b) = (
                (a >> (TrigTableT::BITS - Self::BITS)) as Self,
                (b >> (TrigTableT::BITS - Self::BITS)) as Self,
            );
            (a, b)
        }
    }
    
    )* };
}

arithmetic_float!(f32 f64);
arithmetic_int!((i16, i32) (i32, i64));