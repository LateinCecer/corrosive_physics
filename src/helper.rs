use nalgebra::ComplexField;
use nalgebra::RealField;
use nalgebra::Scalar;
use nalgebra::SimdComplexField;
use nalgebra::SimdRealField;
use num::{One, Zero};
use crate::helper::mat::{Half, Two};

pub mod mat;
pub mod separated_axis;

pub trait BaseFloat : Scalar + ComplexField + RealField + SimdComplexField + SimdRealField
    + Zero + One + Two + Half + Copy
{
    const MIN: Self;
    const MAX: Self;
}

impl BaseFloat for f64 {
    const MIN: Self = f64::MIN;
    const MAX: Self = f64::MAX;
}
impl BaseFloat for f32 {
    const MIN: Self = f32::MIN;
    const MAX: Self = f32::MAX;
}

fn test<T: BaseFloat>() {
    let d = T::simd_sqrt(T::one());
}
