use nalgebra::ClosedAdd;
use nalgebra::Matrix4;
use nalgebra::UnitQuaternion;
use nalgebra::Vector3;
use num::One;
use crate::helper::BaseFloat;

pub trait Two {
    /// Returns the additional double value of the `one` value.
    fn two() -> Self;
}

impl<T> Two for T
where T: One + ClosedAdd<Output=T> {
    fn two() -> Self {
        T::one() + T::one()
    }
}

pub trait Half {
    /// returns 0.5 for the number type
    fn half() -> Self;
}

impl Half for f32 {
    fn half() -> Self {
        0.5
    }
}

impl Half for f64 {
    fn half() -> Self {
        0.5
    }
}


/// Index of quaternion entry `i`
const I: usize = 0;
/// Index of quaternion entry `j`
const J: usize = 1;
/// Index of quaternion entry `k`
const K: usize = 2;
/// Index of quaternion entry `w`
const W: usize = 3;


/// Generates a 4x4 translation matrix from a given 3d point.
pub fn init_translation<T>(pos: &Vector3<T>) -> Matrix4<T>
where T: BaseFloat {
    Matrix4::new(
        T::one(),  T::zero(), T::zero(), pos[0],
        T::zero(), T::one(),  T::zero(), pos[1],
        T::zero(), T::zero(), T::one(),  pos[2],
        T::zero(), T::zero(), T::zero(), T::one(),
    )
}

/// Generates an inverse 4x4 translation matrix from a given 3d point.
pub fn init_inverse_translation<T>(pos: &Vector3<T>) -> Matrix4<T>
where T: BaseFloat {
    Matrix4::new(
        T::one(),  T::zero(), T::zero(), -pos[0],
        T::zero(), T::one(),  T::zero(), -pos[1],
        T::zero(), T::zero(), T::one(),  -pos[2],
        T::zero(), T::zero(), T::zero(), T::one(),
    )
}

/// Generates a 4x4 rotation matrix from a given 3d rotation unit quaternion.
pub fn init_rotation<T>(rot: &UnitQuaternion<T>) -> Matrix4<T>
where T: BaseFloat {

    Matrix4::new(
        // right
        T::one() - T::two() * (rot.coords[J] * rot.coords[J] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[I] * rot.coords[J] - rot.coords[W] * rot.coords[K]),
        T::two() * (rot.coords[W] * rot.coords[J] + rot.coords[I] * rot.coords[K]),
        T::zero(),
        // up
        T::two() * (rot.coords[W] * rot.coords[K] + rot.coords[I] * rot.coords[J]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[J] * rot.coords[K] - rot.coords[W] * rot.coords[I]),
        T::zero(),
        // forward
        T::two() * (rot.coords[I] * rot.coords[K] - rot.coords[W] * rot.coords[J]),
        T::two() * (rot.coords[W] * rot.coords[I] + rot.coords[J] * rot.coords[K]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[J] * rot.coords[J]),
        T::zero(),
        //
        T::zero(), T::zero(), T::zero(), T::one(),
    )
}

/// Generates an inverse 4x4 rotation matrix from a given 3d rotation unit quaternion.
pub fn init_inverse_rotation<T>(rot: &UnitQuaternion<T>) -> Matrix4<T>
where T: BaseFloat {

    Matrix4::new(
        // right
        T::one() - T::two() * (rot.coords[J] * rot.coords[J] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[I] * rot.coords[J] + rot.coords[W] * rot.coords[K]),
        T::two() * (rot.coords[I] * rot.coords[K] - rot.coords[W] * rot.coords[J]),
        T::zero(),
        // up
        T::two() * (rot.coords[I] * rot.coords[J] - rot.coords[W] * rot.coords[K]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[J] * rot.coords[K] + rot.coords[W] * rot.coords[I]),
        T::zero(),
        // forward
        T::two() * (rot.coords[I] * rot.coords[K] + rot.coords[W] * rot.coords[J]),
        T::two() * (rot.coords[J] * rot.coords[K] - rot.coords[W] * rot.coords[I]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[J] * rot.coords[J]),
        T::zero(),
        //
        T::zero(), T::zero(), T::zero(), T::one(),
    )
}

/// Generates a 4x4 scale matrix from a 3d scale vector.
pub fn init_scale<T>(scale: &Vector3<T>) -> Matrix4<T>
where T: BaseFloat {
    Matrix4::new(
        scale[0], T::zero(), T::zero(), T::zero(),
        T::zero(), scale[1], T::zero(), T::zero(),
        T::zero(), T::zero(), scale[2], T::zero(),
        T::zero(), T::zero(), T::zero(), T::one()
    )
}

/// Generates an inverse 4x4 scale matrix from a 3d scale vector.
pub fn init_inverse_scale<T>(scale: &Vector3<T>) -> Matrix4<T>
where T: BaseFloat {
    Matrix4::new(
        T::one() / scale[0], T::zero(), T::zero(), T::zero(),
        T::zero(), T::one() / scale[1], T::zero(), T::zero(),
        T::zero(), T::zero(), T::one() / scale[2], T::zero(),
        T::zero(), T::zero(), T::zero(), T::one(),
    )
}

/// Returns the vector pointing in the right direction from the rotation unit-quaternion.
pub fn right<T>(rot: &UnitQuaternion<T>) -> Vector3<T>
where T: BaseFloat {

    Vector3::new(
        T::one() - T::two() * (rot.coords[J] * rot.coords[J] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[W] * rot.coords[K] + rot.coords[I] * rot.coords[J]),
        T::two() * (rot.coords[I] * rot.coords[K] - rot.coords[W] * rot.coords[J]),
    )
}

/// Returns the vector pointing in the upward direction from the rotation unit-quaternion.
pub fn up<T>(rot: &UnitQuaternion<T>) -> Vector3<T>
where T: BaseFloat {

    Vector3::new(
        T::two() * (rot.coords[I] * rot.coords[J] - rot.coords[W] * rot.coords[K]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[K] * rot.coords[K]),
        T::two() * (rot.coords[W] * rot.coords[I] + rot.coords[J] * rot.coords[K]),
    )
}

/// Returns the vector pointing in the forward direction from the rotation unit-quaternion.
pub fn forward<T>(rot: &UnitQuaternion<T>) -> Vector3<T>
where T: BaseFloat {

    Vector3::new(
        T::two() * (rot.coords[W] * rot.coords[J] + rot.coords[I] * rot.coords[K]),
        T::two() * (rot.coords[J] * rot.coords[K] - rot.coords[W] * rot.coords[I]),
        T::one() - T::two() * (rot.coords[I] * rot.coords[I] + rot.coords[J] * rot.coords[J])
    )
}
