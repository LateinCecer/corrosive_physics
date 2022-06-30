use nalgebra::{SVector, Vector2, Vector3};
use num::Signed;
use crate::helper::BaseFloat;


macro_rules! abs {
    ($T:ty, $val:expr) => (
        <$T as Signed>::abs(&$val)
    )
}

macro_rules! add {
    ($x:expr) => ($x);
    ($x:expr, $($y:expr),+) => (
        $x + add!($($y),+)
    )
}

macro_rules! intersect_projection {
    ($T:ty, $axis:expr, $($fac:ident * |$p:ident|),+ $(+ $add:expr)*) => (
        abs!($T, $axis) > add!($($add, )* $($fac * abs!($T, $p)),+)
    );
}


/// OBB-OBB non-intersection test using the separating axis theorem in three spacial dimensions.
pub fn intersects_obb_obb<T: BaseFloat>(
    a0: &Vector3<T>, a1: &Vector3<T>, a2: &Vector3<T>,
    b0: &Vector3<T>, b1: &Vector3<T>, b2: &Vector3<T>,
    rel: &Vector3<T>,
    sa0: T, sa1: T, sa2: T,
    sb0: T, sb1: T, sb2: T
) -> bool {


    // -- axis A0
    let c00 = a0.dot(b0);
    let c01 = a0.dot(b1);
    let c02 = a0.dot(b2);

    let a0d = a0.dot(rel);
    if intersect_projection!(T, a0d, sb0 * |c00|, sb1 * |c01|, sb2 * |c02| + sa0) {
        return false;
    }

    // -- axis A1
    let c10 = a1.dot(b0);
    let c11 = a1.dot(b1);
    let c12 = a1.dot(b2);

    let a1d = a1.dot(rel);
    if intersect_projection!(T, a1d, sb0 * |c10|, sb1 * |c11|, sb2 * |c12| + sa1) {
        return false;
    }

    // -- axis A2
    let c20 = a2.dot(b0);
    let c21 = a2.dot(b1);
    let c22 = a2.dot(b2);

    let a2d = a2.dot(rel);
    if intersect_projection!(T, a2d, sb0 * |c20|, sb1 * |c21|, sb2 * |c22| + sa2) {
        return false;
    }


    // -- axis B0, B1, B2
           !intersect_projection!(T, b0.dot(rel), sa0 * |c00|, sa1 * |c10|, sa2 * |c20| + sb0)
        && !intersect_projection!(T, b1.dot(rel), sa0 * |c01|, sa1 * |c11|, sa2 * |c21| + sb1)
        && !intersect_projection!(T, b2.dot(rel), sa0 * |c02|, sa1 * |c12|, sa2 * |c22| + sb2)
    // -- axis A0 x B0, A0 x B1, A0 x B2
        && !intersect_projection!(T, c10 * a2d - c20 * a1d, sa1 * |c20|, sa2 * |c10|, sb1 * |c02|, sb2 * |c01|)
        && !intersect_projection!(T, c11 * a2d - c21 * a1d, sa1 * |c21|, sa2 * |c11|, sb0 * |c02|, sb2 * |c00|)
        && !intersect_projection!(T, c12 * a2d - c22 * a1d, sa1 * |c22|, sa2 * |c12|, sb0 * |c01|, sb1 * |c00|)
    // -- axis A1 x B0, A1 x B1, A1 x B2
        && !intersect_projection!(T, c20 * a0d - c00 * a2d, sa0 * |c20|, sa2 * |c00|, sb1 * |c12|, sb2 * |c11|)
        && !intersect_projection!(T, c21 * a0d - c01 * a2d, sa0 * |c21|, sa2 * |c01|, sb0 * |c12|, sb2 * |c10|)
        && !intersect_projection!(T, c22 * a0d - c02 * a2d, sa0 * |c22|, sb2 * |c02|, sb0 * |c11|, sb1 * |c10|)
    // -- axis A2 x B0, A2 x B1, A2 x B2
        && !intersect_projection!(T, c00 * a1d - c10 * a0d, sa0 * |c10|, sa1 * |c00|, sb1 * |c22|, sb2 * |c21|)
        && !intersect_projection!(T, c01 * a1d - c11 * a0d, sa0 * |c11|, sa1 * |c01|, sb0 * |c22|, sb2 * |c20|)
        && !intersect_projection!(T, c02 * a1d - c12 * a0d, sa0 * |c12|, sa1 * |c02|, sb0 * |c21|, sb1 * |c20|)
}

/// OBB-OBB non-intersection test using the separation axis theorem in two spacial dimensions.
pub fn intersects_obb_obb_2d<T: BaseFloat>(
    a0: &Vector2<T>, a1: &Vector2<T>,
    b0: &Vector2<T>, b1: &Vector2<T>,
    rel: &Vector2<T>,
    sa0: T, sa1: T,
    sb0: T, sb1: T,
) -> bool {

    // -- axis A0
    let c00 = a0.dot(b0);
    let c01 = a0.dot(b1);

    let a0d = a0.dot(rel);
    if intersect_projection!(T, a0d, sb0 * |c00|, sb1 * |c01| + sa0) {
        return false;
    }

    // -- axis A1
    let c10 = a1.dot(b0);
    let c11 = a1.dot(b1);

    let a1d = a1.dot(rel);
    if intersect_projection!(T, a1d, sb0 * |c10|, sb1 * |c11| + sa1) {
        return false;
    }

    // -- axis B0, B1, B2
           !intersect_projection!(T, b0.dot(rel), sa0 * |c00|, sa1 * |c10| + sb0)
        && !intersect_projection!(T, b1.dot(rel), sa0 * |c01|, sa1 * |c11| + sb1)
    // -- axis A2 x B0, A2 x B1
        && !intersect_projection!(T, c00 * a1d - c10 * a0d, sa0 * |c10|, sa1 * |c00|)
        && !intersect_projection!(T, c01 * a1d - c11 * a0d, sa0 * |c11|, sa1 * |c10|)
}


/// OBB-AABB non-intersection test using the separation axis theorem in three spacial dimensions.
pub fn intersects_obb_aabb<T: BaseFloat>(
    a0: &Vector3<T>, a1: &Vector3<T>, a2: &Vector3<T>,
    rel: &Vector3<T>,
    sa0: T, sa1: T, sa2: T,
    sb0: T, sb1: T, sb2: T,
) -> bool {

    // -- axis A0
    let c00 = a0.x;
    let c01 = a0.y;
    let c02 = a0.z;

    let a0d = a0.dot(rel);
    if intersect_projection!(T, a0d, sb0 * |c00|, sb1 * |c01|, sb2 * |c02| + sa0) {
        return false;
    }

    // -- axis A1
    let c10 = a1.x;
    let c11 = a1.y;
    let c12 = a1.z;

    let a1d = a1.dot(rel);
    if intersect_projection!(T, a1d, sb0 * |c10|, sb1 * |c11|, sb2 * |c12| + sa1) {
        return false;
    }

    // -- axis A2
    let c20 = a2.x;
    let c21 = a2.y;
    let c22 = a2.z;

    let a2d = a2.dot(rel);
    if intersect_projection!(T, a2d, sb0 * |c20|, sb1 * |c21|, sb2 * |c22| + sa2) {
        return false;
    }


    // -- axis B0, B1, B2
           !intersect_projection!(T, rel.x, sa0 * |c00|, sa1 * |c10|, sa2 * |c20| + sb0)
        && !intersect_projection!(T, rel.y, sa0 * |c01|, sa1 * |c11|, sa2 * |c21| + sb1)
        && !intersect_projection!(T, rel.z, sa0 * |c02|, sa1 * |c12|, sa2 * |c22| + sb2)
    // -- axis A0 x B0, A0 x B1, A0 x B2
        && !intersect_projection!(T, c10 * a2d - c20 * a1d, sa1 * |c20|, sa2 * |c10|, sb1 * |c02|, sb2 * |c01|)
        && !intersect_projection!(T, c11 * a2d - c21 * a1d, sa1 * |c21|, sa2 * |c11|, sb0 * |c02|, sb2 * |c00|)
        && !intersect_projection!(T, c12 * a2d - c22 * a1d, sa1 * |c22|, sa2 * |c12|, sb0 * |c01|, sb1 * |c00|)
    // -- axis A1 x B0, A1 x B1, A1 x B2
        && !intersect_projection!(T, c20 * a0d - c00 * a2d, sa0 * |c20|, sa2 * |c00|, sb1 * |c12|, sb2 * |c11|)
        && !intersect_projection!(T, c21 * a0d - c01 * a2d, sa0 * |c21|, sa2 * |c01|, sb0 * |c12|, sb2 * |c10|)
        && !intersect_projection!(T, c22 * a0d - c02 * a2d, sa0 * |c22|, sb2 * |c02|, sb0 * |c11|, sb1 * |c10|)
    // -- axis A2 x B0, A2 x B1, A2 x B2
        && !intersect_projection!(T, c00 * a1d - c10 * a0d, sa0 * |c10|, sa1 * |c00|, sb1 * |c22|, sb2 * |c21|)
        && !intersect_projection!(T, c01 * a1d - c11 * a0d, sa0 * |c11|, sa1 * |c01|, sb0 * |c22|, sb2 * |c20|)
        && !intersect_projection!(T, c02 * a1d - c12 * a0d, sa0 * |c12|, sa1 * |c02|, sb0 * |c21|, sb1 * |c20|)
}

/// OBB-AABB non-intersection test using the separation axis theorem in two spacial dimensions.
pub fn intersects_obb_aabb_2d<T: BaseFloat>(
    a0: &Vector2<T>, a1: &Vector2<T>,
    rel: &Vector2<T>,
    sa0: T, sa1: T,
    sb0: T, sb1: T,
) -> bool {

    // -- axis A0
    let c00 = a0.x;
    let c01 = a0.y;

    let a0d = a0.dot(rel);
    if intersect_projection!(T, a0d, sb0 * |c00|, sb1 * |c01| + sa0) {
        return false;
    }

    // -- axis A1
    let c10 = a1.x;
    let c11 = a1.y;

    let a1d = a1.dot(rel);
    if intersect_projection!(T, a1d, sb0 * |c10|, sb1 * |c11| + sa1) {
        return false;
    }

    // -- axis B0, B1, B2
    !intersect_projection!(T, rel.x, sa0 * |c00|, sa1 * |c10| + sb0)
        && !intersect_projection!(T, rel.y, sa0 * |c01|, sa1 * |c11| + sb1)

        && !intersect_projection!(T, c00 * a1d - c10 * a0d, sa0 * |c10|, sa1 * |c00|)
        && !intersect_projection!(T, c01 * a1d - c11 * a0d, sa0 * |c11|, sa1 * |c10|)
}

/// AABB-AABB non-intersection test using the separation axis theorem in arbitrary spacial
/// dimensions.
pub fn intersects_aabb_aabb<T: BaseFloat, const DIM: usize>(
    min0: &SVector<T, DIM>, max0: &SVector<T, DIM>,
    min1: &SVector<T, DIM>, max1: &SVector<T, DIM>,
) -> bool {
    for i in 0..DIM {
        if max1[i] < min0[i] || min1[i] > max0[i] {
            return false;
        }
    }
    true

    // max1.x >= min0.x && min1.x <= max0.x
    //     && max1.y >= min0.y && min1.y <= max0.y
    //     && max1.z >= min0.z && min1.z <= max0.z
}
