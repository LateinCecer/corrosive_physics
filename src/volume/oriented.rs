use nalgebra::{SVector, Vector3};
use crate::helper::{BaseFloat, separated_axis};
use crate::system::inertia::Transformer;
use crate::volume::aabb::AABB;
use crate::volume::{BoundingVolume, BVIntersector};

/// An implementation for an oriented bounding box
pub struct OBB<T> {
    half_size: Vector3<T>,
    pub transform: Transformer<T>
}

impl<T: BaseFloat> BoundingVolume<T, 3> for OBB<T> {
    fn center(&self) -> Vector3<T> {
        self.transform.pos + self.transform.trafo_vec(&self.transform.offset)
    }

    fn area(&self) -> T {
        self.half_size.x * self.half_size.y
            + self.half_size.y * self.half_size.z
            + self.half_size.z * self.half_size.x
    }

    fn min(&self) -> Vector3<T> {
        let min = self.transform.trafo_point(&self.half_size);
        let max = self.transform.trafo_point(&(-self.half_size));
        Vector3::new(
            T::min(min.x, max.x),
            T::min(min.y, max.y),
            T::min(min.z, max.z),
        )
    }

    fn max(&self) -> Vector3<T> {
        let min = self.transform.trafo_point(&self.half_size);
        let max = self.transform.trafo_point(&(-self.half_size));
        Vector3::new(
            T::max(min.x, max.x),
            T::max(min.y, max.y),
            T::max(min.z, max.z),
        )
    }

    fn size(&self) -> Vector3<T> {
        self.half_size * T::two()
    }

    fn half_size(&self) -> Vector3<T> {
        self.half_size
    }
}

impl<T: BaseFloat> BVIntersector<T, OBB<T>, 3> for OBB<T> {
    fn intersects(&self, other: &OBB<T>) -> bool {
        separated_axis::intersects_obb_obb(
            &self.transform.right(),
            &self.transform.up(),
            &self.transform.forward(),
            &other.transform.right(),
            &other.transform.up(),
            &other.transform.forward(),
            &(other.center() - self.center()),
            self.half_size.x, self.half_size.y, self.half_size.z,
            other.half_size.x, other.half_size.y, self.half_size.z
        )
    }
}

impl<T: BaseFloat> BVIntersector<T, AABB<T, 3>, 3> for OBB<T> {
    fn intersects(&self, other: &AABB<T, 3>) -> bool {
        let other_half_size = other.half_size();
        separated_axis::intersects_obb_aabb(
            &self.transform.right(),
            &self.transform.up(),
            &self.transform.forward(),
            &(other.center() - self.center()),
            self.half_size.x, self.half_size.y, self.half_size.z,
            other_half_size.x, other_half_size.y, other_half_size.z
        )
    }
}

impl<T: BaseFloat> BVIntersector<T, SVector<T, 3>, 3> for OBB<T> {
    fn intersects(&self, other: &SVector<T, 3>) -> bool {
        // transform the point into the reference system of the obb
        let rel = self.transform.inv_trafo_point(other);
        for i in 0..3 {
            if rel[i] < -self.half_size[i] || rel[i] > self.half_size[i] {
                return false;
            }
        }
        true
    }
}
