use nalgebra::SVector;
use crate::helper::BaseFloat;
use crate::volume::{BoundingVolume, BVIntersector};
use crate::volume::aabb::AABB;
use crate::volume::oriented::OBB;

impl<T: BaseFloat, const DIM: usize> BoundingVolume<T, DIM> for SVector<T, DIM> {
    fn center(&self) -> SVector<T, DIM> {
        self.clone()
    }

    fn area(&self) -> T {
        T::zero()
    }

    fn min(&self) -> SVector<T, DIM> {
        self.clone()
    }

    fn max(&self) -> SVector<T, DIM> {
        self.clone()
    }

    fn size(&self) -> SVector<T, DIM> {
        SVector::zeros()
    }

    fn half_size(&self) -> SVector<T, DIM> {
        SVector::zeros()
    }
}

impl<T: BaseFloat, const DIM: usize> BVIntersector<T, AABB<T, DIM>, DIM> for SVector<T, DIM> {
    fn intersects(&self, other: &AABB<T, DIM>) -> bool {
        for i in 0..DIM {
            if self[i] < other.min[i] || self[i] > other.max[i] {
                return false;
            }
        }
        return true;
    }
}

impl<T: BaseFloat> BVIntersector<T, OBB<T>, 3> for SVector<T, 3> {
    fn intersects(&self, other: &OBB<T>) -> bool {
        // point-OBB intersections are already implemented for the OBB-struct. Use that
        // implementation here to avoid duplications
        other.intersects(self)
    }
}
