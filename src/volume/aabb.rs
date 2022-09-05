use nalgebra::{DimMin, SVector, Vector3};
use crate::helper::{BaseFloat, separated_axis};
use crate::volume::{BoundingVolume, BVIntersector};
use crate::volume::oriented::OBB;

/// Axis aligned bounding box.
#[derive(Clone, Copy, Debug)]
pub struct AABB<T, const DIM: usize> {
    pub min: SVector<T, DIM>,
    pub max: SVector<T, DIM>
}


impl<T, const DIM: usize> AABB<T, DIM>
where T: BaseFloat {
    /// Creates a new AABB instance, where the min point is set the the base float's max and the
    /// max point is set to the base float's min value. Initializing the boundaries this way makes
    /// it easier to fit the AABB to some boundary volume.
    pub fn new() -> Self {
        AABB {
            min: SVector::repeat(T::MAX),
            max: SVector::repeat(T::MIN),
        }
    }

    /// Resets the min and max values of the AABB to the base float's max-, min values (see
    /// `new()` for more information).
    pub fn reset(&mut self) {
        self.min = SVector::repeat(T::MAX);
        self.max = SVector::repeat(T::MIN);
    }

    /// Adjusts the boundaries of the AABB to wrap the two specified AABBs.
    pub fn adjust(&mut self, left: &AABB<T, DIM>, right: &AABB<T, DIM>) {
        for i in 0..DIM {
            self.min[i] = T::min(left.min[i], right.min[i]);
            self.max[i] = T::max(left.max[i], right.max[i]);
        }
    }

    /// Grows the size of this AABB to wrap the specified `other` AABB. As the name of this method
    /// implies, this process can only grow the AABB, not shrink it to any extend.
    pub fn grow_other(&mut self, other: &AABB<T, DIM>) {
        if T::is_finite(&other.min[0]) {
            for i in 0..DIM {
                self.min[i] = T::min(self.min[i], other.min[i]);
                self.max[i] = T::max(self.max[i], other.max[i]);
            }
        }
    }

    /// Grows the size of the AABB to wrap the specified point `p`. As the name of this method
    /// implies, this process can only grow the AABB, not shrink it to any extend.
    pub fn grow(&mut self, p: &SVector<T, DIM>) {
        for i in 0..DIM {
            self.min[i] = T::min(self.min[i], p[i]);
            self.max[i] = T::max(self.max[i], p[i]);
        }
    }

    /// Grows the `min` bounds of this AABB to fit the specified point. If the point lies to the
    /// positive side of the center of the AABB, this method will not change the AABB and the point
    /// will not be included.
    pub fn grow_min(&mut self, p: &SVector<T, DIM>) {
        for i in 0..DIM {
            self.min[i] = T::min(self.min[i], p[i]);
        }
    }

    /// Grows the `max` bounds of this AABB to fit the specified point. If the point lies to the
    /// negative side of the center of the AABB, this method will not change the AABB and the point
    /// will not be included.
    pub fn grow_max(&mut self, p: &SVector<T, DIM>) {
        for i in 0..DIM {
            self.max[i] = T::max(self.max[i], p[i]);
        }
    }
}

impl<T: BaseFloat, const DIM: usize> BoundingVolume<T, DIM> for AABB<T, DIM> {
    fn center(&self) -> SVector<T, DIM> {
        (self.min + self.max) * T::half()
    }

    fn area(&self) -> T {
        let size = self.max - self.min;
        let mut sum = T::zero();
        for i in 0..DIM {
            sum += size[i] * size[(i + 1) % DIM];
        }
        sum
    }

    fn min(&self) -> SVector<T, DIM> {
        self.min.clone()
    }

    fn max(&self) -> SVector<T, DIM> {
        self.max.clone()
    }

    fn size(&self) -> SVector<T, DIM> {
        self.max - self.min
    }

    fn half_size(&self) -> SVector<T, DIM> {
        (self.max - self.min) * T::half()
    }
}

impl<T: BaseFloat, const DIM: usize> BVIntersector<T, AABB<T, DIM>, DIM> for AABB<T, DIM> {
    fn intersects(&self, other: &AABB<T, DIM>) -> bool {
        separated_axis::intersects_aabb_aabb(
            &self.min, &self.max,
            &other.min, &other.max
        )
    }
}

impl<T: BaseFloat> BVIntersector<T, OBB<T>, 3> for AABB<T, 3> {
    fn intersects(&self, other: &OBB<T>) -> bool {
        // AABB-OBB intersections are already implemented for the OBB struct. Use that
        // implementation here to avoid duplications.
        other.intersects(self)
    }
}

impl<T: BaseFloat, const DIM: usize> BVIntersector<T, SVector<T, DIM>, DIM> for AABB<T, DIM> {
    fn intersects(&self, other: &SVector<T, DIM>) -> bool {
        // AABB-point intersections are already implemented for the SVector struct. Use that
        // implementation here to avoid duplications
        other.intersects(self)
    }
}
