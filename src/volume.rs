use nalgebra::SVector;

pub mod aabb;
pub mod tlas;
pub mod bvh;
pub mod bvh_splitting;
pub mod oriented;
pub mod point;


pub trait BoundingVolume<T, const DIM: usize> {
    /// Returns the center point of the bounding volume.
    ///
    /// This center point is to be seen as the position reference point for the volume. It's
    /// definition may be arbitrary for different bounding volumes.
    fn center(&self) -> SVector<T, DIM>;

    /// Returns a value representative for the surface area of the bounding volume.
    ///
    /// This method does not have to return the actual surface area, but only a value that is
    /// representative for it. This way, computational cost may be saved by returning an
    /// approximation for the total area. The return value of this method may be used, for example,
    /// for a traversal-cost analysis using BVH construction.
    fn area(&self) -> T;

    /// Returns the minimal euclidean x, y & z coordinate of this bounding volume.
    ///
    /// This may be used to find a AABB fit around the volume.
    fn min(&self) -> SVector<T, DIM>;

    /// Returns the maximal euclidean x, y & z coordinate of this bounding volume.
    ///
    /// This may be used to find a AABB fit around the volume.
    fn max(&self) -> SVector<T, DIM>;

    /// Returns the total size of the bounding volume. For more complex geometric shapes, this
    /// method may return an approximate wrapping box size around the center of the BV.
    fn size(&self) -> SVector<T, DIM>;

    /// Returns half of the total size of the bounding volume. A call to this method is equivalent
    /// to `size() * 0.5` but that may be more inefficient in certain contexts, like OBBs, which
    /// store the half size of the box directly.
    fn half_size(&self) -> SVector<T, DIM>;
}

pub trait BVIntersector<T, O: BoundingVolume<T, DIM>, const DIM: usize> {
    /// Returns true, if there is an overlap between the implementation of this trait and the
    /// specified bounding volume.
    fn intersects(&self, other: &O) -> bool;
}
