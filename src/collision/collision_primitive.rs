use nalgebra::SVector;
use crate::collision::intersection::Ray;
use crate::collision::model::VertexBuffer;
use crate::volume::aabb::AABB;

pub type Edge = (usize, usize);

pub trait CollisionPrimitive<T, const DIM: usize> {
    fn indices(&self) -> &[usize];
    fn edges(&self) -> &[Edge];

    fn centroid(&self, id: usize, vbo: &VertexBuffer<T, DIM>) -> SVector<T, DIM>;
    fn wrap(&self, id: usize, vbo: &VertexBuffer<T, DIM>) -> AABB<T, DIM>;
    fn intersect_ray(&self, id: usize, vbo: &VertexBuffer<T, DIM>, ray: &mut Ray<T, DIM>);
}


