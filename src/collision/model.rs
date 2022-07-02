use std::ops::{Index, IndexMut};
use nalgebra::SVector;
use crate::collision::collision_primitive::CollisionPrimitive;
use crate::helper::BaseFloat;
use crate::system::inertia::Transformer;

pub struct VertexBuffer<T, const DIM: usize> {
    vertices: Vec<SVector<T, DIM>>
}

impl<T> VertexBuffer<T, 3>
where T: BaseFloat {
    pub fn transformed(&self, transform: &Transformer<T>) -> Self {
        VertexBuffer {
            vertices: self.vertices.iter()
                .map(|d| transform.trafo_point(d)).collect()
        }
    }

    pub fn transform_mut(&mut self, transform: &Transformer<T>) {
        for v in self.vertices.iter_mut() {
            *v = transform.trafo_point(v);
        }
    }
}

impl<T, const DIM: usize> Index<usize> for VertexBuffer<T, DIM> {
    type Output = SVector<T, DIM>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.vertices[index]
    }
}

impl<T, const DIM: usize> IndexMut<usize> for VertexBuffer<T, DIM> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.vertices[index]
    }
}

pub struct IndexBuffer {
    indices: Vec<usize>
}

impl Index<usize> for IndexBuffer {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.indices[index]
    }
}

pub struct PhysicsMesh<T, Primitive: CollisionPrimitive<T, DIM>, const DIM: usize> {
    vbo: VertexBuffer<T, DIM>,
    ibo: IndexBuffer,
    prim: Primitive
}

impl<T, Primitive: CollisionPrimitive<T, DIM>, const DIM: usize> PhysicsMesh<T, Primitive, DIM> {
    /// Returns the vertex corresponding to the specified index id. The corresponding inner call
    /// structure is
    /// ``
    /// vbo[ibo[idx]]
    /// ``
    pub fn vertex(&self, idx: usize) -> &SVector<T, DIM> {
        &self.vbo[self.ibo[idx]]
    }
}
