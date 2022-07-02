use nalgebra::SVector;

pub struct RayIntersection<T, const DIM: usize> {
    pub pos: SVector<T, DIM>,
    pub normal: SVector<T, DIM>,
    pub prim_id: usize
}

pub struct Ray<T, const DIM: usize> {
    pub d: T,
    pub origin: SVector<T, DIM>,
    pub dir: SVector<T, DIM>,
    pub intersection: Option<RayIntersection<T, DIM>>,
}
