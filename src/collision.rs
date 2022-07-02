mod collider;
mod collision_primitive;
mod intersection;
mod model;

use nalgebra::{UnitQuaternion, Vector3};
use crate::helper::BaseFloat;
use crate::volume::{BoundingVolume, BVIntersector};
use crate::volume::aabb::AABB;
use crate::volume::oriented::OBB;

pub enum MovementTrigger<T> {
    Translation(Vector3<T>),
    Rotation(UnitQuaternion<T>),
    Scale(Vector3<T>),
}

pub trait Collider<T, const DIM: usize> {
    fn wrap(&self) -> &dyn BoundingVolume<T, DIM>;
}

pub struct ColliderVolume<'a, T> {
    trigger: MovementTrigger<T>,
    collider: &'a dyn Collider<T, 3>,
}


impl<'a, T: BaseFloat> BVIntersector<T, AABB<T, 3>, 3> for ColliderVolume<'a, T> {
    fn intersects(&self, other: &AABB<T, 3>) -> bool {
        todo!()
    }
}

impl<'a, T: BaseFloat> BVIntersector<T, OBB<T>, 3> for ColliderVolume<'a, T> {
    fn intersects(&self, other: &OBB<T>) -> bool {
        todo!()
    }
}

impl<'a, T: BaseFloat> BVIntersector<T, Vector3<T>, 3> for ColliderVolume<'a, T> {
    fn intersects(&self, other: &Vector3<T>) -> bool {
        todo!()
    }
}
