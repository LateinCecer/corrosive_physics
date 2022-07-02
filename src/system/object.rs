use crate::helper::BaseFloat;
use crate::system::inertia::{IS, MassDistribution, Transformer};
use crate::volume::aabb::AABB;
use crate::volume::oriented::OBB;
use crate::volume::tlas::TLASElement;
use bevy::prelude::{Component, Res, Time};
use nalgebra::Vector3;
use crate::volume::BoundingVolume;



#[derive(Clone, PartialEq, Component)]
pub struct PhyEntityID {
    pub world_id: u8,
    pub chunk_id: usize,
    pub entity_id: usize,
}


pub struct PhyEntity<T: BaseFloat> {
    pub id: PhyEntityID,
    pub is: IS<T>,
    collider_id: usize,
    obb: OBB<T>,
}

impl<T: BaseFloat> PhyEntity<T> {
    pub fn cube(id: PhyEntityID, size: Vector3<T>) -> Self {
        PhyEntity {
            id,
            is: IS::new(Vector3::zeros(), Vector3::zeros(), Transformer::default(), MassDistribution::default()),
            collider_id: 0,
            obb: OBB { half_size: size.scale(T::half()), transform: Transformer::default() }
        }
    }

    pub fn sync(&mut self) {
        self.is.sync();
        self.obb.transform = self.is.state.clone();
    }

    pub fn tick(&mut self, time: &Res<Time>) {
        // TODO

    }
}

impl<T: BaseFloat> TLASElement<T, 3> for PhyEntity<T> {
    type BV = OBB<T>;

    fn wrap(&self) -> AABB<T, 3> {
        AABB {
            min: self.obb.min(),
            max: self.obb.max(),
        }
    }

    fn bounding_volume(&self) -> &Self::BV {
        &self.obb
    }
}
