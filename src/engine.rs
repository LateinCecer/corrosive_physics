use std::collections::HashMap;
use std::ops::{Index, IndexMut};
use std::sync::Arc;
use crate::collision::Collider;
use crate::helper::BaseFloat;
use crate::system::object::{PhyEntity, PhyEntityID};
use crate::volume::bvh::VecPool;
use crate::volume::tlas::{TLAS, TLASElement, TLASNode};
use parking_lot::{RawRwLock, RwLock};
use parking_lot::lock_api::{RwLockReadGuard, RwLockWriteGuard};


pub struct PERef<T: BaseFloat> {
    arc: Option<Arc<RwLock<PhysicsEngine<T>>>>
}

impl<T: BaseFloat> PERef<T> {
    pub fn new(engine: PhysicsEngine<T>) -> Self {
        PERef {
            arc: Some(Arc::new(RwLock::new(engine)))
        }
    }

    pub fn lock(&self) -> RwLockReadGuard<'_, RawRwLock, PhysicsEngine<T>> {
        match self.arc.as_ref() {
            Some(a) => a.read(),
            None => panic!("Physics Engine is not initiated")
        }
    }

    pub fn lock_mut(&self) -> RwLockWriteGuard<'_, RawRwLock, PhysicsEngine<T>> {
        match self.arc.as_ref() {
            Some(a) => a.write(),
            None => panic!("Physics Engine is not initialed")
        }
    }
}

impl<T: BaseFloat> Default for PERef<T> {
    fn default() -> Self {
        PERef {
            arc: None
        }
    }
}




pub static mut PHYSICS_ENGINE : PERef<f64> = PERef { arc: None };


pub struct PhysicsEngine<T: BaseFloat> {
    collider: HashMap<usize, Box<dyn Collider<T, 3>>>,
    pub world: TLAS<T, PhyEntity<T>, VecPool<TLASNode<T, 3>>, VecPool<PhyEntity<T>>, 3>
}

impl<T: BaseFloat> PhysicsEngine<T> {
    pub fn new() -> Self {
        PhysicsEngine {
            collider: HashMap::new(),
            world: TLAS::new(64),
        }
    }

    pub fn query_colliders(&self, id: PhyEntityID) -> Vec<&PhyEntity<T>> {
        let header = &self.world.blas()[id.entity_id];
        let colliders = self.world
            .intersect(header.bounding_volume(), 0);
        colliders
    }
}

impl<T: BaseFloat> Index<PhyEntityID> for PhysicsEngine<T> {
    type Output = PhyEntity<T>;

    fn index(&self, index: PhyEntityID) -> &Self::Output {
        &self.world.blas()[index.entity_id]
    }
}

impl<T: BaseFloat> IndexMut<PhyEntityID> for PhysicsEngine<T> {
    fn index_mut(&mut self, index: PhyEntityID) -> &mut Self::Output {
        &mut self.world.blas_mut()[index.entity_id]
    }
}

impl PhysicsEngine<f64> {
    pub unsafe fn init_global(engine : PhysicsEngine<f64>) {
        PHYSICS_ENGINE = PERef::new(engine);
    }

    pub fn global() -> RwLockReadGuard<'static, RawRwLock, PhysicsEngine<f64>> {
        unsafe {
            PHYSICS_ENGINE.lock()
        }
    }

    pub fn global_mut() -> RwLockWriteGuard<'static, RawRwLock, PhysicsEngine<f64>> {
        unsafe {
            PHYSICS_ENGINE.lock_mut()
        }
    }
}
