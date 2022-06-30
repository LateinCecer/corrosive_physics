use std::marker::PhantomData;
use std::ops::{Index, IndexMut};
use nalgebra::{SVector, Vector3};
use crate::helper::BaseFloat;
use crate::volume::aabb::AABB;
use crate::volume::bvh::VecPool;

#[derive(Clone, Debug)]
pub struct TLASNode<T: BaseFloat, const DIM: usize> {
    min: SVector<T, DIM>,
    max: SVector<T, DIM>,
    left_right: u32,
    blas: u32,
}

impl<T: BaseFloat, const DIM: usize> TLASNode<T, DIM> {

    /// Creates a new TLAS node instance, where every value is initiated with zero.
    pub fn new() -> Self {
        TLASNode {
            min: SVector::zeros(),
            max: SVector::zeros(),
            left_right: 0,
            blas: 0,
        }
    }

    /// Copies all values from the specified `other` TLAS node.
    pub fn copy_from(&mut self, other: &Self) {
        self.min = other.min.clone();
        self.max = other.max.clone();
        self.left_right = other.left_right.clone();
        self.blas = other.blas.clone();
    }

    /// Returns true, only if the node is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left_right == 0u32
    }

    /// Returns the pool index of the left child of this node
    pub fn get_left_child(&self) -> u16 {
        (&self.left_right >> 16) as u16
    }

    /// Returns the pool index of the right child of this node.
    pub fn get_right_child(&self) -> u16 {
        (&self.left_right & 0xFFFF) as u16
    }
}



pub trait TLASPool<T: Sized> : Index<usize, Output=T> + IndexMut<usize, Output=T> {
    /// Pushes an element `el` to the TLASPool. The pushed element is appended to the back of the
    /// pool vector.
    fn push(&mut self, el: T);

    /// Pops and returns the last element of the pool. If the pool is empty, `None` is returned.
    fn pop(&mut self) -> Option<T>;

    /// Returns the amount of elements that is currently storged in the pool.
    fn size(&self) -> usize;

    /// Returns the total capacity of the pool.
    fn capacity(&self) -> usize;

    /// Trims the pool to the specified target length.
    fn trim(&mut self, target_len: usize);

    /// Returns a shared reference to the first element in the pool. If the pool is emtpy, `None` is
    /// returned.
    fn front(&self) -> Option<&T>;

    /// Returns a mutable reference to the first element in the pool. If the pool is empty, `None`
    /// is returned.
    fn front_mut(&mut self) -> Option<&mut T>;

    /// Returns a shared reference to the last element of the pool. If the pool is empty, `None`
    /// is returned.
    fn back(&self) -> Option<&T>;

    /// Returns a mutable reference to the last element of the pool. If the pool is empty, `None`
    /// is returned.
    fn back_mut(&mut self) -> Option<&mut T>;
}

pub trait TLASElement<T: BaseFloat, const DIM: usize> {
    /// Returns an AABB that fully wraps around the TLAS element.
    fn wrap(&self) -> AABB<T, DIM>;
}


impl<T: Sized> TLASPool<T> for VecPool<T> {
    fn push(&mut self, el: T) {
        self.vec.push(el);
    }

    fn pop(&mut self) -> Option<T> {
        self.vec.pop()
    }

    fn size(&self) -> usize {
        self.vec.len()
    }

    fn capacity(&self) -> usize {
        self.vec.capacity()
    }

    fn trim(&mut self, target_len: usize) {
        while self.vec.len() > target_len {
            self.vec.pop();
        }
    }

    fn front(&self) -> Option<&T> {
        self.vec.first()
    }

    fn front_mut(&mut self) -> Option<&mut T> {
        self.vec.first_mut()
    }

    fn back(&self) -> Option<&T> {
        self.vec.last()
    }

    fn back_mut(&mut self) -> Option<&mut T> {
        self.vec.last_mut()
    }
}



/// A TLAS (short for Top Level Acceleration Structure) is a BVH-like hierarchical tree structure
/// that contains instances of `BoundingVolume`. The tree structure can be traversed in O(log n) to
/// query for elements that intersect an appropriate intersector instance. This is useful for
/// rending (specially but not limited to ray tracing), culling and collision detection.
pub struct TLAS<T: BaseFloat, B: Sized, NodePool: TLASPool<TLASNode<T, DIM>>, BlasPool: TLASPool<B>, const DIM: usize> {
    nodes: NodePool,
    blas: BlasPool,

    _t: PhantomData<T>,
    _b: PhantomData<B>,
}

impl<T, B, NodePool, BlasPool, const DIM: usize> TLAS<T, B, NodePool, BlasPool, DIM>
where T: BaseFloat,
      B: TLASElement<T, DIM> + Sized,
      NodePool: TLASPool<TLASNode<T, DIM>>,
      BlasPool: TLASPool<B> {

    /// Returns a shared reference to the `TLASPool` instance that contains the TLAS nodes.
    pub fn nodes(&self) -> &NodePool {
        &self.nodes
    }

    /// Returns a shared reference to the 'TLASPool' instance that contains the BLAS elements for
    /// the TLAS.
    pub fn blas(&self) -> &BlasPool {
        &self.blas
    }

    /// Rebuilds the TLAS bottom up.
    pub fn build(&mut self) {
        let mut node_idx = Vec::<usize>::with_capacity(self.blas.size());
        let node_indices = self.blas.size();

        // set leaf nodes
        self.nodes.trim(1);
        for i in 0..self.blas.size() {
            node_idx.push(self.nodes.size());
            let bounds = self.blas[i].wrap();
            self.nodes.push(TLASNode {
                min: bounds.min,
                max: bounds.max,
                blas: i as u32,
                left_right: 0,
            });
        }

        // use agglomerative clustering to build the TLAS (bottom-to-top)
        let mut a = 0;
        let mut b = match self.find_best_match(&node_idx, node_indices, a) {
            Some(b) => b,
            None => return
        };
        while node_indices > 1 {
            let c = match self.find_best_match(&node_idx, node_indices, b) {
                Some(c) => c,
                None => break
            };
            if a == c {
                let node_idx_a = node_idx[a];
                let node_idx_b = node_idx[b];

                let node_a = &self.nodes[node_idx_a];
                let node_b = &self.nodes[node_idx_b];
                node_idx[a] = self.nodes.size();
                node_idx[b] = node_idx[node_indices - 1];


                let mut min = SVector::<T, DIM>::zeros();
                let mut max = SVector::<T, DIM>::zeros();
                for i in 0..DIM {
                    min[i] = T::min(node_a.min[i], node_b.min[i]);
                    max[i] = T::max(node_a.max[i], node_b.max[i]);
                }

                self.nodes.push(TLASNode {
                    left_right: node_idx_a as u32 | ((node_idx_b as u32) << 16),
                    min,
                    max,
                    blas: 0
                });

                b = match self.find_best_match(&node_idx, node_indices, a) {
                    Some(b) => b,
                    None => break
                };
            } else {
                a = b;
                b = c;
            }
        }
        // set root node
        self.nodes[0] = self.nodes[node_idx[a]].clone();
    }

    /// Finds the most cost-effective clustering partner for the node with id `list[a]`. For this,
    /// the `n` first entries in `list` are considered.
    fn find_best_match(&self, list: &Vec<usize>, n: usize, a: usize) -> Option<usize> {
        let mut smallest = T::MAX;
        let mut best_b = None;

        for b in 0..n {
            if b == a {
                break;
            }

            let a_node = &self.nodes[list[a]];
            let b_node = &self.nodes[list[b]];

            // calc wrapping node sizes
            let mut size = SVector::<T, DIM>::zeros();
            for i in 0..DIM {
                size[i] = T::max(a_node.max[i], b_node.max[i])
                    - T::min(a_node.min[i], b_node.min[i]);
            }

            // calc surface area estimate for cost analysis
            let mut surface_area = T::zero();
            for i in 0..DIM {
                surface_area += size[i] * size[(i + 1) % DIM];
            }
            if surface_area < smallest {
                smallest = surface_area;
                best_b = Some(b);
            }
        }
        return best_b;
    }
}
