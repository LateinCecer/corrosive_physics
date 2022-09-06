use std::marker::PhantomData;
use std::mem;
use std::ops::{Index, IndexMut};
use nalgebra::{SVector};
use crate::helper::BaseFloat;
use crate::volume::aabb::AABB;
use crate::volume::bvh::VecPool;
use crate::volume::{BoundingVolume, BVIntersector};

#[derive(Clone, Debug)]
pub struct TLASNode<T: BaseFloat, const DIM: usize> {
    aabb: AABB<T, DIM>,
    left_right: u32,
    blas: u32,
}

impl<T: BaseFloat, const DIM: usize> TLASNode<T, DIM> {

    /// Creates a new TLAS node instance, where every value is initiated with zero.
    pub fn new() -> Self {
        TLASNode {
            aabb: AABB::new(),
            left_right: 0,
            blas: 0,
        }
    }

    /// Copies all values from the specified `other` TLAS node.
    pub fn copy_from(&mut self, other: &Self) {
        self.aabb = other.aabb.clone();
        self.left_right = other.left_right.clone();
        self.blas = other.blas.clone();
    }

    pub fn aabb(&self) -> &AABB<T, DIM> {
        &self.aabb
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
    type BV: BoundingVolume<T, DIM>;

    /// Returns an AABB that fully wraps around the TLAS element.
    fn wrap(&self) -> AABB<T, DIM>;
    fn bounding_volume(&self) -> &Self::BV;
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
        assert_eq!(target_len, 1);

        let first = self.vec.remove(0);
        self.vec.clear();
        self.vec.push(first);
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

impl<T, B, const DIM: usize> TLAS<T, B, VecPool<TLASNode<T, DIM>>, VecPool<B>, DIM>
where T: BaseFloat,
      B: TLASElement<T, DIM> + Sized {

    pub fn new(cap: usize) -> Self {
        let mut tlas = TLAS {
            nodes: VecPool::with_capacity(cap * 2),
            blas: VecPool::with_capacity(cap),
            _t: PhantomData::default(),
            _b: PhantomData::default(),
        };
        tlas.nodes.push(TLASNode {
            aabb: AABB::new(),
            blas: 0,
            left_right: 0
        });

        tlas
    }
}


#[derive(Clone, Copy, Debug)]
enum CondIndex {
    S(usize),
    N
}

impl PartialEq for CondIndex {
    fn eq(&self, other: &Self) -> bool {
        match self {
            CondIndex::S(i) => {
                match other {
                    CondIndex::S(j) => *i == *j,
                    CondIndex::N => false
                }
            },
            CondIndex::N => false
        }
    }
}

impl CondIndex {
    fn matches_index(&self, i: usize) -> bool {
        match self {
            CondIndex::S(j) => *j == i,
            CondIndex::N => false
        }
    }

    fn unwrap(&self) -> &usize {
        match self {
            CondIndex::S(i) => i,
            CondIndex::N => panic!("Invalid index")
        }
    }
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

    pub fn blas_mut(&mut self) -> &mut BlasPool {
        &mut self.blas
    }

    pub fn refit(&mut self) {
        // since a parent node is always further to the back of the tree, we can loop through here
        // front-to-back
        for i in 1..self.nodes.size() {
            let node = &self.nodes[i];
            if node.is_leaf() {
                self.nodes[i].aabb = self.blas[node.blas as usize].wrap();
            } else {
                let left_child = &self.nodes[node.get_left_child() as usize].aabb;
                let right_child = &self.nodes[node.get_right_child() as usize].aabb;

                let mut aabb = node.aabb.clone();
                aabb.adjust(left_child, right_child);
                self.nodes[i].aabb = aabb;
            }
        }
    }

    /// Rebuilds the TLAS bottom up.
    pub fn build(&mut self) {
        let mut node_idx = Vec::<usize>::with_capacity(self.blas.size());
        let mut node_indices = self.blas.size();

        // set leaf nodes
        self.nodes.trim(1);
        for i in 0..self.blas.size() {
            node_idx.push(self.nodes.size());
            self.nodes.push(TLASNode {
                aabb: self.blas[i].wrap(),
                blas: i as u32,
                left_right: 0,
            });
        }

        // eprintln!("init node len: {}", self.nodes.size());

        // use agglomerative clustering to build the TLAS (bottom-to-top)
        let mut a = 0_i32;
        let mut b = self.find_best_match(&node_idx, node_indices, a);
        while node_indices > 1 {
            let c = self.find_best_match(&node_idx, node_indices, b);
            if a == c {
                let node_idx_a = node_idx[a as usize];
                let node_idx_b = node_idx[b as usize];

                let node_a = &self.nodes[node_idx_a];
                let node_b = &self.nodes[node_idx_b];
                node_idx[a as usize] = self.nodes.size();
                node_idx[b as usize] = node_idx[node_indices - 1];


                let mut aabb = AABB::new();
                aabb.adjust(&node_a.aabb, &node_b.aabb);
                self.nodes.push(TLASNode {
                    left_right: node_idx_a as u32 + ((node_idx_b as u32) << 16),
                    aabb,
                    blas: 0
                });

                node_indices -= 1;
                b = self.find_best_match(&node_idx, node_indices, a);
            } else {
                a = b;
                b = c;
            }
        }
        // eprintln!("nodes.len() = {}", self.nodes.size());

        // set root node
        self.nodes[0] = self.nodes[node_idx[a as usize]].clone();
        // eprintln!("nodes:");
        // for i in 0..self.nodes.size() {
        //     eprintln!("  [{}]: {:?}     >>   left={},    >>   right={}",
        //               i, self.nodes[i],
        //               self.nodes[i].get_left_child(),
        //               self.nodes[i].get_right_child());
        // }
    }

    /// Finds the most cost-effective clustering partner for the node with id `list[a]`. For this,
    /// the `n` first entries in `list` are considered.
    fn find_best_match(&self, list: &Vec<usize>, n: usize, a: i32) -> i32 {
        let mut smallest = T::MAX;
        let mut best_b = -1_i32;

        for b in 0..n {
            if b as i32 == a {
                continue;
            }

            let a_node = &self.nodes[list[a as usize]];
            let b_node = &self.nodes[list[b]];

            // calc wrapping node sizes
            let mut size = SVector::<T, DIM>::zeros();
            for i in 0..DIM {
                size[i] = T::max(a_node.aabb.max[i], b_node.aabb.max[i])
                    - T::min(a_node.aabb.min[i], b_node.aabb.min[i]);
            }

            // calc surface area estimate for cost analysis
            let mut surface_area = T::zero();
            for i in 0..DIM {
                surface_area += size[i] * size[(i + 1) % DIM];
            }


            if surface_area < smallest {
                smallest = surface_area;
                best_b = b as i32;
            }
        }
        return best_b;
    }

    pub fn intersect<I: BVIntersector<T, B::BV, DIM> + BVIntersector<T, AABB<T, DIM>, DIM>>(
        &self, intersector: &I, node_idx: usize
    ) -> Vec<&B> {

        let mut v = Vec::<&B>::with_capacity(64);

        let mut node = &self.nodes[node_idx];
        let mut stack = [node; 64];
        let mut stack_ptr = 0usize;

        loop {
            if node.is_leaf() {
                if intersector.intersects(self.blas[node.blas as usize].bounding_volume()) {
                    v.push(&self.blas[node.blas as usize]);
                }

                if stack_ptr == 0 {
                    break;
                } else {
                    stack_ptr -= 1;
                    node = stack[stack_ptr];
                }
            } else {
                let mut child1 = &self.nodes[node.get_left_child() as usize];
                let mut child2 = &self.nodes[node.get_right_child() as usize];

                let mut inter1 = intersector.intersects(&child1.aabb);
                let mut inter2 = intersector.intersects(&child2.aabb);
                // -(for ray intersections, do ray sorting here)
                if !inter1 {
                    // if child 1 does not intersect the intersector, swap with child 2
                    mem::swap(&mut child1, &mut child2);
                    mem::swap(&mut inter1, &mut inter2);
                }

                if !inter1 {
                    // both children do not intersect. Checkout stack
                    if stack_ptr == 0 {
                        break;
                    } else {
                        stack_ptr -= 1;
                        node = stack[stack_ptr];
                    }
                } else {
                    node = child1;
                    // checkout child 1 first and save child 2 for later
                    if inter2 {
                        stack[stack_ptr] = child2;
                        stack_ptr += 1;
                    }
                }
            }
        }
        v
    }
}
