use std::marker::PhantomData;
use std::mem;
use std::ops::{Index, IndexMut};
use nalgebra::SVector;
use crate::helper::BaseFloat;
use crate::volume::aabb::AABB;
use crate::volume::{BoundingVolume, BVIntersector};
use crate::volume::bvh_splitting::BVHSplitting;


/// Basic data structure for a BVH node.
pub struct BVHNode<T, const DIM: usize> {
    aabb: AABB<T, DIM>,
    left_first: usize,
    num_prims: usize,
}

impl<T, const DIM: usize> BVHNode<T, DIM> {
    pub fn left_child(&self) -> usize {
        self.left_first
    }

    /// Returns the right child of the bvh node.
    pub fn right_child(&self) -> usize {
        self.left_first + 1
    }

    /// Returns true, if the node is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.num_prims > 0
    }

    pub fn num_prims(&self) -> &usize {
        &self.num_prims
    }

    pub fn aabb(&self) -> &AABB<T, DIM> {
        &self.aabb
    }
}

impl<T, const DIM: usize> BVHNode<T, DIM>
where T: BaseFloat {
    pub fn new() -> Self {
        BVHNode {
            aabb: AABB::new(),
            left_first: 0,
            num_prims: 0,
        }
    }
}



pub trait BVHPool<T, const DIM: usize> : Index<usize, Output=BVHNode<T, DIM>> + IndexMut<usize, Output=BVHNode<T, DIM>> {
    /// Returns the capacity of the BVH pool
    fn capacity(&self) -> usize;
}

pub trait BVHElement<T, const DIM: usize> : BoundingVolume<T, DIM> {
    /// Returns the geometric centroid for the element
    fn centroid(&self) -> SVector<T, DIM>;

    /// Wraps an AABB box around the element.
    fn wrap(&self) -> AABB<T, DIM>;
}

pub trait BVHElementPool<T, ElementType: BVHElement<T, DIM>, const DIM: usize> : Index<usize, Output=ElementType>
    + IndexMut<usize, Output=ElementType> {

    /// Returns the capacity of the BVH Element pool
    fn capacity(&self) -> usize;

    /// Returns the amount of elements within the element pool
    fn len(&self) -> usize;

    /// Swaps the element at index `i` with the element at index `j`.
    fn swap(&mut self, i: usize, j: usize);
}



/// A `VecPool` is a memory pool implementation based on an `alloc::vec::Vec`.
pub struct VecPool<T: Sized> {
    pub vec: Vec<T>,
}

impl<T: Sized> Index<usize> for VecPool<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.vec[index]
    }
}

impl<T: Sized> IndexMut<usize> for VecPool<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.vec[index]
    }
}

impl<T: Sized, const DIM: usize> BVHPool<T, DIM> for VecPool<BVHNode<T, DIM>> {
    fn capacity(&self) -> usize {
        self.vec.capacity()
    }
}

impl<T: Sized> VecPool<T> {
    pub fn with_capacity(n: usize) -> Self {
        VecPool {
            vec: Vec::with_capacity(n)
        }
    }

    pub fn new() -> Self {
        VecPool {
            vec: Vec::new()
        }
    }

    pub fn push(&mut self, item: T) {
        self.vec.push(item);
    }

    pub fn remove(&mut self, idx: usize) -> T {
        self.vec.remove(idx)
    }

    pub fn clear(&mut self) {
        self.vec.clear();
    }
}

impl<T: Sized, E: BVHElement<T, DIM>, const DIM: usize> BVHElementPool<T, E, DIM> for VecPool<E> {
    fn capacity(&self) -> usize {
        self.vec.capacity()
    }

    fn len(&self) -> usize {
        self.vec.len()
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.vec.swap(i, j);
    }
}




pub struct BVH<T, E, NodePool, ElementPool, const DIM: usize>
where
    E: BVHElement<T, DIM>,
    NodePool: BVHPool<T, DIM>,
    ElementPool: BVHElementPool<T, E, DIM>, {

    pub pool: NodePool,
    pub elements: ElementPool,
    root: usize,
    nodes_in_use: usize,


    _t: PhantomData<T>,
    _e: PhantomData<E>,
}

impl<T, E, ElementPool, const DIM: usize> BVH<T, E, VecPool<BVHNode<T, DIM>>, ElementPool, DIM>
where T: BaseFloat + From<u32>,
      E: BVHElement<T, DIM>,
      ElementPool: BVHElementPool<T, E, DIM> {

    /// Creates a new BVH from the specified element pool.
    ///
    /// This function will only construct the basic data structure of the BVH. It will not attempt
    /// to construct it. A BVH-tree constructed from this function may be build using
    /// ``
    /// let mut bvh = BVH::new(elements);
    /// bvh.rebuild<BVHSplitting>();
    /// ``
    pub fn new(elements: ElementPool) -> Self {
        let mut pool = VecPool::with_capacity(elements.capacity() * 2 - 1);
        for _ in 0..pool.vec.capacity() {
            pool.push(BVHNode::new());
        }

        BVH {
            pool,
            elements,
            root: 0,
            nodes_in_use: 1,

            _t: PhantomData::default(),
            _e: PhantomData::default(),
        }
    }
}

impl<T, E, NodePool, ElementPool, const DIM: usize> BVH<T, E, NodePool, ElementPool, DIM>
where T: BaseFloat + From<u32>,
      E: BVHElement<T, DIM>,
      NodePool: BVHPool<T, DIM>,
      ElementPool: BVHElementPool<T, E, DIM> {

    /// Rebuilds the BVH-tree using the specified splitting function `SF`.
    pub fn rebuild<SF: BVHSplitting<T, E, NodePool, ElementPool, DIM>>(&mut self) {
        self.nodes_in_use = 1;
        let root = &mut self.pool[self.root];
        root.left_first = 0;
        root.num_prims = self.elements.len();

        self.update_bounds(self.root);
        self.subdivide::<SF>(self.root);
    }

    /// Refits the BVH-tree to the current state of the tree nodes.
    pub fn refit(&mut self) {
        for i in 0..self.nodes_in_use {
            let id = self.nodes_in_use - i - 1;
            let node = &self.pool[id];
            if node.is_leaf() {
                self.update_bounds(id);
            } else {
                let left_child = self.pool[node.left_first].aabb.clone();
                let right_child = self.pool[node.right_child()].aabb.clone();

                self.pool[id].aabb.adjust(&left_child, &right_child);
            }
        }
    }

    /// Updates the bounds for the node with the specified `node_id`.
    pub fn update_bounds(&mut self, node_id: usize) {
        let node = &mut self.pool[node_id];
        node.aabb.reset();

        let first = node.left_first;
        for i in 0..node.num_prims {
            let element = &self.elements[first + i];
            node.aabb.grow_other(&element.wrap());
        }
    }

    /// Subdivides the node specified by `node_id` by using the specified splitting function.
    pub fn subdivide<SF: BVHSplitting<T, E, NodePool, ElementPool, DIM>>(
        &mut self, node_id: usize
    ) {
        let node = &self.pool[node_id];

        // split plane axis and position
        let split = SF::find(self, node);
        if split.cost >= Self::calc_node_cost(node) {
            return; // not splitting is more cost-effective
        }

        // split the group in two halves
        let mut i = node.left_first;
        let mut j = i + node.num_prims - 1;
        while i <= j {
            if self.elements[i].centroid()[split.axis] < split.pos {
                // element is to the left of the split
                i += 1;
            } else {
                // element is to the right of the split
                self.elements.swap(i, j);
                j -= 1;
            }
        }

        // create child nodes for each half
        let left_count = i - node.left_first;
        if left_count == 0 || left_count == node.num_prims {
            return;
        }

        let left_child_idx = self.nodes_in_use;
        self.nodes_in_use += 1;
        let right_child_idx = self.nodes_in_use;
        self.nodes_in_use += 1;

        let left_first = node.left_first;
        let num_prims = node.num_prims;

        let left_child = &mut self.pool[left_child_idx];
        left_child.left_first = left_first;
        left_child.num_prims = left_count;
        let right_child = &mut self.pool[right_child_idx];
        right_child.left_first = i;
        right_child.num_prims = num_prims - left_count;

        let node = &mut self.pool[node_id];
        node.num_prims = 0;
        node.left_first = left_child_idx;


        // update child bounds
        self.update_bounds(left_child_idx);
        self.update_bounds(right_child_idx);
        // try to recursively subdivide the children
        self.subdivide::<SF>(left_child_idx);
        self.subdivide::<SF>(right_child_idx);
    }

    /// Returns the SAH evaluation for the specified `node` with the specified splitting `pos` along
    /// the specified splitting `axis`. The return value of this method by be used as an
    /// approximation for the cost of splitting the node at the specified split when traversing the
    /// tree during an intersection search. This is used by the different splitting functions.
    pub fn eval_sah(&self, node: &BVHNode<T, DIM>, axis: usize, pos: T) -> T {
        // determine element counts and bounds for this split candidate
        let mut leftbox = AABB::<T, DIM>::new();
        let mut rightbox = AABB::<T, DIM>::new();
        let mut left_count = 0usize;
        let mut right_count = 0usize;
        for i in 0..node.num_prims {
            let element = &self.elements[node.left_first + i];
            if element.centroid()[axis] < pos {
                left_count += 1;
                leftbox.grow_other(&element.wrap());
            } else {
                right_count += 1;
                rightbox.grow_other(&element.wrap());
            }
        }
        let cost = T::from(left_count as u32) * leftbox.area() + T::from(right_count as u32) * rightbox.area();
        if cost > T::zero() {
            cost
        } else {
            T::MAX
        }
    }

    /// Returns a cost approximation for searching the specified node.
    fn calc_node_cost(node: &BVHNode<T, DIM>) -> T {
        T::from(node.num_prims as u32) * node.aabb.area()
    }

    /// Returns a `Vec` to references of the member elements of this tree that intersect the
    /// specified intersector. Since intersection tests from the side of the tree are done in the
    /// BVH's frame of reference, the `intersector` instance should be transformed into the
    /// reference frame of the BVH *before* this method is called.
    pub fn intersect<I: BVIntersector<T, E, DIM> + BVIntersector<T, AABB<T, DIM>, DIM>>(
        &self, intersector: &I, node_idx: usize) -> Vec<&E> {

        let mut v = Vec::<&E>::with_capacity(64);

        let mut node = &self.pool[node_idx];
        let mut stack = [node; 64];
        let mut stack_ptr = 0usize;

        loop {
            if node.is_leaf() {
                for i in 0..node.num_prims {
                    if intersector.intersects(&self.elements[node.left_first + i]) {
                        v.push(&self.elements[node.left_first + i]);
                    }
                }

                if stack_ptr == 0 {
                    break;
                } else {
                    stack_ptr -= 1;
                    node = stack[stack_ptr];
                }
            } else {
                let mut child1 = &self.pool[node.left_first];
                let mut child2 = &self.pool[node.right_child()];

                let mut inter1 = intersector.intersects(&child1.aabb);
                let mut inter2 = intersector.intersects(&child2.aabb);
                // -(for ray intersections, do ray sorting here)-
                if !inter1 {
                    // if child 1 does not intersect the intersector, swap with child 2
                    mem::swap(&mut child1, &mut child2);
                    mem::swap(&mut inter1, &mut inter2);
                }

                if !inter1 {
                    // both children do not intersect the intersector. Checkout stack
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




#[cfg(test)]
mod test {
    use nalgebra::SVector;
    use crate::volume::aabb::AABB;
    use crate::volume::{BoundingVolume, bvh_splitting};
    use crate::volume::bvh::{BVH, BVHElement, BVHNode, VecPool};

    struct Test<const DIM: usize> {
        bounds: AABB<f64, DIM>
    }

    impl<const DIM: usize> BoundingVolume<f64, DIM> for Test<DIM> {
        fn center(&self) -> SVector<f64, DIM> {
            self.bounds.center()
        }

        fn area(&self) -> f64 {
            self.bounds.area()
        }

        fn min(&self) -> SVector<f64, DIM> {
            self.bounds.min.clone()
        }

        fn max(&self) -> SVector<f64, DIM> {
            self.bounds.max.clone()
        }

        fn size(&self) -> SVector<f64, DIM> {
            self.bounds.size()
        }

        fn half_size(&self) -> SVector<f64, DIM> {
            self.bounds.half_size()
        }
    }

    impl<const DIM: usize> BVHElement<f64, DIM> for Test<DIM> {
        fn centroid(&self) -> SVector<f64, DIM> {
            self.bounds.center()
        }

        fn wrap(&self) -> AABB<f64, DIM> {
            self.bounds.clone()
        }
    }

    #[test]
    fn test() {
        let mut elements = VecPool::<Test<2>>::with_capacity(10);

        let mut bvh = BVH::<f64, Test<2>, VecPool<BVHNode<f64, 2>>, VecPool<Test<2>>, 2>::new(elements);
        bvh.rebuild::<bvh_splitting::BinnedSAHSplit<8>>();
    }
}

