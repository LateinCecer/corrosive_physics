use crate::helper::BaseFloat;
use crate::volume::aabb::AABB;
use crate::volume::BoundingVolume;
use crate::volume::bvh::{BVH, BVHElement, BVHElementPool, BVHNode, BVHPool};

pub struct BVHSplit<T> {
    pub cost: T,
    pub pos: T,
    pub axis: usize,
}

pub trait BVHSplitting<T, E, NPool, EPool, const DIM: usize>
where E: BVHElement<T, DIM>,
      NPool: BVHPool<T, DIM>,
      EPool: BVHElementPool<T, E, DIM> {

    /// Finds a split for the specified BVH-tree and -node.
    fn find(
        bvh: &BVH<T, E, NPool, EPool, DIM>,
        node: &BVHNode<T, DIM>
    ) -> BVHSplit<T>;
}




pub struct FullSAHSplit {}
impl<T: BaseFloat + From<u32>, E, NPool, EPool, const DIM: usize> BVHSplitting<T, E, NPool, EPool, DIM>
for FullSAHSplit
where E: BVHElement<T, DIM>,
      NPool: BVHPool<T, DIM>,
      EPool: BVHElementPool<T, E, DIM> {

    fn find(bvh: &BVH<T, E, NPool, EPool, DIM>, node: &BVHNode<T, DIM>) -> BVHSplit<T> {
        // split plane axis and position
        let mut split_pos = T::zero();
        let mut best_axis = 0usize;

        // determine split axis using SAH
        let mut best_cost = T::MAX;
        for i in 0..*node.num_prims() {
            let element = &bvh.elements[node.left_child() + i];
            for axis in 0..DIM {
                // println!("      searching axis {axis}");
                let candidate_pos = element.centroid()[axis];
                let cost = bvh.eval_sah(node, axis, candidate_pos);
                if cost < best_cost {
                    split_pos = candidate_pos;
                    best_axis = axis;
                    best_cost = cost;
                }
            }
        }

        BVHSplit {
            cost: best_cost,
            pos: split_pos,
            axis: best_axis,
        }
    }
}


macro_rules! axis_min_max {
    ($T:ty, $bvh:expr, $node:expr, $axis:expr) => {{
        let mut bounds_min = <$T>::MAX;
        let mut bounds_max = <$T>::MIN;
        for i in 0..*$node.num_prims() {
            let centroid = $bvh.elements[$node.left_child() + i].centroid();
            bounds_min = <$T>::min(bounds_min, centroid[$axis]);
            bounds_max = <$T>::max(bounds_max, centroid[$axis]);
        }
        (bounds_min, bounds_max)
    }}
}


pub struct MidpointSAHSplit {}
impl<T: BaseFloat + From<u32>, E, NPool, EPool, const DIM: usize> BVHSplitting<T, E, NPool, EPool, DIM>
for MidpointSAHSplit
where E: BVHElement<T, DIM>,
      NPool: BVHPool<T, DIM>,
      EPool: BVHElementPool<T, E, DIM> {

    fn find(bvh: &BVH<T, E, NPool, EPool, DIM>, node: &BVHNode<T, DIM>) -> BVHSplit<T> {
        let mut best_cost = T::MAX;
        let mut split_pos = T::zero();
        let mut best_axis = 0usize;

        // try every axis
        for axis in 0..DIM {
            let (bounds_min, bounds_max) = axis_min_max!(T, bvh, node, axis);

            if bounds_min != bounds_max {
                let candidate_pos = (bounds_max + bounds_min) * T::half();
                let cost = bvh.eval_sah(node, axis, candidate_pos);
                if cost < best_cost {
                    split_pos = candidate_pos;
                    best_axis = axis;
                    best_cost = cost;
                }
            }
        }

        BVHSplit {
            cost: best_cost,
            pos: split_pos,
            axis: best_axis,
        }
    }
}



pub struct PartialSAHSplit<const NUM_PLANES: usize> {}
impl<T: BaseFloat + From<u32>, E, NPool, EPool, const NUM_PLANES: usize, const DIM: usize>
BVHSplitting<T, E, NPool, EPool, DIM>
for PartialSAHSplit<NUM_PLANES>
where E: BVHElement<T, DIM>,
      NPool: BVHPool<T, DIM>,
      EPool: BVHElementPool<T, E, DIM> {

    fn find(bvh: &BVH<T, E, NPool, EPool, DIM>, node: &BVHNode<T, DIM>) -> BVHSplit<T> {
        let r_num_planes = T::one() / T::from(NUM_PLANES as u32);

        let mut best_cost = T::MAX;
        let mut split_pos = T::zero();
        let mut best_axis = 0usize;


        // loop through axis
        for axis in 0..DIM {
            let (bounds_min, bounds_max) = axis_min_max!(T, bvh, node, axis);

            if bounds_min != bounds_max {
                let scale = (bounds_max - bounds_min) * r_num_planes;
                for i in 1..NUM_PLANES {
                    let candidate_pos = bounds_min + T::from(i as u32) * scale;
                    let cost = bvh.eval_sah(node, axis, candidate_pos);
                    if cost < best_cost {
                        split_pos = candidate_pos;
                        best_axis = axis;
                        best_cost = cost;
                    }
                }
            }
        }

        BVHSplit {
            pos: split_pos,
            axis: best_axis,
            cost: best_cost,
        }
    }
}




#[derive(Clone, Copy)]
struct Bin<T: BaseFloat, const DIM: usize> {
    aabb: AABB<T, DIM>,
    prime_count: usize,
}
impl<T: BaseFloat, const DIM: usize> Bin<T, DIM> {
    pub fn zero() -> Self {
        Bin {
            aabb: AABB::new(),
            prime_count: 0
        }
    }

    pub fn reset(&mut self) {
        self.aabb.reset();
        self.prime_count = 0;
    }
}

pub struct BinnedSAHSplit<const NUM_BINS: usize> {}

impl<T: BaseFloat + From<u32>, E, NPool, EPool, const NUM_BINS: usize, const DIM: usize>
BVHSplitting<T, E, NPool, EPool, DIM>
for BinnedSAHSplit<NUM_BINS>
where E: BVHElement<T, DIM>,
      NPool: BVHPool<T, DIM>,
      EPool: BVHElementPool<T, E, DIM> {

    fn find(bvh: &BVH<T, E, NPool, EPool, DIM>, node: &BVHNode<T, DIM>) -> BVHSplit<T> {
        let r_num_bins = T::one() / T::from(NUM_BINS as u32);

        let mut best_cost = T::MAX;
        let mut split_pos = T::zero();
        let mut best_axis = 0usize;


        let mut bins = [Bin::<T, DIM>::zero(); NUM_BINS];
        let mut left_area = [T::zero(); NUM_BINS];
        let mut right_area = [T::zero(); NUM_BINS];
        let mut left_count = [0usize; NUM_BINS];
        let mut right_count = [0usize; NUM_BINS];
        let mut leftbox = AABB::<T, DIM>::new();
        let mut rightbox = AABB::<T, DIM>::new();

        for axis in 0..DIM {
            let (bounds_min, bounds_max) = axis_min_max!(T, bvh, node, axis);
            if bounds_min == bounds_max {
                continue;
            }


            // reset base
            bins.iter_mut().for_each(Bin::<T, DIM>::reset);
            // populate bins
            let mut scale = T::from(NUM_BINS as u32) / (bounds_max - bounds_min);
            for i in 0..*node.num_prims() {
                let element = &bvh.elements[node.left_child() + i];
                let bin_idx = usize::min(
                    NUM_BINS - 1,
                    T::floor_to_u32((element.centroid()[axis] - bounds_min) * scale) as usize);
                bins[bin_idx].prime_count += 1;
                bins[bin_idx].aabb.grow_other(&element.wrap());
            }

            // gather data for the `NUM_BINS - 1` planes between the `NUM_BINS` bins
            leftbox.reset();
            rightbox.reset();
            let mut left_sum = 0usize;
            let mut right_sum = 0usize;

            for i in 0..(NUM_BINS - 1) {
                left_sum += bins[i].prime_count;
                left_count[i] = left_sum;
                leftbox.grow_other(&bins[i].aabb);
                left_area[i] = leftbox.area();

                right_sum += bins[NUM_BINS - 1 - i].prime_count;
                right_count[NUM_BINS - 2 - i] = right_sum;
                rightbox.grow_other(&bins[NUM_BINS - 1 - i].aabb);
                right_area[NUM_BINS - 2 - i] = rightbox.area();
            }
            // calculate SAH cost for the planes
            scale = (bounds_max - bounds_min) * r_num_bins;
            for i in 0..(NUM_BINS - 1) {
                let plane_cost = T::from(left_count[i] as u32) * left_area[i]
                    + T::from(right_count[i] as u32) * right_area[i];

                if plane_cost < best_cost {
                    best_axis = axis;
                    split_pos = bounds_min + scale * (T::from(i as u32) + T::one());
                    best_cost = plane_cost;
                }
            }
        }

        BVHSplit {
            axis: best_axis,
            cost: best_cost,
            pos: split_pos,
        }
    }
}
