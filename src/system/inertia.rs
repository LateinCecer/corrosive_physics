use std::mem;
use std::ops::{AddAssign, Neg, SubAssign};
use nalgebra::{ClosedAdd, ClosedMul, ComplexField, Const, DefaultAllocator, Dim, Matrix, Matrix3, Matrix4, OMatrix, RealField, Scalar, Storage, UnitQuaternion, Vector3};
use nalgebra::allocator::Allocator;
use num::{One, Zero};
use crate::helper::{BaseFloat, mat};

/// The base error type for the error classes used by the physics engine core.
pub enum ErrorType {
    /// The math error enum type is used for all algebraic errors, like, for example, when
    /// dividing by a zero-value or trying to invert a non-invertible matrix.
    MathError,
    /// The physics error enum type is used for all errors, that originate from errors in the
    /// physical model of the world. For example, the inertia tensor for a 3d-object should always
    /// be a 3x3 invertible matrix.
    PhysicsError,
}

/// Base error structure. An error consists of an error base type and an optional error message.
/// To generate an error, the build-in `err!` macro should be used.
pub struct Error {
    msg: Option<String>,
    ty: ErrorType,
}

impl Error {
    /// Constructor for a n error type.
    pub fn new(ty: ErrorType, msg: Option<String>) -> Self {
        Error {
            msg,
            ty
        }
    }
}

macro_rules! err {
    (math) => (
        Error::new(ErrorType::MathError, None)
    );
    (math $msg:expr) => (
        Error::new(ErrorType::MathError, Some(String::from($msg)))
    );
    (physics) => (
        Error::new(ErrorType::PhysicsError, None)
    );
    (physics $msg:expr) => (
        Error::new(ErrorType::PhysicsError, Some(String::from($msg)))
    );
}
pub(crate) use err;



/// Data structure for a physical inertial system.
///
/// A physical inertial system is defined by the momentum (translational), angular momentum
/// (rotational), the transformer state (position, scale, rotation, ...) and the mass distribution
/// of the system.
/// This structure may be packaged into component data structures together with children objects,
/// mesh-data, and other components.
#[derive(Clone, Debug)]
pub struct IS<T> {
    pub momentum: Vector3<T>,
    pub angular_mom: Vector3<T>,
    pub state: Transformer<T>,
    pub mass: MassDistribution<T>,
}

/// Data structure for the mass distributions of an inertial system.
#[derive(Clone, Debug)]
pub struct MassDistribution<T> {
    mass: T,
    center_of_mass: Vector3<T>,
    inertia: Matrix3<T>,
    inv_inertia: Matrix3<T>,
}

/// Data structure for a transformer state.
#[derive(Clone, Debug)]
pub struct Transformer<T> {
    pub pos: Vector3<T>,
    pub offset: Vector3<T>,
    pub scale: Vector3<T>,
    pub rot: UnitQuaternion<T>,

    /// Transformation matrix for transforming points and vectors into the laboratory frame
    mat: Matrix4<T>,
    /// Transformation matrix for transforming points and vectors into the inertial reference frame
    inv_mat: Matrix4<T>,
}


pub trait Inertia<T>
where T: Scalar + Copy + ClosedMul<T> + ClosedAdd<T> + AddAssign<T> + Neg<Output=T> {
    /// Adds a mass point to the inertia system. The mass point is specified by a point vector `r`
    /// and a scalar `mass`.
    fn add_mass_point(&mut self, r: &Vector3<T>, mass: T);
    /// Subs a mass point to the inertia system. The mass point is specified by a point vector `r`
    /// and a scalar `mass`.
    fn sub_mass_point(&mut self, r: &Vector3<T>, mass: T);
}

macro_rules! assignop_inertia {
    ($in:expr, $mass:expr, $r:expr, $op:ident) => {{
        unsafe {
            $in.get_unchecked_mut((0,0)).$op($mass * ($r[1] * $r[1] + $r[2] * $r[2]))
        };
        unsafe {
            $in.get_unchecked_mut((1,0)).$op(- $mass * $r[1] * $r[0])
        };
        unsafe {
            $in.get_unchecked_mut((2,0)).$op(- $mass * $r[2] * $r[0])
        };

        unsafe {
            *$in.get_unchecked_mut((0,1)) = $in.get_unchecked((1,0)).clone()
        };
        unsafe {
            $in.get_unchecked_mut((1,1)).$op($mass * ($r[0] * $r[0] + $r[2] * $r[2]))
        };
        unsafe {
            $in.get_unchecked_mut((2,1)).$op(- $mass * $r[2] * $r[1])
        };

        unsafe {
            *$in.get_unchecked_mut((0,2)) = $in.get_unchecked((2,0)).clone()
        };
        unsafe {
            *$in.get_unchecked_mut((1,2)) = $in.get_unchecked((2,1)).clone()
        };
        unsafe {
            $in.get_unchecked_mut((2,2)).$op($mass * ($r[0] * $r[0] + $r[1] * $r[1]))
        };
    }}
}

impl<T> Inertia<T> for Matrix3<T>
where T: Scalar + Copy + ClosedMul<T> + ClosedAdd<T> + AddAssign<T> + SubAssign<T> + Neg<Output=T> {
    fn add_mass_point(&mut self, r: &Vector3<T>, mass: T) {
        assignop_inertia!(self, mass, r, add_assign)
    }

    fn sub_mass_point(&mut self, r: &Vector3<T>, mass: T) {
        assignop_inertia!(self, mass, r, sub_assign)
    }
}


impl<T> IS<T> {
    /// Constructor for an inertial system.
    pub fn new(
        mom: Vector3<T>,
        angular_mom: Vector3<T>,
        state: Transformer<T>,
        mass: MassDistribution<T>,
    ) -> Self {
        IS {
            momentum: mom,
            angular_mom,
            state,
            mass,
        }
    }
}

impl<T> IS<T>
where T: BaseFloat {

    /// Returns the velocity of a single point within the inertial system. The specified point
    /// position and the velocity are specified as within the reference frame of this inertial
    /// system.
    ///
    /// To get the point velocity from outside of the inertial system, all values have to be
    /// transformed. This could look something like this:
    ///
    /// ``
    /// is.trafo_outof(
    ///     &is.get_point_vel(
    ///         &is.trafo_into(&point)
    ///     )
    /// )
    /// ``
    pub fn get_point_vel(&self, point: &Vector3<T>) -> Vector3<T> {
        self.get_angular_vel().cross(point)
    }

    /// Returns the angular velocity of the inertial system within the reference frame of the
    /// inertial system.
    pub fn get_angular_vel(&self) -> Vector3<T> {
        self.mass.inv_inertia * self.angular_mom
    }

    /// Applies an impulse to a specified point of the inertial system. All values are to be
    /// provided from the reference frame of the inertial system.
    pub fn apply_impulse(&mut self, imp: &Vector3<T>, point: &Vector3<T>) {
        self.momentum += imp;
        self.angular_mom += point.cross(imp);
    }

    pub fn integrate(&mut self, t: T) {
        self.state.pos += self.momentum.scale(t / self.mass.mass);
        let rot = UnitQuaternion::new(self.get_angular_vel().scale(t));
        self.state.rot = rot * self.state.rot;
    }

    pub fn sync(&mut self) {
        self.state.update_transformation();
    }

    /// Transforms a matrix value from the laboratory frame into the reference frame of the
    /// inertial system.
    pub fn trafo_into<C, ST>(&self, vec: &Matrix<T, Const<4>, C, ST>) -> OMatrix<T, Const<4>, C>
    where C: Dim,
          ST: Storage<T, Const<4>, C>,
          DefaultAllocator: Allocator<T, Const<4>, C> {
        self.state.inv_mat * vec
    }

    /// Transforms a matrix value from the reference frame of the inertial system into the
    /// laboratory frame.
    pub fn trafo_outof<C, ST>(&self, vec: &Matrix<T, Const<4>, C, ST>) -> OMatrix<T, Const<4>, C>
    where C: Dim,
          ST: Storage<T, Const<4>, C>,
          DefaultAllocator: Allocator<T, Const<4>, C> {
        self.state.mat * vec
    }

    /// Transforms a 3d-vector value from the laboratory frame into the reference frame of the
    /// inertial system.
    pub fn trafo_vec_into(&self, vec: &Vector3<T>) -> Vector3<T> {
        self.state.inv_trafo_vec(vec)
    }

    /// Transforms a 3d-point value from the laboratory frame into the reference frame of the
    /// inertial system.
    pub fn trafo_point_into(&self, point: &Vector3<T>) -> Vector3<T> {
        self.state.inv_trafo_point(point)
    }

    /// Transforms a unit quaternion from the laboratory frame into the reference frame of the
    /// inertial system.
    pub fn trafo_rot_into(&self, rot: &UnitQuaternion<T>) -> UnitQuaternion<T> {
        self.state.inv_trafo_rot(rot)
    }

    /// Transforms a transformer state from the laboratory frame into the reference frame of the
    /// inertial system.
    pub fn trafo_state_into(&self, state: &Transformer<T>) -> Transformer<T> {
        self.state.inv_trafo(state)
    }

    /// Mutably transforms a transformer state from the laboratory frame into the reference frame of
    /// the inertial system.
    pub fn trafo_state_into_mut(&self, state: &mut Transformer<T>) {
        self.state.inv_trafo_mut(state)
    }

    /// Transforms a 3d-vector value from the reference frame of the inertial system into the
    /// laboratory frame.
    pub fn trafo_vec_outof(&self, vec: &Vector3<T>) -> Vector3<T> {
        self.state.trafo_vec(vec)
    }

    /// Transforms a 3d-point value from the reference frame of the inertial system into the
    /// laboratory frame.
    pub fn trafo_point_outof(&self, point: &Vector3<T>) -> Vector3<T> {
        self.state.trafo_point(point)
    }

    /// Transforms a unit quaternion from the reference frame of the inertial system into the
    /// laboratory frame.
    pub fn trafo_rot_outof(&self, rot: &UnitQuaternion<T>) -> UnitQuaternion<T> {
        self.state.trafo_rot(rot)
    }

    /// Transforms a transformer state from the reference frame of the inertial system into the
    /// laboratory frame.
    pub fn trafo_state_outof(&self, state: &Transformer<T>) -> Transformer<T> {
        self.state.trafo(state)
    }

    /// Mutably transforms a transformer state from the reference frame of the inertial system into
    /// the laboratory frame.
    pub fn trafo_state_outof_mut(&self, state: &mut Transformer<T>) {
        self.state.trafo_mut(state)
    }
}





impl<T> Default for MassDistribution<T>
where T: Scalar + Zero + One {
    fn default() -> Self {
        MassDistribution {
            mass: T::one(),
            center_of_mass: Vector3::zeros(),
            inertia: Matrix3::identity(),
            inv_inertia: Matrix3::identity(),
        }
    }
}

impl<T> MassDistribution<T>
where T: Scalar + ComplexField {
    /// Builds a new mass distribution structure from the point mass `mass`, the center of mass
    /// `com` and the inertia tensor `inertia` of the inertial system in question. Since the
    /// inertia tensor is stored in both the default and the inverted state, this method will
    /// attempt to attain the inverse of the inertia matrix. If this fails, this means that
    ///
    /// - a: the system is unphysical and that
    /// - b: the mass distribution cannot be generated.
    ///
    /// In this case, an error type is returned.
    pub fn new(mass: T, com: Vector3<T>, inertia: Matrix3<T>) -> Result<Self, Error> {
        // try to invert inertia and build mass distribution from there
        Ok(MassDistribution {
            mass,
            center_of_mass: com,
            inv_inertia: inertia.clone().try_inverse()
                .ok_or(err!(physics "Failed to invert inertia tensor"))?,
            inertia,
        })
    }
}

impl<T> MassDistribution<T> {
    /// Returns the total mass of the mass distribution.
    pub fn mass(&self) -> &T {
        &self.mass
    }

    /// Returns the center of mass of the mass distribution.
    ///
    /// In an inertia system (`IS`) the center of mass will usually be used as the `offset` in the
    /// transformation state of the system. This ensures, that a free object rotates around it's
    /// center of mass.
    /// A non-free object on the other hand, is forced to rotate around an other reference point.
    /// In this case, the offset of the translation state is this rotational reference point.
    pub fn center_of_mass(&self) -> &Vector3<T> {
        &self.center_of_mass
    }

    /// Returns the inertia tensor for the mass distribution.
    ///
    /// Since the inertia tensor is usually defined within the moving reference frame, the tensor
    /// does not change when the transformation of the reference system changes.
    pub fn inertia(&self) -> &Matrix3<T> {
        &self.inertia
    }

    /// Returns the inverse of the inertia tensor.
    ///
    /// Since the inertia tensor only changes when the mass distribution of an object is changed,
    /// it makes sense to cache the inverted inertia tensor to save on the computational cost of
    /// inverting the matrix every time the inverse is needed.
    pub fn inv_inertia(&self) -> &Matrix3<T> {
        &self.inv_inertia
    }
}





impl<T> Default for Transformer<T>
where T: Scalar + Zero + One + RealField {
    fn default() -> Self {
        Transformer {
            mat: Matrix4::identity(),
            inv_mat: Matrix4::identity(),
            pos: Vector3::zeros(),
            offset: Vector3::zeros(),
            rot: UnitQuaternion::identity(),
            scale: Vector3::repeat(T::one())
        }
    }
}

impl<T> Transformer<T>
where T: BaseFloat {

    pub fn new(pos: Vector3<T>, rot: UnitQuaternion<T>, scale: Vector3<T>, offset: Vector3<T>) -> Self {
        Transformer {
            mat: Self::gen_mat(&pos, &rot, &scale, &offset),
            inv_mat: Self::gen_inv_mat(&pos, &rot, &scale, &offset),
            pos,
            rot,
            scale,
            offset,
        }
    }

    /// Updates the transformation matrices of this transformer.
    pub fn update_transformation(&mut self) {
        self.mat = Self::gen_mat(&self.pos, &self.rot, &self.scale, &self.offset);
        self.inv_mat = Self::gen_inv_mat(&self.pos, &self.rot, &self.scale, &self.offset);
    }

    /// Generates a transformation matrix for the specified transformer state.
    fn gen_mat(pos: &Vector3<T>, rot: &UnitQuaternion<T>, scale: &Vector3<T>, offset: &Vector3<T>) -> Matrix4<T> {
        mat::init_translation(pos)
            * mat::init_rotation(rot)
            * mat::init_scale(scale)
            * mat::init_translation(offset)
    }

    /// Generates an inverse transformation matrix for the specified transformation matrix.
    fn gen_inv_mat(pos: &Vector3<T>, rot: &UnitQuaternion<T>, scale: &Vector3<T>, offset: &Vector3<T>) -> Matrix4<T> {
        mat::init_inverse_translation(offset)
            * mat::init_inverse_scale(scale)
            * mat::init_inverse_rotation(rot)
            * mat::init_inverse_translation(pos)
    }

    /// Returns the transformation matrix for this transformer.
    pub fn tsro(&self) -> &Matrix4<T> {
        &self.mat
    }

    /// Returns the inverse transformation matrix for this transformer.
    pub fn inv_tsro(&self) -> &Matrix4<T> {
        &self.inv_mat
    }
}

macro_rules! mat_vec_mul_row {
    ($mat:expr, point $point:expr, ($row:tt)) => (
        unsafe {
            $mat.get_unchecked(($row,0)).clone() * $point.get_unchecked(0).clone()
                + $mat.get_unchecked(($row,1)).clone() * $point.get_unchecked(1).clone()
                + $mat.get_unchecked(($row,2)).clone() * $point.get_unchecked(2).clone()
                + $mat.get_unchecked(($row,3)).clone()
        }
    );
    ($mat:expr, vec $vec:expr, ($row:tt)) => (
        unsafe {
            $mat.get_unchecked(($row,0)).clone() * $vec.get_unchecked(0).clone()
                + $mat.get_unchecked(($row,1)).clone() * $vec.get_unchecked(1).clone()
                + $mat.get_unchecked(($row,2)).clone() * $vec.get_unchecked(2).clone()
        }
    );
}

impl<T> Transformer<T>
where T: BaseFloat {
    pub fn trafo_point(&self, point: &Vector3<T>) -> Vector3<T> {
        Vector3::new(
            mat_vec_mul_row!(self.mat, point point, (0)),
            mat_vec_mul_row!(self.mat, point point, (1)),
            mat_vec_mul_row!(self.mat, point point, (2)),
        )
    }

    pub fn trafo_vec(&self, vec: &Vector3<T>) -> Vector3<T> {
        Vector3::new(
            mat_vec_mul_row!(self.mat, vec vec, (0)),
            mat_vec_mul_row!(self.mat, vec vec, (1)),
            mat_vec_mul_row!(self.mat, vec vec, (2)),
        )
    }

    pub fn inv_trafo_point(&self, point: &Vector3<T>) -> Vector3<T> {
        Vector3::new(
            mat_vec_mul_row!(self.inv_mat, point point, (0)),
            mat_vec_mul_row!(self.inv_mat, point point, (1)),
            mat_vec_mul_row!(self.inv_mat, point point, (2)),
        )
    }

    pub fn inv_trafo_vec(&self, vec: &Vector3<T>) -> Vector3<T> {
        Vector3::new(
            mat_vec_mul_row!(self.inv_mat, vec vec, (0)),
            mat_vec_mul_row!(self.inv_mat, vec vec, (1)),
            mat_vec_mul_row!(self.inv_mat, vec vec, (2)),
        )
    }

    pub fn trafo_rot(&self, rot: &UnitQuaternion<T>) -> UnitQuaternion<T> {
        self.rot * rot
    }

    pub fn inv_trafo_rot(&self, rot: &UnitQuaternion<T>) -> UnitQuaternion<T> {
        self.rot.conjugate() * rot
    }

    pub fn trafo(&self, trafo: &Transformer<T>) -> Transformer<T> {
        Transformer {
            pos: self.trafo_point(&trafo.pos),
            offset: self.trafo_vec(&trafo.offset),
            rot: self.trafo_rot(&trafo.rot),
            scale: self.scale.component_mul(&trafo.scale),

            mat: self.mat * trafo.mat,
            inv_mat: trafo.inv_mat * self.inv_mat,
        }
    }

    pub fn trafo_mut(&self, trafo: &mut Transformer<T>) {
        trafo.pos = self.trafo_point(&trafo.pos);
        trafo.offset = self.trafo_vec(&trafo.offset);
        trafo.rot = self.trafo_rot(&trafo.rot);
        trafo.scale.component_mul_assign(&self.scale);
        trafo.mat = self.mat * trafo.mat;
        trafo.inv_mat = trafo.inv_mat * self.inv_mat;
    }

    pub fn inv_trafo(&self, trafo: &Transformer<T>) -> Transformer<T> {
        Transformer {
            pos: self.inv_trafo_point(&trafo.pos),
            offset: self.inv_trafo_vec(&trafo.offset),
            rot: self.inv_trafo_rot(&trafo.rot),
            scale: trafo.scale.component_div(&self.scale),

            mat: self.inv_mat * trafo.mat,
            inv_mat: trafo.inv_mat * self.mat,
        }
    }

    pub fn inv_trafo_mut(&self, trafo: &mut Transformer<T>) {
        trafo.pos = self.inv_trafo_point(&trafo.pos);
        trafo.offset = self.inv_trafo_vec(&trafo.offset);
        trafo.rot = self.inv_trafo_rot(&trafo.rot);
        trafo.scale.component_div_assign(&self.scale);
        trafo.mat = self.inv_mat * trafo.mat;
        trafo.inv_mat = trafo.inv_mat * self.mat;
    }

    /// Generates an inverted copy of the transformation state.
    pub fn inverse(&self) -> Transformer<T> {
        Transformer {
            pos: -self.pos,
            offset: -self.offset,
            rot: self.rot.conjugate(),
            scale: Vector3::repeat(T::one()).component_div(&self.scale),
            mat: self.inv_mat,
            inv_mat: self.mat,
        }
    }

    /// Inverts the current transformation state instance.
    pub fn inverse_mut(&mut self) {
        self.pos = -self.pos;
        self.offset = -self.offset;
        self.rot.conjugate_mut();
        self.scale = Vector3::repeat(T::one()).component_div(&self.scale);
        mem::swap(&mut self.inv_mat, &mut self.mat);
    }

    /// Returns the vector pointing to the 'right' in the laboratory frame for the current transformer
    /// state. In a right-handed euclidean coordinate system, the 'right' is defined as the unit
    /// vector pointing in _positive x_ direction.
    pub fn right(&self) -> Vector3<T> {
        mat::right(&self.rot)
    }

    /// Returns the vector pointing to the 'left' in the laboratory frame for the current transformer
    /// state. In a right-handing euclidean coordinate system, the 'left' is defined as the unit
    /// vector pointing in _negative x_ direction.
    pub fn left(&self) -> Vector3<T> {
        -mat::right(&self.rot)
    }

    /// Returns the vector pointing 'upwards' in the laboratory frame for the current transformer
    /// state. In a right-handed euclidean coordinate system, 'upwards' is defined as the unit
    /// vector pointing in _positive y_ direction.
    pub fn up(&self) -> Vector3<T> {
        mat::up(&self.rot)
    }

    /// Returns the vector pointing 'downwards' in the laboratory frame for the current transformer
    /// state. In a right-handed euclidean coordinate system, 'downwards' is defined as the unit
    /// vector pointing in _negative y_ direction.
    pub fn down(&self) -> Vector3<T> {
        -mat::up(&self.rot)
    }

    /// Returns the vector pointing 'forward' in the laboratory frame for the current transformer
    /// state. In a right-handed euclidean coordinate system, 'forwards' is defined as the unit
    /// vector pointing in _positive z_ direction.
    pub fn forward(&self) -> Vector3<T> {
        mat::forward(&self.rot)
    }

    /// Returns the vector pointing 'backwards' in the laboratory frame for the current transformer
    /// state. In a right-handed euclidean coordinate system, 'backwards' is defined as the unit
    /// vector pointing in _negative z_ direction.
    pub fn backward(&self) -> Vector3<T> {
        -mat::forward(&self.rot)
    }
}
