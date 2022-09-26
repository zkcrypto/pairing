use crate::ff::*;
use crate::*;
use std::fmt;

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait GenericCurveProjective:
    PartialEq
    + Eq
    + Sized
    + Copy
    + Clone
    + Send
    + Sync
    + fmt::Debug
    + fmt::Display
    + rand::Rand
    + 'static
{
    type Scalar: PrimeField;
    type Base: SqrtField;
    type Affine: GenericCurveAffine<Projective = Self, Scalar = Self::Scalar, Base = Self::Base>;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point is the point at infinity.
    fn is_zero(&self) -> bool;

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

    /// Checks if the point is already "normalized" so that
    /// cheap affine conversion is possible.
    fn is_normalized(&self) -> bool;

    /// Doubles this element.
    fn double(&mut self);

    /// Adds another element to this element.
    fn add_assign(&mut self, other: &Self);

    /// Subtracts another element from this element.
    fn sub_assign(&mut self, other: &Self) {
        let mut tmp = *other;
        tmp.negate();
        self.add_assign(&tmp);
    }

    /// Adds an affine element to this element.
    fn add_assign_mixed(&mut self, other: &Self::Affine);

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element.
    fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S);

    /// Converts this element into its affine representation.
    fn into_affine(&self) -> Self::Affine;

    /// Recommends a wNAF window table size given a scalar. Always returns a number
    /// between 2 and 22, inclusive.
    fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize;

    /// Recommends a wNAF window size given the number of scalars you intend to multiply
    /// a base by. Always returns a number between 2 and 22, inclusive.
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize;

    /// Returns references to underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn as_xyz(&self) -> (&Self::Base, &Self::Base, &Self::Base) {
        unimplemented!("default implementation does not exist for this function")
    }
    
    /// Returns underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn into_xyz_unchecked(self) -> (Self::Base, Self::Base, Self::Base) {
        unimplemented!("default implementation does not exist for this function")
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is NOT performed
    fn from_xyz_unchecked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Self {
        unimplemented!("default implementation does not exist for this function")
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is performed
    fn from_xyz_checked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Result<Self, GroupDecodingError> {
        unimplemented!("default implementation does not exist for this function")
    }
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait GenericCurveAffine:
    Copy + Clone + Sized + Send + Sync + fmt::Debug + fmt::Display + PartialEq + Eq + 'static
{
    type Scalar: PrimeField;
    type Base: SqrtField;
    type Projective: GenericCurveProjective<Affine = Self, Scalar = Self::Scalar, Base = Self::Base>;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point represents the point at infinity; the
    /// additive identity.
    fn is_zero(&self) -> bool;

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element with mixed addition.
    fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, other: S) -> Self::Projective;

    /// Converts this element into its affine representation.
    fn into_projective(&self) -> Self::Projective;

    /// Returns references to underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn as_xy(&self) -> (&Self::Base, &Self::Base);
    
    /// Returns underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn into_xy_unchecked(self) -> (Self::Base, Self::Base);

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is NOT performed
    fn from_xy_unchecked(x: Self::Base, y: Self::Base) -> Self;

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is performed
    fn from_xy_checked(x: Self::Base, y: Self::Base) -> Result<Self, GroupDecodingError>;

    /// returns A coefficient for a short Weierstrass form
    fn a_coeff() -> Self::Base;

    /// returns B coefficient for a short Weierstrass form
    fn b_coeff() -> Self::Base;
}
pub trait GenericUncompressedEncodable<const N: usize>: GenericCurveAffine {
    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity.
    fn into_uncompressed(&self) -> EncodingBytes<Self, N>;

    /// Converts an uncompressed encoding into the curve point
    fn from_uncompressed(encoding: EncodingBytes<Self, N>) -> Result<Self, GroupDecodingError>;
}
pub trait GenericCompressedEncodable<const N: usize>: GenericCurveAffine {
    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity.
    fn into_compressed(&self) -> (EncodingBytes<Self, N>, bool);

    /// Converts an uncompressed encoding into the curve point
    fn from_compressed(encoding: EncodingBytes<Self, N>, parity: bool) -> Result<Self, GroupDecodingError>;
}
pub trait GenericRawEncodable<const N: usize>: GenericUncompressedEncodable<N> {
    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity. Leaves coordinates in Montgommery form
    fn into_raw_uncompressed_le(&self) -> [u8; N];

    /// Creates a point from raw encoded coordinates without checking on curve
    fn from_raw_uncompressed_le_unchecked(encoded: &[u8; N], infinity: bool) -> Result<Self, GroupDecodingError>;

    /// Creates a point from raw encoded coordinates
    fn from_raw_uncompressed_le(encoded: &[u8; N], infinity: bool) -> Result<Self, GroupDecodingError>;
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct EncodingBytes<G: GenericCurveAffine, const N: usize>{
    bytes: [u8; N],
    _marker: std::marker::PhantomData<G>
}

impl<G: GenericCurveAffine, const N: usize> AsRef<[u8]> for EncodingBytes<G, N> {
    fn as_ref(&self) -> &[u8] {
        &self.bytes[..]
    }
}

impl<G: GenericCurveAffine, const N: usize> AsMut<[u8]> for EncodingBytes<G, N> {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.bytes[..]
    }
}

impl<G: GenericCurveAffine, const N: usize> EncodingBytes<G, N> {
    /// Creates an empty representation.
    pub fn empty() -> Self {
        Self {
            bytes: [0u8; N],
            _marker: std::marker::PhantomData
        }
    }

    /// Returns the number of bytes consumed by this representation.
    pub fn size() -> usize {
        N
    }

    /// Transforms into raw bytes without a type marker
    pub fn into_bytes(self) -> [u8; N] {
        self.bytes
    }

    /// Transforms from raw bytes by adding a type marker
    pub fn from_bytes(bytes: [u8; N]) -> Self {
        Self {
            bytes,
            _marker: std::marker::PhantomData
        }
    }
}

impl<G: CurveAffine> GenericCurveAffine for G {
    type Scalar = <Self as CurveAffine>::Scalar;
    type Base = <Self as CurveAffine>::Base;
    type Projective = <Self as CurveAffine>::Projective;

    /// Returns the additive identity.
    fn zero() -> Self {
        <Self as CurveAffine>::zero()
    }

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self {
        <Self as CurveAffine>::one()
    }

    /// Determines if this point is the point at infinity.
    fn is_zero(&self) -> bool {
        <Self as CurveAffine>::is_zero(&self)
    }

    /// Negates this element.
    fn negate(&mut self) {
        <Self as CurveAffine>::negate(self)
    }

    /// Performs scalar multiplication of this element with mixed addition.
    fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, other: S) -> Self::Projective {
        <Self as CurveAffine>::mul(self, other)
    }

    /// Converts this element into its affine representation.
    fn into_projective(&self) -> Self::Projective {
        <Self as CurveAffine>::into_projective(&self)
    }

    /// Returns references to underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn as_xy(&self) -> (&Self::Base, &Self::Base) {
        <Self as CurveAffine>::as_xy(&self)
    }
    
    /// Returns underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn into_xy_unchecked(self) -> (Self::Base, Self::Base) {
        <Self as CurveAffine>::into_xy_unchecked(self)
    }

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is NOT performed
    fn from_xy_unchecked(x: Self::Base, y: Self::Base) -> Self {
        <Self as CurveAffine>::from_xy_unchecked(x, y)
    }

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is performed
    fn from_xy_checked(x: Self::Base, y: Self::Base) -> Result<Self, GroupDecodingError> {
        <Self as CurveAffine>::from_xy_checked(x, y)
    }

    /// returns A coefficient for a short Weierstrass form
    fn a_coeff() -> Self::Base {
        <Self as CurveAffine>::a_coeff()
    }

    /// returns B coefficient for a short Weierstrass form
    fn b_coeff() -> Self::Base {
        <Self as CurveAffine>::b_coeff()
    }
}

impl<G: CurveProjective> GenericCurveProjective for G {
    type Scalar = <Self as CurveProjective>::Scalar;
    type Base = <Self as CurveProjective>::Base;
    type Affine = <Self as CurveProjective>::Affine;

    /// Returns the additive identity.
    fn zero() -> Self {
        <Self as CurveProjective>::zero()
    }

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self {
        <Self as CurveProjective>::one()
    }

    /// Determines if this point is the point at infinity.
    fn is_zero(&self) -> bool {
        <Self as CurveProjective>::is_zero(&self)
    }

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]){
        <Self as CurveProjective>::batch_normalization(v)
    }

    /// Checks if the point is already "normalized" so that
    /// cheap affine conversion is possible.
    fn is_normalized(&self) -> bool {
        <Self as CurveProjective>::is_normalized(&self)
    }

    /// Doubles this element.
    fn double(&mut self) {
        <Self as CurveProjective>::double(self)
    }

    /// Adds another element to this element.
    fn add_assign(&mut self, other: &Self) {
        <Self as CurveProjective>::add_assign(self, other)
    }

    /// Subtracts another element from this element.
    fn sub_assign(&mut self, other: &Self) {
        <Self as CurveProjective>::sub_assign(self, other)
    }

    /// Adds an affine element to this element.
    fn add_assign_mixed(&mut self, other: &Self::Affine) {
        <Self as CurveProjective>::add_assign_mixed(self, other);
    }

    /// Negates this element.
    fn negate(&mut self) {
        <Self as CurveProjective>::negate(self);
    }

    /// Performs scalar multiplication of this element.
    fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S) {
        <Self as CurveProjective>::mul_assign(self, other);
    }

    /// Converts this element into its affine representation.
    fn into_affine(&self) -> Self::Affine {
        <Self as CurveProjective>::into_affine(self)
    }

    /// Recommends a wNAF window table size given a scalar. Always returns a number
    /// between 2 and 22, inclusive.
    fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize {
        <Self as CurveProjective>::recommended_wnaf_for_scalar(scalar)
    }

    /// Recommends a wNAF window size given the number of scalars you intend to multiply
    /// a base by. Always returns a number between 2 and 22, inclusive.
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        <Self as CurveProjective>::recommended_wnaf_for_num_scalars(num_scalars)
    }

    /// Returns references to underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn as_xyz(&self) -> (&Self::Base, &Self::Base, &Self::Base) {
        <Self as CurveProjective>::as_xyz(self)
    }
    
    /// Returns underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn into_xyz_unchecked(self) -> (Self::Base, Self::Base, Self::Base) {
        <Self as CurveProjective>::into_xyz_unchecked(self)
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is NOT performed
    fn from_xyz_unchecked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Self {
        <Self as CurveProjective>::from_xyz_unchecked(_x, _y, _z)
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is performed
    fn from_xyz_checked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Result<Self, GroupDecodingError> {
        <Self as CurveProjective>::from_xyz_checked(_x, _y, _z)
    }
}