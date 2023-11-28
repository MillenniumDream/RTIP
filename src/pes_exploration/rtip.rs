//! This is the core module, containing the transformation and function for the roto-translational invariant potential (RTIP)

use crate::common::error::*;
use crate::pes_exploration::system::System;
use crate::pes_exploration::traits::PES;
use ndarray::{Array1, Array2, Axis, array};
use ndarray_linalg::{EighInto, UPLO};





/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the minimum roto-translational invariant (RTI) distance between the two atomic coordinates
/// under corresponding transformation (rot, tran), which makes ( coord2 - (coord1 * rot + tran) ) minimum.
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of System 1 (natom * 3)
/// coord2: the atomic coordinates of System 2 (natom * 3)
/// rti_dist: the output minimum roto-translational invariant (RTI) distance between coord1 and coord2
/// ```
pub fn rti_dist(coord1: &Array2<f64>, coord2: &Array2<f64>) -> f64
{
    // Move the geometric centers of the two systems to the origin
    let o1: Array1<f64> = coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
    let o2: Array1<f64> = coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
    let coord1: Array2<f64> = coord1 - &o1;
    let coord2: Array2<f64> = coord2 - &o2;

    // Find out the RTI distance using the quaternion method (i.e. the root of the minimum eigenvalue)
    let add: Array2<f64> = &coord1 + &coord2;
    let sub: Array2<f64> = &coord1 - &coord2;
    let mut s: Array2<f64> = Array2::zeros((4, 4));
    let mut a: Array2<f64>;
    for i in 0..coord1.len_of(Axis(0))
    {
        a = array!
        [
            [         0.0,  sub[[i,0]],  sub[[i,1]],  sub[[i,2]] ],
            [ -sub[[i,0]],         0.0, -add[[i,2]],  add[[i,1]] ],
            [ -sub[[i,1]],  add[[i,2]],         0.0, -add[[i,0]] ],
            [ -sub[[i,2]], -add[[i,1]],  add[[i,0]],         0.0 ],
        ];
        s += &a.t().dot(&a);
    }
    let (eigvals, _p) = s.eigh_into(UPLO::Upper).expect(&error_none_value("(eigvals, eigvecs)"));

    eigvals[0].sqrt()
}

/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the total four roto-translational invariant (RTI) distances between the two atomic coordinates.
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of System 1 (natom * 3)
/// coord2: the atomic coordinates of System 2 (natom * 3)
/// rti_dist: the output total four roto-translational invariant (RTI) distances between coord1 and coord2
/// ```
pub fn rti_dists(coord1: &Array2<f64>, coord2: &Array2<f64>) -> Vec<f64>
{
    // Move the geometric centers of the two systems to the origin
    let o1: Array1<f64> = coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
    let o2: Array1<f64> = coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
    let coord1: Array2<f64> = coord1 - &o1;
    let coord2: Array2<f64> = coord2 - &o2;

    // Find out the RTI distance using the quaternion method (i.e. the root of the minimum eigenvalue)
    let add: Array2<f64> = &coord1 + &coord2;
    let sub: Array2<f64> = &coord1 - &coord2;
    let mut s: Array2<f64> = Array2::zeros((4, 4));
    let mut a: Array2<f64>;
    for i in 0..coord1.len_of(Axis(0))
    {
        a = array!
        [
            [         0.0,  sub[[i,0]],  sub[[i,1]],  sub[[i,2]] ],
            [ -sub[[i,0]],         0.0, -add[[i,2]],  add[[i,1]] ],
            [ -sub[[i,1]],  add[[i,2]],         0.0, -add[[i,0]] ],
            [ -sub[[i,2]], -add[[i,1]],  add[[i,0]],         0.0 ],
        ];
        s += &a.t().dot(&a);
    }
    let (eigvals, _p) = s.eigh_into(UPLO::Upper).expect(&error_none_value("(eigvals, eigvecs)"));

    vec![eigvals[0].sqrt(), eigvals[1].sqrt(), eigvals[2].sqrt(), eigvals[3].sqrt()]
}





/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the minimum roto-translational invariant (RTI) distance between the two atomic coordinates
/// under corresponding transformation (rot, tran), which makes ( coord2 - (coord1 * rot + tran) ) minimum,
/// output the corresponding vector from coord1 to coord2
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of System 1 (natom * 3)
/// coord2: the atomic coordinates of System 2 (natom * 3)
/// rti_dist: the output minimum roto-translational invariant (RTI) distance between coord1 and coord2
/// rti_vec: the output roto-translational invariant (RTI) vector from coord1 to coord2
/// ```
pub fn rti_dist_vec(coord1: &Array2<f64>, coord2: &Array2<f64>) -> (f64, Array2<f64>)
{
    assert_eq!(coord1.len_of(Axis(0)), coord2.len_of(Axis(0)));
    assert_eq!(coord1.len_of(Axis(1)), 3);
    assert_eq!(coord2.len_of(Axis(1)), 3);

    // Move the geometric centers of the two systems to the origin
    let o1: Array1<f64> = coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
    let o2: Array1<f64> = coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
    let coord1: Array2<f64> = coord1 - &o1;
    let coord2: Array2<f64> = coord2 - &o2;

    // Find out the corresponding rotation for the RTI distance using the quaternion method (i.e. finding the minimum eigenvalue and the corresponding eigenvector)
    let add: Array2<f64> = &coord1 + &coord2;
    let sub: Array2<f64> = &coord1 - &coord2;
    let mut s: Array2<f64> = Array2::zeros((4, 4));
    let mut a: Array2<f64>;
    for i in 0..coord1.len_of(Axis(0))
    {
        a = array!
        [
            [         0.0,  sub[[i,0]],  sub[[i,1]],  sub[[i,2]] ],
            [ -sub[[i,0]],         0.0, -add[[i,2]],  add[[i,1]] ],
            [ -sub[[i,1]],  add[[i,2]],         0.0, -add[[i,0]] ],
            [ -sub[[i,2]], -add[[i,1]],  add[[i,0]],         0.0 ],
        ];
        s += &a.t().dot(&a);
    }
    let (eigvals, p) = s.eigh_into(UPLO::Upper).expect(&error_none_value("(eigvals, eigvecs)"));

    // Transform the quaternion representation (i.e. the eigenvector of the minimum eigenvalue) to the 3*3 rotational matrix
    let rot: Array2<f64> = array!
    [
        [ p[[0,0]]*p[[0,0]] + p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * (p[[1,0]]*p[[2,0]] + p[[0,0]]*p[[3,0]]), 2.0 * (p[[1,0]]*p[[3,0]] - p[[0,0]]*p[[2,0]]) ],
        [ 2.0 * (p[[1,0]]*p[[2,0]] - p[[0,0]]*p[[3,0]]), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] + p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * (p[[2,0]]*p[[3,0]] + p[[0,0]]*p[[1,0]]) ],
        [ 2.0 * (p[[1,0]]*p[[3,0]] + p[[0,0]]*p[[2,0]]), 2.0 * (p[[2,0]]*p[[3,0]] - p[[0,0]]*p[[1,0]]), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] + p[[3,0]]*p[[3,0]] ],
    ];
    let rti_vec: Array2<f64> = &coord2 - &coord1.dot(&rot);

    return (eigvals[0].sqrt(), rti_vec)
}

/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the total four roto-translational invariant (RTI) distances between the two atomic coordinates,
/// output the corresponding vectors from coord1 to coord2
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of System 1 (natom * 3)
/// coord2: the atomic coordinates of System 2 (natom * 3)
/// rti_dist: the output total four roto-translational invariant (RTI) distances between coord1 and coord2
/// rti_vec: the output total four roto-translational invariant (RTI) vectors from coord1 to coord2
/// ```
pub fn rti_dists_vecs(coord1: &Array2<f64>, coord2: &Array2<f64>) -> (Vec<f64>, Vec<Array2<f64>>)
{
    assert_eq!(coord1.len_of(Axis(0)), coord2.len_of(Axis(0)));
    assert_eq!(coord1.len_of(Axis(1)), 3);
    assert_eq!(coord2.len_of(Axis(1)), 3);

    // Move the geometric centers of the two systems to the origin
    let o1: Array1<f64> = coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
    let o2: Array1<f64> = coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
    let coord1: Array2<f64> = coord1 - &o1;
    let coord2: Array2<f64> = coord2 - &o2;

    // Find out the corresponding rotation for the RTI distance using the quaternion method (i.e. finding the minimum eigenvalue and the corresponding eigenvector)
    let add: Array2<f64> = &coord1 + &coord2;
    let sub: Array2<f64> = &coord1 - &coord2;
    let mut s: Array2<f64> = Array2::zeros((4, 4));
    let mut a: Array2<f64>;
    for i in 0..coord1.len_of(Axis(0))
    {
        a = array!
        [
            [         0.0,  sub[[i,0]],  sub[[i,1]],  sub[[i,2]] ],
            [ -sub[[i,0]],         0.0, -add[[i,2]],  add[[i,1]] ],
            [ -sub[[i,1]],  add[[i,2]],         0.0, -add[[i,0]] ],
            [ -sub[[i,2]], -add[[i,1]],  add[[i,0]],         0.0 ],
        ];
        s += &a.t().dot(&a);
    }
    let (eigvals, p) = s.eigh_into(UPLO::Upper).expect(&error_none_value("(eigvals, eigvecs)"));

    // Transform the quaternion representation (i.e. the eigenvector of the minimum eigenvalue) to the 3*3 rotational matrix
    let rot0: Array2<f64> = array!
    [
        [ p[[0,0]]*p[[0,0]] + p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * (p[[1,0]]*p[[2,0]] + p[[0,0]]*p[[3,0]]), 2.0 * (p[[1,0]]*p[[3,0]] - p[[0,0]]*p[[2,0]]) ],
        [ 2.0 * (p[[1,0]]*p[[2,0]] - p[[0,0]]*p[[3,0]]), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] + p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * (p[[2,0]]*p[[3,0]] + p[[0,0]]*p[[1,0]]) ],
        [ 2.0 * (p[[1,0]]*p[[3,0]] + p[[0,0]]*p[[2,0]]), 2.0 * (p[[2,0]]*p[[3,0]] - p[[0,0]]*p[[1,0]]), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] + p[[3,0]]*p[[3,0]] ],
    ];
    let rti_vec0: Array2<f64> = &coord2 - &coord1.dot(&rot0);
    let rot1: Array2<f64> = array!
    [
        [ p[[0,1]]*p[[0,1]] + p[[1,1]]*p[[1,1]] - p[[2,1]]*p[[2,1]] - p[[3,1]]*p[[3,1]], 2.0 * (p[[1,1]]*p[[2,1]] + p[[0,1]]*p[[3,1]]), 2.0 * (p[[1,1]]*p[[3,1]] - p[[0,1]]*p[[2,1]]) ],
        [ 2.0 * (p[[1,1]]*p[[2,1]] - p[[0,1]]*p[[3,1]]), p[[0,1]]*p[[0,1]] - p[[1,1]]*p[[1,1]] + p[[2,1]]*p[[2,1]] - p[[3,1]]*p[[3,1]], 2.0 * (p[[2,1]]*p[[3,1]] + p[[0,1]]*p[[1,1]]) ],
        [ 2.0 * (p[[1,1]]*p[[3,1]] + p[[0,1]]*p[[2,1]]), 2.0 * (p[[2,1]]*p[[3,1]] - p[[0,1]]*p[[1,1]]), p[[0,1]]*p[[0,1]] - p[[1,1]]*p[[1,1]] - p[[2,1]]*p[[2,1]] + p[[3,1]]*p[[3,1]] ],
    ];
    let rti_vec1: Array2<f64> = &coord2 - &coord1.dot(&rot1);
    let rot2: Array2<f64> = array!
    [
        [ p[[0,2]]*p[[0,2]] + p[[1,2]]*p[[1,2]] - p[[2,2]]*p[[2,2]] - p[[3,2]]*p[[3,2]], 2.0 * (p[[1,2]]*p[[2,2]] + p[[0,2]]*p[[3,2]]), 2.0 * (p[[1,2]]*p[[3,2]] - p[[0,2]]*p[[2,2]]) ],
        [ 2.0 * (p[[1,2]]*p[[2,2]] - p[[0,2]]*p[[3,2]]), p[[0,2]]*p[[0,2]] - p[[1,2]]*p[[1,2]] + p[[2,2]]*p[[2,2]] - p[[3,2]]*p[[3,2]], 2.0 * (p[[2,2]]*p[[3,2]] + p[[0,2]]*p[[1,2]]) ],
        [ 2.0 * (p[[1,2]]*p[[3,2]] + p[[0,2]]*p[[2,2]]), 2.0 * (p[[2,2]]*p[[3,2]] - p[[0,2]]*p[[1,2]]), p[[0,2]]*p[[0,2]] - p[[1,2]]*p[[1,2]] - p[[2,2]]*p[[2,2]] + p[[3,2]]*p[[3,2]] ],
    ];
    let rti_vec2: Array2<f64> = &coord2 - &coord1.dot(&rot2);
    let rot3: Array2<f64> = array!
    [
        [ p[[0,3]]*p[[0,3]] + p[[1,3]]*p[[1,3]] - p[[2,3]]*p[[2,3]] - p[[3,3]]*p[[3,3]], 2.0 * (p[[1,3]]*p[[2,3]] + p[[0,3]]*p[[3,3]]), 2.0 * (p[[1,3]]*p[[3,3]] - p[[0,3]]*p[[2,3]]) ],
        [ 2.0 * (p[[1,3]]*p[[2,3]] - p[[0,3]]*p[[3,3]]), p[[0,3]]*p[[0,3]] - p[[1,3]]*p[[1,3]] + p[[2,3]]*p[[2,3]] - p[[3,3]]*p[[3,3]], 2.0 * (p[[2,3]]*p[[3,3]] + p[[0,3]]*p[[1,3]]) ],
        [ 2.0 * (p[[1,3]]*p[[3,3]] + p[[0,3]]*p[[2,3]]), 2.0 * (p[[2,3]]*p[[3,3]] - p[[0,3]]*p[[1,3]]), p[[0,3]]*p[[0,3]] - p[[1,3]]*p[[1,3]] - p[[2,3]]*p[[2,3]] + p[[3,3]]*p[[3,3]] ],
    ];
    let rti_vec3: Array2<f64> = &coord2 - &coord1.dot(&rot3);

    return (vec![eigvals[0].sqrt(), eigvals[1].sqrt(), eigvals[2].sqrt(), eigvals[3].sqrt()], vec![rti_vec0, rti_vec1, rti_vec2, rti_vec3])
}





/// Input the atomic coordinates of two systems (coord1, coord2),
/// find out the corresponding transformation (rot, tran) for the roto-translational invariant (RTI) distance between the two atomic coordinates,
/// which makes ( coord2 - (coord1 * rot + tran) ) minimum
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of System 1 (natom * 3)
/// coord2: the atomic coordinates of System 2 (natom * 3)
/// rot: the output rotational matrix (3 * 3) for rti_dist
/// tran: the output translational vector (1 * 3) for rti_dist
/// ```
pub fn rti_rot_tran(coord1: &Array2<f64>, coord2: &Array2<f64>) -> (Array2<f64>, Array1<f64>)
{
    assert_eq!(coord1.len_of(Axis(0)), coord2.len_of(Axis(0)));
    assert_eq!(coord1.len_of(Axis(1)), 3);
    assert_eq!(coord2.len_of(Axis(1)), 3);

    // Move the geometric centers of the two systems to the origin
    let o1: Array1<f64> = coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
    let o2: Array1<f64> = coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
    let coord1: Array2<f64> = coord1 - &o1;
    let coord2: Array2<f64> = coord2 - &o2;

    // Find out the corresponding rotation for the RTI distance using the quaternion method (i.e. finding the minimun eigenvalue and the corresponding eigenvector)
    let add: Array2<f64> = &coord1 + &coord2;
    let sub: Array2<f64> = &coord1 - &coord2;
    let mut s: Array2<f64> = Array2::zeros((4, 4));
    let mut a: Array2<f64>;
    for i in 0..coord1.len_of(Axis(0))
    {
        a = array!
        [
            [         0.0,  sub[[i,0]],  sub[[i,1]],  sub[[i,2]] ],
            [ -sub[[i,0]],         0.0, -add[[i,2]],  add[[i,1]] ],
            [ -sub[[i,1]],  add[[i,2]],         0.0, -add[[i,0]] ],
            [ -sub[[i,2]], -add[[i,1]],  add[[i,0]],         0.0 ],
        ];
        s += &a.t().dot(&a);
    }
    let (_eigvals, p) = s.eigh_into(UPLO::Upper).expect(&error_none_value("(eigvals, eigvecs)"));

    // Transform the quaternion representation (the eigenvector of the minimum eigenvalue) to the 3*3 rotational matrix, and obtain the translational vector
    let rot: Array2<f64> = array!
    [
        [ p[[0,0]]*p[[0,0]] + p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * ( p[[1,0]]*p[[2,0]] + p[[0,0]]*p[[3,0]] ), 2.0 * ( p[[1,0]]*p[[3,0]] - p[[0,0]]*p[[2,0]] ) ],
        [ 2.0 * ( p[[1,0]]*p[[2,0]] - p[[0,0]]*p[[3,0]] ), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] + p[[2,0]]*p[[2,0]] - p[[3,0]]*p[[3,0]], 2.0 * ( p[[2,0]]*p[[3,0]] + p[[0,0]]*p[[1,0]] ) ],
        [ 2.0 * ( p[[1,0]]*p[[3,0]] + p[[0,0]]*p[[2,0]] ), 2.0 * ( p[[2,0]]*p[[3,0]] - p[[0,0]]*p[[1,0]] ), p[[0,0]]*p[[0,0]] - p[[1,0]]*p[[1,0]] - p[[2,0]]*p[[2,0]] + p[[3,0]]*p[[3,0]] ],
    ];
    let tran: Array1<f64> = o2 - o1.dot(&rot);

    (rot, tran)
}





fn f(x: f64) -> f64
{
    1.0 / x.powi(7)
}

fn df(x: f64) -> f64
{
    -7.0 / x.powi(8)
}

/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the roto-translational invariant (RTI) potential energy between coord1 and coord2
/// e_rtip = a * exp( -r^2 / (2*sigma^2) ), a Gaussian function where r is the RTI distance between coord1 and coord2
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of system 1 (natom * 3)
/// coord2: the atomic coordinates of system 2 (natom * 3)
/// a: the height of the Gaussian function
/// sigma: the width of the Gaussian function
/// pot_rti: the output RTI potential energy
/// ```
fn rti_pot(coord1: &Array2<f64>, coord2: &Array2<f64>, a: f64, sigma: f64) -> f64
{
    let rti_dists: Vec<f64> = rti_dists(coord1, coord2);

    let weight0: f64 = f(rti_dists[0]);
    let weight1: f64 = f(rti_dists[1]);
    let weight2: f64 = f(rti_dists[2]);
    let weight3: f64 = f(rti_dists[3]);
    let sum: f64 = weight0 + weight1 + weight2 + weight3;

    let pot_rti0: f64 = a * (weight0/sum) * ( -(rti_dists[0]*rti_dists[0]) / (2.0*sigma*sigma) ).exp();
    let pot_rti1: f64 = a * (weight1/sum) * ( -(rti_dists[1]*rti_dists[1]) / (2.0*sigma*sigma) ).exp();
    let pot_rti2: f64 = a * (weight2/sum) * ( -(rti_dists[2]*rti_dists[2]) / (2.0*sigma*sigma) ).exp();
    let pot_rti3: f64 = a * (weight3/sum) * ( -(rti_dists[3]*rti_dists[3]) / (2.0*sigma*sigma) ).exp();

    pot_rti0 + pot_rti1 + pot_rti2 + pot_rti3
}

/// Input the atomic coordinates of two systems (coord1, coord2),
/// output the roto-translational invariant (RTI) potential energy between coord1 and coord2,
/// and the corresponding atomic forces on coord2.
/// e_rtip = a * exp( -r^2 / (2*sigma^2) ), a Gaussian function where r is the RTI distance between coord1 and coord2
///
/// # Parameters
/// ```
/// coord1: the atomic coordinates of system 1 (natom * 3)
/// coord2: the atomic coordinates of system 2 (natom * 3)
/// a: the height of the Gaussian function
/// sigma: the width of the Gaussian function
/// pot_rti: the output RTI potential energy
/// force_rti: the output RTI atomic forces on coord2
/// ```
fn rti_pot_force(coord1: &Array2<f64>, coord2: &Array2<f64>, a: f64, sigma: f64) -> (f64, Array2<f64>)
{
    let (rti_dists, rti_vecs) = rti_dists_vecs(coord1, coord2);

    let weight0: f64 = f(rti_dists[0]);
    let weight1: f64 = f(rti_dists[1]);
    let weight2: f64 = f(rti_dists[2]);
    let weight3: f64 = f(rti_dists[3]);
    let sum: f64 = weight0 + weight1 + weight2 + weight3;

    let dweight0: f64 = df(rti_dists[0]);
    let dweight1: f64 = df(rti_dists[1]);
    let dweight2: f64 = df(rti_dists[2]);
    let dweight3: f64 = df(rti_dists[3]);

    let u0: f64 = a * ( -(rti_dists[0]*rti_dists[0]) / (2.0*sigma*sigma) ).exp();
    let u1: f64 = a * ( -(rti_dists[1]*rti_dists[1]) / (2.0*sigma*sigma) ).exp();
    let u2: f64 = a * ( -(rti_dists[2]*rti_dists[2]) / (2.0*sigma*sigma) ).exp();
    let u3: f64 = a * ( -(rti_dists[3]*rti_dists[3]) / (2.0*sigma*sigma) ).exp();

    let pot_rti0: f64 = (weight0/sum) * u0;
    let pot_rti1: f64 = (weight1/sum) * u1;
    let pot_rti2: f64 = (weight2/sum) * u2;
    let pot_rti3: f64 = (weight3/sum) * u3;

    let force_rti0: Array2<f64> = &rti_vecs[0] * ( pot_rti0/(sigma*sigma) + dweight0/(sum*sum*rti_dists[0]) * ( weight1*(u1-u0) + weight2*(u2-u0) + weight3*(u3-u0) ) );
    let force_rti1: Array2<f64> = &rti_vecs[1] * ( pot_rti1/(sigma*sigma) + dweight1/(sum*sum*rti_dists[1]) * ( weight0*(u0-u1) + weight2*(u2-u1) + weight3*(u3-u1) ) );
    let force_rti2: Array2<f64> = &rti_vecs[2] * ( pot_rti2/(sigma*sigma) + dweight2/(sum*sum*rti_dists[2]) * ( weight0*(u0-u2) + weight1*(u1-u2) + weight3*(u3-u2) ) );
    let force_rti3: Array2<f64> = &rti_vecs[3] * ( pot_rti3/(sigma*sigma) + dweight3/(sum*sum*rti_dists[3]) * ( weight0*(u0-u3) + weight1*(u1-u3) + weight2*(u2-u3) ) );

    return (pot_rti0 + pot_rti1 + pot_rti2 + pot_rti3, force_rti0 + force_rti1 + force_rti2 + force_rti3)
}





/// The structure containing the parameters for a roto-translational invariant (RTI) repulsive potential energy surface (PES)
///
/// # Fields
/// ```
/// local_min: a local minimum structure, where the RTI repulsive potential (a Gaussian function) centered
/// a: the height of the Gaussian function
/// sigma: the width of the Gaussian function
/// ```
pub struct Rtip0PES<'a>
{
    pub local_min: &'a System,
    pub nearby_ts: &'a Vec<System>,
    pub a_min: f64,
    pub a_ts: f64,
    pub sigma_min: f64,
    pub sigma_ts: Vec<f64>,
}

impl PES for Rtip0PES<'_>
{
    /// Input a structure,
    /// calculate and output its potential energy according to the roto-translational invariant (RTI) repulsive potential energy surface (PES)
    fn get_energy(&self, s: &System) -> f64
    {
        match (&self.local_min.atom_add_pot, &s.atom_add_pot)
        {
            // If the local minimum and input structure both do not specify the atoms to be added RTI potential
            (None, None) =>
            {
                let mut pot_rtip: f64 = rti_pot(&self.local_min.coord, &s.coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_rtip += rti_pot(&self.nearby_ts[i].coord, &s.coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_rtip
            },

            // If the local minimum does not specify the atoms to be added RTI potential but the input structure specifies
            (None, Some(atom_add_pot_s)) =>
            {
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot_s.len(), 3));
                for i in 0..atom_add_pot_s.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot_s[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot_s[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot_s[i], 2]];
                }

                let mut pot_rtip: f64 = rti_pot(&self.local_min.coord, &s_fragment_coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_rtip += rti_pot(&self.nearby_ts[i].coord, &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_rtip
            },

            // If the local minimum specifies the atoms to be added RTI potential but the input structure does not specify
            (Some(atom_add_pot_min), None) =>
            {
                let mut local_min_coord: Array2<f64> = Array2::zeros((atom_add_pot_min.len(), 3));
                let mut nearby_ts_coord: Vec< Array2<f64> > = vec![local_min_coord.clone(); self.nearby_ts.len()];
                for i in 0..atom_add_pot_min.len()
                {
                    local_min_coord[[i, 0]] = self.local_min.coord[[atom_add_pot_min[i], 0]];
                    local_min_coord[[i, 1]] = self.local_min.coord[[atom_add_pot_min[i], 1]];
                    local_min_coord[[i, 2]] = self.local_min.coord[[atom_add_pot_min[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot_min.len()
                    {
                        nearby_ts_coord[j][[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 0]];
                        nearby_ts_coord[j][[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 1]];
                        nearby_ts_coord[j][[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 2]];
                    }
                }

                let mut pot_rtip: f64 = rti_pot(&local_min_coord, &s.coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_rtip += rti_pot(&nearby_ts_coord[i], &s.coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_rtip
            },

            // If the local minimum and input structure both specify the atoms to be added RTI potential
            (Some(atom_add_pot_min), Some(atom_add_pot_s)) =>
            {
                let mut local_min_coord: Array2<f64> = Array2::zeros((atom_add_pot_min.len(), 3));
                let mut nearby_ts_coord: Vec< Array2<f64> > = vec![local_min_coord.clone(); self.nearby_ts.len()];
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot_s.len(), 3));
                for i in 0..atom_add_pot_min.len()
                {
                    local_min_coord[[i, 0]] = self.local_min.coord[[atom_add_pot_min[i], 0]];
                    local_min_coord[[i, 1]] = self.local_min.coord[[atom_add_pot_min[i], 1]];
                    local_min_coord[[i, 2]] = self.local_min.coord[[atom_add_pot_min[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot_min.len()
                    {
                        nearby_ts_coord[j][[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 0]];
                        nearby_ts_coord[j][[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 1]];
                        nearby_ts_coord[j][[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 2]];
                    }
                }
                for i in 0..atom_add_pot_s.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot_s[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot_s[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot_s[i], 2]];
                }

                let mut pot_rtip: f64 = rti_pot(&local_min_coord, &s_fragment_coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_rtip += rti_pot(&nearby_ts_coord[i], &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_rtip
            },
        }
    }

    /// Input a structure,
    /// calculate and output its potential energy and atomic forces according to the roto-translational invariant (RTI) repulsive potential energy surface (PES)
    fn get_energy_force(&self, s: &System) -> (f64, Array2<f64>)
    {
        match (&self.local_min.atom_add_pot, &s.atom_add_pot)
        {
            // If the local minimum and input structure both do not specify the atoms to be added RTI potential
            (None, None) =>
            {
                let (mut pot_rtip, mut force_rtip) = rti_pot_force(&self.local_min.coord, &s.coord, self.a_min, self.sigma_min);
                let mut pot_rtip_ts: f64;
                let mut force_rtip_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_rtip_ts, force_rtip_ts) = rti_pot_force(&self.nearby_ts[i].coord, &s.coord, self.a_ts, self.sigma_ts[i]);
                    pot_rtip += pot_rtip_ts;
                    force_rtip += &force_rtip_ts;
                }

                return (pot_rtip, force_rtip)
            },

            // If the local minimum does not specify the atoms to be added RTI potential but the input structure specifies
            (None, Some(atom_add_pot_s)) =>
            {
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot_s.len(), 3));
                for i in 0..atom_add_pot_s.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot_s[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot_s[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot_s[i], 2]];
                }

                let (mut pot_rtip, mut force_rtip_fragment) = rti_pot_force(&self.local_min.coord, &s_fragment_coord, self.a_min, self.sigma_min);
                let mut pot_rtip_ts: f64;
                let mut force_rtip_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_rtip_ts, force_rtip_ts) = rti_pot_force(&self.nearby_ts[i].coord, &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                    pot_rtip += pot_rtip_ts;
                    force_rtip_fragment += &force_rtip_ts;
                }

                let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
                for i in 0..atom_add_pot_s.len()
                {
                    force_rtip[[atom_add_pot_s[i], 0]] = force_rtip_fragment[[i, 0]];
                    force_rtip[[atom_add_pot_s[i], 1]] = force_rtip_fragment[[i, 1]];
                    force_rtip[[atom_add_pot_s[i], 2]] = force_rtip_fragment[[i, 2]];
                }

                return (pot_rtip, force_rtip)
            },

            // If the local minimum specifies the atoms to be added RTI potential but the input structure does not specify
            (Some(atom_add_pot_min), None) =>
            {
                let mut local_min_coord: Array2<f64> = Array2::zeros((atom_add_pot_min.len(), 3));
                let mut nearby_ts_coord: Vec< Array2<f64> > = vec![local_min_coord.clone(); self.nearby_ts.len()];
                for i in 0..atom_add_pot_min.len()
                {
                    local_min_coord[[i, 0]] = self.local_min.coord[[atom_add_pot_min[i], 0]];
                    local_min_coord[[i, 1]] = self.local_min.coord[[atom_add_pot_min[i], 1]];
                    local_min_coord[[i, 2]] = self.local_min.coord[[atom_add_pot_min[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot_min.len()
                    {
                        nearby_ts_coord[j][[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 0]];
                        nearby_ts_coord[j][[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 1]];
                        nearby_ts_coord[j][[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 2]];
                    }
                }

                let (mut pot_rtip, mut force_rtip) = rti_pot_force(&local_min_coord, &s.coord, self.a_min, self.sigma_min);
                let mut pot_rtip_ts: f64;
                let mut force_rtip_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_rtip_ts, force_rtip_ts) = rti_pot_force(&nearby_ts_coord[i], &s.coord, self.a_ts, self.sigma_ts[i]);
                    pot_rtip += pot_rtip_ts;
                    force_rtip += &force_rtip_ts;
                }

                return (pot_rtip, force_rtip)
            },

            // If the local minimum and input structure both specify the atoms to be added RTI potential
            (Some(atom_add_pot_min), Some(atom_add_pot_s)) =>
            {
                let mut local_min_coord: Array2<f64> = Array2::zeros((atom_add_pot_min.len(), 3));
                let mut nearby_ts_coord: Vec< Array2<f64> > = vec![local_min_coord.clone(); self.nearby_ts.len()];
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot_s.len(), 3));
                for i in 0..atom_add_pot_min.len()
                {
                    local_min_coord[[i, 0]] = self.local_min.coord[[atom_add_pot_min[i], 0]];
                    local_min_coord[[i, 1]] = self.local_min.coord[[atom_add_pot_min[i], 1]];
                    local_min_coord[[i, 2]] = self.local_min.coord[[atom_add_pot_min[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot_min.len()
                    {
                        nearby_ts_coord[j][[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 0]];
                        nearby_ts_coord[j][[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 1]];
                        nearby_ts_coord[j][[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot_min[i], 2]];
                    }
                }
                for i in 0..atom_add_pot_s.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot_s[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot_s[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot_s[i], 2]];
                }

                let (mut pot_rtip, mut force_rtip_fragment) = rti_pot_force(&local_min_coord, &s_fragment_coord, self.a_min, self.sigma_min);
                let mut pot_rtip_ts: f64;
                let mut force_rtip_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_rtip_ts, force_rtip_ts) = rti_pot_force(&nearby_ts_coord[i], &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                    pot_rtip += pot_rtip_ts;
                    force_rtip_fragment += &force_rtip_ts;
                }

                let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
                for i in 0..atom_add_pot_s.len()
                {
                    force_rtip[[atom_add_pot_s[i], 0]] = force_rtip_fragment[[i, 0]];
                    force_rtip[[atom_add_pot_s[i], 1]] = force_rtip_fragment[[i, 1]];
                    force_rtip[[atom_add_pot_s[i], 2]] = force_rtip_fragment[[i, 2]];
                }

                return (pot_rtip, force_rtip)
            },
        }
    }
}










