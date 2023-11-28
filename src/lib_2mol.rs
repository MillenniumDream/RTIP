//! RTIP
//!
//! RTIP (roto-translational invariant potential) is a biased potential appended to the real PES
//! (potential energy surface) to drive the molucule (aperiodic structure) escaping from local minimum
//! along the most flat directions.It aims at intelligent pathway searching for chemical and biological
//! reactions, like protein folding and enzyme catalysis.

//extern crate libc;
//extern crate mpi;

pub mod logo;
pub mod constants;
pub mod input;
pub mod output;
pub mod error;
pub mod matrix;
pub mod traits;
pub mod cp2k;
pub mod system;
pub mod potential;
pub mod rtip;
pub mod idwm;
pub mod optimization;
pub mod pathway_sampling;
pub mod md;
pub mod synthesis;

use mpi::traits::*;
use crate::cp2k::*;
use crate::system::System;
use crate::potential::*;
use crate::constants::*;
//use crate::rtip::*;
use crate::input::Para;
use crate::traits::*;
use crate::synthesis::*;
//use mpi::topology::UserCommunicator;
use mpi::topology::Color;
use ndarray::prelude::*;
//use ndarray_rand::RandomExt;
//use ndarray_rand::rand_distr::Uniform;
//use ndarray_linalg::{SVD, Determinant, EighInto, UPLO};
//use std::ffi::CString;





#[no_mangle]
pub extern fn main() -> i32
{
    println!("Begin Initialization");
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    cp2k_init_without_mpi();




    let s1 = System::read_xyz("1.xyz");
    let s2 = System::read_xyz("2.xyz");
//    let s3 = System::read_xyz("3.xyz");
//    let natom: usize = s1.natom + s2.natom + s3.natom;
    let natom: usize = s1.natom + s2.natom;
    let mut index1: Vec<usize> = vec![0; s1.natom];
    let mut index2: Vec<usize> = vec![0; s2.natom];
//    let mut index3: Vec<usize> = vec![0; s3.natom];

    let mut atom_type: Vec<Element> = vec![Element::H; natom];
    let mut coord: Array2<f64> = Array2::zeros((natom, 3));
    let mut index: usize = 0;

    for i in 0..s1.natom
    {
        index1[i] = index;
        atom_type[index] = s1.atom_type.as_ref().unwrap()[i].clone();
        coord[[index, 0]] = s1.coord[[i, 0]];
        coord[[index, 1]] = s1.coord[[i, 1]];
        coord[[index, 2]] = s1.coord[[i, 2]];
        index += 1;
    }
    for i in 0..s2.natom
    {
        index2[i] = index;
        atom_type[index] = s2.atom_type.as_ref().unwrap()[i].clone();
        coord[[index, 0]] = s2.coord[[i, 0]];
        coord[[index, 1]] = s2.coord[[i, 1]];
        coord[[index, 2]] = s2.coord[[i, 2]];
        index += 1;
    }
/*
    for i in 0..s3.natom
    {
        index3[i] = index;
        atom_type[index] = s3.atom_type.as_ref().unwrap()[i].clone();
        coord[[index, 0]] = s3.coord[[i, 0]];
        coord[[index, 1]] = s3.coord[[i, 1]];
        coord[[index, 2]] = s3.coord[[i, 2]];
        index += 1;
    }
*/

//    let mol_index: Vec<Vec<usize>> = vec![index1, index2, index3];
    let mol_index: Vec<Vec<usize>> = vec![index1, index2];
    let mut s = System
    {
        natom,
        coord,
        cell: None,
        atom_type: Some(atom_type),
        atom_add_pot: None,
        mutable: None,
        pot: 0.0,
    };





    {
        let my_comm = world.split_by_color(Color::with_value(world.rank() % 1)).unwrap();
        let (str_output_file, output_file) = output::output_rtip(&my_comm, Some(world.rank() % 1));
        let cp2k_output_file = output::output_cp2k(Some(world.rank() % 1));





        let para = Para::new();





        synthesis(&my_comm, &mut s, &mol_index, 5.0);
        if world.rank() == 0
        {
            s.write_xyz("IS.xyz", true, 0);
        }
        let synthesis_pot = SynthesisPot
        {
            initial_state: s,
            mol_index,
            para: &para,
            str_output_file,
            output_file,
        };
        let cp2k_pes = Cp2kPES::new("cp2k.inp", &cp2k_output_file, 4);
        synthesis_pot.rtip_path_sampling(&my_comm, &cp2k_pes);
    }





    cp2k_finalize_without_mpi();
    println!("Finalization done");
    return 0;
}










