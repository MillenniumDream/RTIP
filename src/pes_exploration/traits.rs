//! About the traits
use crate::pes_exploration::system::System;
use ndarray::Array2;
use mpi::traits::Communicator;





pub trait PES
{
    fn get_energy(&self, s: &System) -> f64;
    fn get_energy_force(&self, s: &System) -> (f64, Array2<f64>);
}



pub trait RtipPathSampling
{
    fn rtip_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P);
}



pub trait IdwmPathSampling
{
    fn idwm_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P);
}


pub trait RtipNVTMD
{
    fn rtip_nvt_md<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P);
}










