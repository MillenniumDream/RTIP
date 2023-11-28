//! An interface to CP2K, which provides a variety of PES (i.e. total energy and atomic forces) for RTIP.

extern crate libc;

//use crate::constants::BOHR_TO_ANGSTROM;
use std::ffi::CString;
use crate::common::error::*;
use crate::pes_exploration::system::System;
use crate::pes_exploration::traits::PES;
use ndarray::Array2;

pub type ForceEnv = i32;





extern
{
//    fn __libcp2k_MOD_cp2k_get_version(version: *const libc::c_char, length: libc::c_int);

    fn __libcp2k_MOD_cp2k_init();
    fn __libcp2k_MOD_cp2k_init_without_mpi();
    fn __libcp2k_MOD_cp2k_finalize();
    fn __libcp2k_MOD_cp2k_finalize_without_mpi();

    fn __libcp2k_MOD_cp2k_create_force_env(new_force_env: &mut ForceEnv, input_file: *const libc::c_char, output_file: *const libc::c_char);
    fn __libcp2k_MOD_cp2k_create_force_env_comm(new_force_env: &mut ForceEnv, input_file: *const libc::c_char, output_file: *const libc::c_char, mpi_comm: libc::c_int);
    fn __libcp2k_MOD_cp2k_destroy_force_env(force_env: ForceEnv);

    fn __libcp2k_MOD_cp2k_get_natom(force_env: ForceEnv, natom: &mut libc::c_int);
    fn __libcp2k_MOD_cp2k_get_nparticle(force_env: ForceEnv, nparticle: &mut libc::c_int);
    fn __libcp2k_MOD_cp2k_get_positions(force_env: ForceEnv, pos: &mut libc::c_double, n_el: libc::c_int);
    fn __libcp2k_MOD_cp2k_get_cell(force_env: ForceEnv, cell: &mut libc::c_double);
    fn __libcp2k_MOD_cp2k_get_qmmm_cell(force_env: ForceEnv, cell: &mut libc::c_double);
    fn __libcp2k_MOD_cp2k_get_potential_energy(force_env: ForceEnv, e_pot: &mut libc::c_double);
    fn __libcp2k_MOD_cp2k_get_forces(force_env: ForceEnv, force: &mut libc::c_double, n_el: libc::c_int);
    fn __libcp2k_MOD_cp2k_get_result(force_env: ForceEnv, description: *const libc::c_char, result: &mut libc::c_double, n_el: libc::c_int);

    fn __libcp2k_MOD_cp2k_set_positions(force_env: ForceEnv, new_pos: *const libc::c_double, n_el: libc::c_int);
    fn __libcp2k_MOD_cp2k_set_cell(force_env: ForceEnv, new_cell: *const libc::c_double);
    fn __libcp2k_MOD_cp2k_set_velocities(force_env: ForceEnv, new_vel: *const libc::c_double, n_el: libc::c_int);

    fn __libcp2k_MOD_cp2k_calc_energy(force_env: ForceEnv);
    fn __libcp2k_MOD_cp2k_calc_energy_force(force_env: ForceEnv);
    fn __libcp2k_MOD_cp2k_run_input(input_file: *const libc::c_char, output_file: *const libc::c_char);
    fn __libcp2k_MOD_cp2k_run_input_comm(input_file: *const libc::c_char, output_file: *const libc::c_char, mpi_comm: libc::c_int);
}





/// Initialize CP2K and MPI.
///
/// Warning: This function could only be called once.
/// You are supposed to call cp2k_finalize() before exiting the program.
pub fn cp2k_init()
{
    unsafe
    {
        __libcp2k_MOD_cp2k_init();
    }
}

/// Initialize CP2K without initializing MPI.
///
/// Warning: This function could only be called after MPI initialization, and could only be called once.
/// You are supposed to call cp2k_finalize() before exiting the program.
pub fn cp2k_init_without_mpi()
{
    unsafe
    {
        __libcp2k_MOD_cp2k_init_without_mpi();
    }
}

/// Finalize CP2K and MPI.
pub fn cp2k_finalize()
{
    unsafe
    {
        __libcp2k_MOD_cp2k_finalize();
    }
}

/// Finalize CP2K without finalizing MPI.
///
/// Warning: You are sopposed to finalize MPI after this function being called.
pub fn cp2k_finalize_without_mpi()
{
    unsafe
    {
        __libcp2k_MOD_cp2k_finalize_without_mpi();
    }
}



/// Read from the input file, create a new force environment, and write the informations to the output file.
///
/// Warning: You are supposed to call cp2k_destroy_force_env() for cleanup before calling cp2k_finalize().
///
/// # Parameters
/// ```
/// new_force_env: the created new force environment
/// input_file: a file containing all the settings for CP2K
/// output_file: a file where CP2K is going to append its output to (created if non-existent)
/// ```
///
/// # Examples
/// ```
/// let mut force_env: ForceEnv = 0;
/// cp2k_create_force_env(&mut force_env, "input.cp2k", "cp2k.out");
/// ```
pub fn cp2k_create_force_env(new_force_env: &mut ForceEnv, input_file: &str, output_file: &str)
{
    let input_file = CString::new(input_file).expect(&error_str_to_cstring(input_file));
    let output_file = CString::new(output_file).expect(&error_str_to_cstring(output_file));

    unsafe
    {
        __libcp2k_MOD_cp2k_create_force_env(new_force_env, input_file.as_ptr(), output_file.as_ptr());
    }
}

/// A function like cp2k_create_force_env(). Create a new force environment in the given MPI communicator.
///
/// Warning: You are supposed to call cp2k_destroy_force_env() for cleanup before calling cp2k_finalize().
///
/// # Parameters
/// ```
/// new_force_env: the created new force environment
/// input_file: a file containing all the settings for CP2K
/// output_file: a file where CP2K is going to append its output to (created if non-existent)
/// mpi_comm: a given MPI communicator for calculation (if MPI is not managed by CP2K)
/// ```
///
/// # Examples
/// ```
/// let mut force_env: ForceEnv = 0;
/// let mpi_comm: i32 = 0;
/// cp2k_create_force_env_comm(force_env, "input.cp2k", "cp2k.out", mpi_comm);
/// ```
pub fn cp2k_create_force_env_comm(new_force_env: &mut ForceEnv, input_file: &str, output_file: &str, mpi_comm: i32)
{
    let input_file = CString::new(input_file).expect(&error_str_to_cstring(input_file));
    let output_file = CString::new(output_file).expect(&error_str_to_cstring(output_file));
    
    unsafe
    {
        __libcp2k_MOD_cp2k_create_force_env_comm(new_force_env, input_file.as_ptr(), output_file.as_ptr(), mpi_comm);
    }
}

/// Destroy a force environment.
///
/// # Parameters
/// ```
/// force_env: the force environment to be destroy
/// ```
pub fn cp2k_destroy_force_env(force_env: ForceEnv)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_destroy_force_env(force_env);
    }
}



/// Get the number of atoms from CP2K.
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// natom: a mutable reference for restoring the number of atoms
/// ```
///
/// # Examples
/// ```
/// let mut natom: i32 = 0;
/// cp2k_get_natom(force_env, &mut natom);
/// ```
pub fn cp2k_get_natom(force_env: ForceEnv, natom: &mut i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_natom(force_env, natom);
    }
}

/// Get the number of particles from CP2K.
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// nparticle: a mutable reference for restoring the number of particles
/// ```
///
/// # Examples
/// ```
/// let mut nparticle: i32 = 0;
/// cp2k_get_nparticle(force_env, &mut nparticle);
/// ```
pub fn cp2k_get_nparticle(force_env: ForceEnv, nparticle: &mut i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_nparticle(force_env, nparticle);
    }
}

/// Get the positions of the particles from CP2K (Unit: Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// pos: a reference to a mutable array (containing at least 3*nparticle elements) for restoring the positions of the particles
/// n_el: nparticle * 3
/// ```
///
/// # Examples
/// ```
/// let mut pos: Vec<f64> = vec![0.0; (nparticle*3).try_into().unwrap()];
/// let n_el: i32 = nparticle * 3;
/// cp2k_get_positions(force_env, &mut pos, n_el);
/// ```
pub fn cp2k_get_positions(force_env: ForceEnv, pos: &mut[f64], n_el: i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_positions(force_env, &mut pos[0], n_el);
    }
//    for i in 0..n_el.try_into().expect(&error_type_transformation("n_el", "i32", "usize"))
//    {
//        pos[i] *= BOHR_TO_ANGSTROM;
//    }
}

/// Get the cell from CP2K (Unit: Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// cell: a reference to a mutable array (containing 9 elements) for restoring the cell
/// ```
///
/// # Examples
/// ```
/// let mut cell: [f64; 9] = [0.0; 9];
/// cp2k_get_cell(force_env, &mut cell);
/// ```
pub fn cp2k_get_cell(force_env: ForceEnv, cell: &mut[f64])
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_cell(force_env, &mut cell[0]);
    }
//    for i in 0..9
//    {
//        cell[i] *= BOHR_TO_ANGSTROM;
//    }
}

/// Get the QM/MM cell from CP2K (Unit: Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// cell: a reference to a mutable array (containing 9 elements) for restoring the QM/MM cell
/// ```
///
/// # Examples
/// ```
/// let mut qmmm_cell: [f64; 9] = [0.0; 9];
/// cp2k_get_qmmm_cell(force_env, &mut qmmm_cell);
/// ```
pub fn cp2k_get_qmmm_cell(force_env: ForceEnv, cell: &mut[f64])
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_qmmm_cell(force_env, &mut cell[0]);
    }
//    for i in 0..9
//    {
//        cell[i] *= BOHR_TO_ANGSTROM;
//    }
}

/// Get the potential energy of the system from CP2K (Unit: Hartree).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// e_pot: a mutable reference for restoring the potential energy
/// ```
///
/// # Examples
/// ```
/// let mut e_pot: f64 = 0.0;
/// cp2k_get_potential_energy(force_env, &mut e_pot);
/// ```
pub fn cp2k_get_potential_energy(force_env: ForceEnv, e_pot: &mut f64)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_potential_energy(force_env, e_pot);
    }
}

/// Get the forces of the particles from CP2K (Unit: Hartree/Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// force: a reference to a mutable array (containing at least nparticle*3 elements) for restoring the forces of the particles
/// n_el: nparticle * 3
/// ```
///
/// # Examples
/// ```
/// let mut force: Vec<f64> = vec![0.0; (nparticle*3).try_into().unwrap()];
/// let n_el: i32 = nparticle * 3;
/// cp2k_get_forces(force_env, &mut force, n_el);
/// ```
pub fn cp2k_get_forces(force_env: ForceEnv, force: &mut[f64], n_el: i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_get_forces(force_env, &mut force[0], n_el);
    }
}

pub fn cp2k_get_result(force_env: ForceEnv, description: &str, result: &mut[f64], n_el: i32)
{
    let description = CString::new(description).expect(&error_str_to_cstring(description));
    unsafe
    {
        __libcp2k_MOD_cp2k_get_result(force_env, description.as_ptr(), &mut result[0], n_el);
    }
}



/// Set the positions of the particles for CP2K (Unit: Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// new_pos: an immutable reference to the array (containing at least nparticle*3 elements) containing the positions of the particles to be set
/// n_el: nparticle * 3
/// ```
///
/// # Examples
/// ```
/// let mut new_pos: Vec<f64> = vec![0.0; (nparticle*3).try_into().unwrap()];
/// let n_el: i32 = nparticle * 3;
/// cp2k_set_positions(force_env, &new_pos, n_el);
///
/// ```
pub fn cp2k_set_positions(force_env: ForceEnv, new_pos: &[f64], n_el: i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_set_positions(force_env, &new_pos[0], n_el);
    }
}

/// Set the cell of the system for CP2K (Unit: Bohr).
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// new_cell: an immutable reference to the array (containing at least 9 elements) containing the cell to be set
/// ```
///
/// # Examples
/// ```
/// let mut new_cell: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
/// cp2k_set_cell(force_env, &new_cell);
/// ```
pub fn cp2k_set_cell(force_env: ForceEnv, new_cell: &[f64])
{
    unsafe
    {
        __libcp2k_MOD_cp2k_set_cell(force_env, &new_cell[0]);
    }
}

/// Set the velocities of the particles for CP2K.
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// new_vel: an immutable reference to the array (containing at least nparticle*3 elements) containing the velocities of the particles to be set
/// n_el: nparticle * 3
/// ```
///
/// # Examples
/// ```
/// let mut new_vel: Vec<f64> = vec![0.0; (nparticle*3).try_into().unwrap()];
/// let n_el: i32 = nparticle * 3;
/// cp2k_set_velocities(force_env, &new_vel, n_el);
/// ```
pub fn cp2k_set_velocities(force_env: ForceEnv, new_vel: &[f64], n_el: i32)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_set_velocities(force_env, &new_vel[0], n_el);
    }
}



/// Calculate the energy of the system.
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// ```
pub fn cp2k_calc_energy(force_env: ForceEnv)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_calc_energy(force_env);
    }
}

/// Calculate the energy and atomic forces.
///
/// # Parameters
/// ```
/// force_env: the previously built force environment
/// ```
pub fn cp2k_calc_energy_force(force_env: ForceEnv)
{
    unsafe
    {
        __libcp2k_MOD_cp2k_calc_energy_force(force_env);
    }
}

/// Perform a CP2K run with the given input file.
///
/// Warning: Clean up the force environment before calling.
///
/// # Parameters
/// ```
/// input_file: a file containing all the settings for CP2K
/// output_file: seem not work (output to the default file mainLog.log)
/// ```
///
/// # Examples
/// ```
/// cp2k_run_input("input.cp2k", "cp2k.out");
/// ```
pub fn cp2k_run_input(input_file: &str, output_file: &str)
{
    let input_file = CString::new(input_file).expect(&error_str_to_cstring(input_file));
    let output_file = CString::new(output_file).expect(&error_str_to_cstring(output_file));

    unsafe
    {
        __libcp2k_MOD_cp2k_run_input(input_file.as_ptr(), output_file.as_ptr());
    }
}

/// Perform a CP2K run with the given input file in the given MPI communicator.
///
/// Warning: Clean up the force environment before calling.
///
/// # Parameters
/// ```
/// input_file: a file containing all the settings for CP2K.
/// output_file: seem not work (output to the default file mainLog.log)
/// mpi_comm: a given MPI communicator for calculation (if MPI is not managed by CP2K)
/// ```
///
/// # Examples
/// ```
/// cp2k_run_input_comm("input.cp2k", "cp2k.out", 0);
/// ```
pub fn cp2k_run_input_comm(input_file: &str, output_file: &str, mpi_comm: i32)
{
    let input_file = CString::new(input_file).expect(&error_str_to_cstring(input_file));
    let output_file = CString::new(output_file).expect(&error_str_to_cstring(output_file));

    unsafe
    {
        __libcp2k_MOD_cp2k_run_input_comm(input_file.as_ptr(), output_file.as_ptr(), mpi_comm);
    }
}





/// The structure containing the parameters for a CP2K potential energy surface (PES)
///
/// # Fields
/// ```
/// force_env: the force environment built in CP2K
/// ```
pub struct Cp2kPES
{
    force_env: ForceEnv,
}

impl Cp2kPES
{
    /// Input the input file name and output file name, create a CP2K potential energy surface (PES)
    ///
    /// # Parameters
    /// ```
    /// input_file: a file containing all the settings for CP2K
    /// output_file: a file where CP2K is going to append its output to (created if non-existent)
    /// ```
    pub fn new(input_file: &str, output_file: &str, mpi_comm: i32) -> Self
    {
        let mut force_env: ForceEnv = 0;
        cp2k_create_force_env_comm(&mut force_env, input_file, output_file, mpi_comm);
        Cp2kPES
        {
            force_env,
        }
    }
}

impl PES for Cp2kPES
{
    /// Input a structure, calculate and output its potential energy according to the CP2K potential energy surface (PES)
    ///
    /// # Parameters
    /// ```
    /// s: the input structure
    /// pot_cp2k: the output CP2K potential energy
    /// ```
    fn get_energy(&self, s: &System) -> f64
    {
        let mut pot_cp2k: f64 = 0.0;

        s.to_cp2k(self.force_env);
        cp2k_calc_energy(self.force_env);
        cp2k_get_potential_energy(self.force_env, &mut pot_cp2k);

        pot_cp2k
    }

    /// Input a structure,
    /// calculate and output its potential energy and atomic forces according to the CP2K potential energy surface (PES)
    ///
    /// # Parameters
    /// ```
    /// s: the input structure
    /// pot_cp2k: the output CP2K potential energy
    /// force_cp2k: the output CP2K atomic forces (s.coord.natom * 3)
    /// ```
    fn get_energy_force(&self, s: &System) -> (f64, Array2<f64>)
    {
        let n_el: i32 = (s.natom * 3).try_into().expect(&error_type_transformation("n_el", "usize", "i32"));
        let mut pot_cp2k: f64 = 0.0;
        let mut force_cp2k: Array2<f64> = Array2::zeros(s.coord.raw_dim());

        s.to_cp2k(self.force_env);
        cp2k_calc_energy_force(self.force_env);
        cp2k_get_potential_energy(self.force_env, &mut pot_cp2k);
        cp2k_get_forces(self.force_env, force_cp2k.as_slice_mut().expect(&error_as_slice("force_cp2k")), n_el);

        (pot_cp2k, force_cp2k)
    }
}

impl Drop for Cp2kPES
{
    /// Destroy the built force environment in the drop() function
    fn drop(&mut self)
    {
        cp2k_destroy_force_env(self.force_env);
    }
}










