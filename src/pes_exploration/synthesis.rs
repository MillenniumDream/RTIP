use crate::common::constants::*;
use crate::common::error::*;
use crate::pes_exploration::system::System;
use crate::matrix;
use mpi::traits::*;
use ndarray::{Array2, Axis};





pub fn synthesis<C: Communicator>(comm: &C, s: &mut System, mol_index: &Vec<Vec<usize>>, dist: f64)
{
    // Parallel processing
    let root_process = comm.process_at_rank(ROOT_RANK);

    match mol_index.len()
    {
        2 =>
        {
            // Copy the atomic coordinates of the two molecules
            let mut coord1: Array2<f64> = Array2::zeros((mol_index[0].len(), 3));
            let mut coord2: Array2<f64> = Array2::zeros((mol_index[1].len(), 3));
            for i in 0..mol_index[0].len()
            {
                coord1[[i, 0]] = s.coord[[mol_index[0][i], 0]];
                coord1[[i, 1]] = s.coord[[mol_index[0][i], 1]];
                coord1[[i, 2]] = s.coord[[mol_index[0][i], 2]];
            }
            for i in 0..mol_index[1].len()
            {
                coord2[[i, 0]] = s.coord[[mol_index[1][i], 0]];
                coord2[[i, 1]] = s.coord[[mol_index[1][i], 1]];
                coord2[[i, 2]] = s.coord[[mol_index[1][i], 2]];
            }

            // Move the geometric centers of the two molecules to the origin
            coord1 -= &coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
            coord2 -= &coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));

            // Construct the rotational matrix
            let mut rot1: Array2<f64> = matrix::rand_rot();
            let mut rot2: Array2<f64> = matrix::rand_rot();
            let rot1_p: &mut [f64] = rot1.as_slice_mut().expect(&error_as_slice("rot1"));
            let rot2_p: &mut [f64] = rot2.as_slice_mut().expect(&error_as_slice("rot2"));
            root_process.broadcast_into(rot1_p);
            root_process.broadcast_into(rot2_p);

            // Rotate the two molecules
            coord1 = coord1.dot(&rot1);
            coord2 = coord2.dot(&rot2);

            // Construct the synthesis potential energy surface (PES)
            for i in 0..mol_index[0].len()
            {
                s.coord[[mol_index[0][i], 0]] = coord1[[i, 0]] - dist;
                s.coord[[mol_index[0][i], 1]] = coord1[[i, 1]];
                s.coord[[mol_index[0][i], 2]] = coord1[[i, 2]];
            }
            for i in 0..mol_index[1].len()
            {
                s.coord[[mol_index[1][i], 0]] = coord2[[i, 0]] + dist;
                s.coord[[mol_index[1][i], 1]] = coord2[[i, 1]];
                s.coord[[mol_index[1][i], 2]] = coord2[[i, 2]];
            }
        },

        3 =>
        {
            // Copy the atomic coordinates of the three molecules
            let mut coord1: Array2<f64> = Array2::zeros((mol_index[0].len(), 3));
            let mut coord2: Array2<f64> = Array2::zeros((mol_index[1].len(), 3));
            let mut coord3: Array2<f64> = Array2::zeros((mol_index[2].len(), 3));
            for i in 0..mol_index[0].len()
            {
                coord1[[i, 0]] = s.coord[[mol_index[0][i], 0]];
                coord1[[i, 1]] = s.coord[[mol_index[0][i], 1]];
                coord1[[i, 2]] = s.coord[[mol_index[0][i], 2]];
            }
            for i in 0..mol_index[1].len()
            {
                coord2[[i, 0]] = s.coord[[mol_index[1][i], 0]];
                coord2[[i, 1]] = s.coord[[mol_index[1][i], 1]];
                coord2[[i, 2]] = s.coord[[mol_index[1][i], 2]];
            }
            for i in 0..mol_index[2].len()
            {
                coord3[[i, 0]] = s.coord[[mol_index[2][i], 0]];
                coord3[[i, 1]] = s.coord[[mol_index[2][i], 1]];
                coord3[[i, 2]] = s.coord[[mol_index[2][i], 2]];
            }

            // Move the geometric centers of the three molecules to the origin
            coord1 -= &coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
            coord2 -= &coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
            coord3 -= &coord3.mean_axis(Axis(0)).expect(&error_none_value("coord3"));

            // Construct the rotational matrix
            let mut rot1: Array2<f64> = matrix::rand_rot();
            let mut rot2: Array2<f64> = matrix::rand_rot();
            let mut rot3: Array2<f64> = matrix::rand_rot();
            let rot1_p: &mut [f64] = rot1.as_slice_mut().expect(&error_as_slice("rot1"));
            let rot2_p: &mut [f64] = rot2.as_slice_mut().expect(&error_as_slice("rot2"));
            let rot3_p: &mut [f64] = rot3.as_slice_mut().expect(&error_as_slice("rot3"));
            root_process.broadcast_into(rot1_p);
            root_process.broadcast_into(rot2_p);
            root_process.broadcast_into(rot3_p);

            // Rotate the three molecules
            coord1 = coord1.dot(&rot1);
            coord2 = coord2.dot(&rot2);
            coord3 = coord3.dot(&rot3);

            // Construct the synthesis potential energy surface (PES)
            for i in 0..mol_index[0].len()
            {
                s.coord[[mol_index[0][i], 0]] = coord1[[i, 0]];
                s.coord[[mol_index[0][i], 1]] = coord1[[i, 1]] + dist;
                s.coord[[mol_index[0][i], 2]] = coord1[[i, 2]];
            }
            for i in 0..mol_index[1].len()
            {
                s.coord[[mol_index[1][i], 0]] = coord2[[i, 0]] + dist * 0.866025403784439_f64;
                s.coord[[mol_index[1][i], 1]] = coord2[[i, 1]] - dist * 0.5_f64;
                s.coord[[mol_index[1][i], 2]] = coord2[[i, 2]];
            }
            for i in 0..mol_index[2].len()
            {
                s.coord[[mol_index[2][i], 0]] = coord3[[i, 0]] - dist * 0.866025403784439_f64;
                s.coord[[mol_index[2][i], 1]] = coord3[[i, 1]] - dist * 0.5_f64;
                s.coord[[mol_index[2][i], 2]] = coord3[[i, 2]];
            }
        },

        4=>
        {
            // Copy the atomic coordinates of the four molecules
            let mut coord1: Array2<f64> = Array2::zeros((mol_index[0].len(), 3));
            let mut coord2: Array2<f64> = Array2::zeros((mol_index[1].len(), 3));
            let mut coord3: Array2<f64> = Array2::zeros((mol_index[2].len(), 3));
            let mut coord4: Array2<f64> = Array2::zeros((mol_index[3].len(), 3));
            for i in 0..mol_index[0].len()
            {
                coord1[[i, 0]] = s.coord[[mol_index[0][i], 0]];
                coord1[[i, 1]] = s.coord[[mol_index[0][i], 1]];
                coord1[[i, 2]] = s.coord[[mol_index[0][i], 2]];
            }
            for i in 0..mol_index[1].len()
            {
                coord2[[i, 0]] = s.coord[[mol_index[1][i], 0]];
                coord2[[i, 1]] = s.coord[[mol_index[1][i], 1]];
                coord2[[i, 2]] = s.coord[[mol_index[1][i], 2]];
            }
            for i in 0..mol_index[2].len()
            {
                coord3[[i, 0]] = s.coord[[mol_index[2][i], 0]];
                coord3[[i, 1]] = s.coord[[mol_index[2][i], 1]];
                coord3[[i, 2]] = s.coord[[mol_index[2][i], 2]];
            }
            for i in 0..mol_index[3].len()
            {
                coord4[[i, 0]] = s.coord[[mol_index[3][i], 0]];
                coord4[[i, 1]] = s.coord[[mol_index[3][i], 1]];
                coord4[[i, 2]] = s.coord[[mol_index[3][i], 2]];
            }

            // Move the geometric centers of the four molecules to the origin
            coord1 -= &coord1.mean_axis(Axis(0)).expect(&error_none_value("coord1"));
            coord2 -= &coord2.mean_axis(Axis(0)).expect(&error_none_value("coord2"));
            coord3 -= &coord3.mean_axis(Axis(0)).expect(&error_none_value("coord3"));
            coord4 -= &coord4.mean_axis(Axis(0)).expect(&error_none_value("coord4"));

            // Construct the rotational matrix
            let mut rot1: Array2<f64> = matrix::rand_rot();
            let mut rot2: Array2<f64> = matrix::rand_rot();
            let mut rot3: Array2<f64> = matrix::rand_rot();
            let mut rot4: Array2<f64> = matrix::rand_rot();
            let rot1_p: &mut [f64] = rot1.as_slice_mut().expect(&error_as_slice("rot1"));
            let rot2_p: &mut [f64] = rot2.as_slice_mut().expect(&error_as_slice("rot2"));
            let rot3_p: &mut [f64] = rot3.as_slice_mut().expect(&error_as_slice("rot3"));
            let rot4_p: &mut [f64] = rot4.as_slice_mut().expect(&error_as_slice("rot4"));
            root_process.broadcast_into(rot1_p);
            root_process.broadcast_into(rot2_p);
            root_process.broadcast_into(rot3_p);
            root_process.broadcast_into(rot4_p);

            // Rotate the four molecules
            coord1 = coord1.dot(&rot1);
            coord2 = coord2.dot(&rot2);
            coord3 = coord3.dot(&rot3);
            coord4 = coord4.dot(&rot4);

            // Construct the synthesis potential energy surface (PES)
            for i in 0..mol_index[0].len()
            {
                s.coord[[mol_index[0][i], 0]] = coord1[[i, 0]];
                s.coord[[mol_index[0][i], 1]] = coord1[[i, 1]];
                s.coord[[mol_index[0][i], 2]] = coord1[[i, 2]] + dist;
            }
            for i in 0..mol_index[1].len()
            {
                s.coord[[mol_index[1][i], 0]] = coord2[[i, 0]];
                s.coord[[mol_index[1][i], 1]] = coord2[[i, 1]] + dist * 0.942809041582063_f64;
                s.coord[[mol_index[1][i], 2]] = coord2[[i, 2]] - dist * 0.333333333333333_f64;
            }
            for i in 0..mol_index[2].len()
            {
                s.coord[[mol_index[2][i], 0]] = coord3[[i, 0]] + dist * 0.816496580927726_f64;
                s.coord[[mol_index[2][i], 1]] = coord3[[i, 1]] - dist * 0.471404520791032_f64;
                s.coord[[mol_index[2][i], 2]] = coord3[[i, 2]] - dist * 0.333333333333333_f64;
            }
            for i in 0..mol_index[3].len()
            {
                s.coord[[mol_index[3][i], 0]] = coord4[[i, 0]] - dist * 0.816496580927726_f64;
                s.coord[[mol_index[3][i], 1]] = coord4[[i, 1]] - dist * 0.471404520791032_f64;
                s.coord[[mol_index[3][i], 2]] = coord4[[i, 2]] - dist * 0.333333333333333_f64;
            }
        },

        _=>
        {
            panic!("The program for synthesis of more than 4 molecules has not finished.");
        },
    }
}










