//! About the pathway sampling
use std::fs::File;
use std::io::Write;
use crate::common::constants::ROOT_RANK;
use crate::common::error::*;
use crate::pes_exploration::potential::*;
use crate::pes_exploration::rtip::*;
use crate::pes_exploration::idwm::*;
use crate::pes_exploration::traits::*;
use crate::pes_exploration::system::System;
use crate::pes_exploration::optimization;
use mpi::traits::*;
use ndarray::{Array1, Array2, Axis, s};
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;





impl<'a> RtipPathSampling for RepulsivePot<'a>
{
    fn rtip_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P)
    {
        // Parallel processing
        let rank = comm.rank();
        let root_process = comm.process_at_rank(ROOT_RANK);

        // Output the local minimum structure to PDB output file, and output the header to RTIP output file
        if rank == ROOT_RANK
        {
            self.local_min.write_pdb(&self.str_output_file, true, 0);
            let mut rtip_output = File::create(&self.output_file).expect(&error_file("creating", &self.output_file));
            rtip_output.write_all(b"  step        rti_dist        pot_real        pot_rtip          f_real          f_rtip\n").expect(&error_file("writing", &self.output_file));
        }



        // Randomly obtain a nearby structure around self.local_min
        let mut s: System = self.local_min.clone();
        match &self.local_min.atom_add_pot
        {
            // If self.local_min.atom_add_pot has none value, disturb all the atoms 
            None =>
            {
                let mut dcoord: Array2<f64> = Array2::random(s.coord.raw_dim(), Uniform::new(0.0, 1.0));
                let dcoord_p: &mut [f64] = dcoord.as_slice_mut().expect(&error_as_slice("dcoord"));
                root_process.broadcast_into(dcoord_p);
                dcoord = &dcoord - dcoord.mean_axis(Axis(0)).unwrap();
                let dcoord_norm: f64 = (&dcoord * &dcoord).sum().sqrt();
                s.coord += &(dcoord * (0.1/dcoord_norm));
            },

            // If self.local_min.atom_add_pot specifies the atoms to be added RTI potential, only disturb these atoms
            Some(atom_add_pot) =>
            {
                let mut dcoord: Array2<f64> = Array2::random((atom_add_pot.len(), 3), Uniform::new(0.0, 1.0));
                let dcoord_p: &mut [f64] = dcoord.as_slice_mut().expect(&error_as_slice("dcoord"));
                root_process.broadcast_into(dcoord_p);
                dcoord = &dcoord - dcoord.mean_axis(Axis(0)).unwrap();
                let dcoord_norm: f64 = (&dcoord * &dcoord).sum().sqrt();
                for i in 0..atom_add_pot.len()
                {
                    s.coord[[atom_add_pot[i], 0]] += dcoord[[i, 0]] * (0.1/dcoord_norm);
                    s.coord[[atom_add_pot[i], 1]] += dcoord[[i, 1]] * (0.1/dcoord_norm);
                    s.coord[[atom_add_pot[i], 2]] += dcoord[[i, 2]] * (0.1/dcoord_norm);
                }
            },
        }



        // Define the variables for pathway sampling
        let mut pot_real: f64;
        let mut force_real: Array2<f64>;
        let mut f_real: f64;

        let mut local_min: System;                  // A fragment of the original local minimum containing the atoms to be added RTI potential
        let mut nearby_ts: Vec<System>;                 // A fragment of the original nearby TS containing the atoms to be added RTI potential
        let mut s_fragment: System;                 // A fragment of the current structure containing the atoms to be added RTI potential
        let mut rtip0_pes = match &self.local_min.atom_add_pot
        {
            // If self.local_min.atom_add_pot has none value, create the RTIP PES from the original local minimum and nearby TS
            None =>
            {
                // Initialize the fragments
                local_min = System
                {
                    natom: 0,
                    coord: Array2::zeros((0, 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };
                nearby_ts = vec![local_min.clone(); self.nearby_ts.len()];
                s_fragment = local_min.clone();

                // Create the RTI potential from the original local minimum and nearby TS
                Rtip0PES
                {
                    local_min: &self.local_min,
                    nearby_ts: &self.nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; self.nearby_ts.len()],
                }
            },

            // If self.local_min.atom_add_pot specifies the atoms to be added RTI potential, create fragments for the local minimum and nearby TS
            Some(atom_add_pot) =>
            {
                // Initialize the fragments
                local_min = System
                {
                    natom: atom_add_pot.len(),
                    coord: Array2::zeros((atom_add_pot.len(), 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };
                nearby_ts = vec![local_min.clone(); self.nearby_ts.len()];
                s_fragment = local_min.clone();

                // Copy the atomic coordinates from the original local minimum and nearby TS to the corresponding fragments
                for i in 0..atom_add_pot.len()
                {
                    local_min.coord[[i, 0]] = self.local_min.coord[[atom_add_pot[i], 0]];
                    local_min.coord[[i, 1]] = self.local_min.coord[[atom_add_pot[i], 1]];
                    local_min.coord[[i, 2]] = self.local_min.coord[[atom_add_pot[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot.len()
                    {
                        nearby_ts[j].coord[[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot[i], 0]];
                        nearby_ts[j].coord[[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot[i], 1]];
                        nearby_ts[j].coord[[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot[i], 2]];
                    }
                }

                // Create the RTIP PES from the fragments
                Rtip0PES
                {
                    local_min: &local_min,
                    nearby_ts: &nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; self.nearby_ts.len()],
                }
            },
        };

        let mut pot_rtip: f64;
        let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut force_rtip_fragment: Array2<f64>;
        let mut f_rtip: f64;

        let mut force_total: Array2<f64>;
        let mut dcoord: Array2<f64>;
        let mut pot_real_max: f64 = -1000000000000000.0;
        let mut pot_real_min: f64 = 1000000000000000.0;
        let mut add_rtip: bool = true;



        // Perform the pathway sampling iteractively
        for i in 1..(self.para.max_step+1)
        {
            // Calculate the real potential energy and atomic forces
            (pot_real, force_real) = real_pes.get_energy_force(&s);
            s.pot = pot_real;
            f_real = (&force_real * &force_real).sum().sqrt();

            // Calculate the RTIP potential energy and atomic forces
            match &self.local_min.atom_add_pot
            {
                // If self.local_min.atom_add_pot has none value, add the RTI potential on all the atoms
                None =>
                {
                    rtip0_pes.sigma_min = rti_dist(&self.local_min.coord, &s.coord);
                    if add_rtip
                    {
                        match self.para.scale_ts_sigma
                        {
                            None =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    rtip0_pes.sigma_ts[j] = rti_dist(&self.nearby_ts[j].coord, &s.coord);
                                }
                            },
                            Some(scale_ts_sigma) =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    rtip0_pes.sigma_ts[j] = (0.5 * scale_ts_sigma) * rti_dist(&self.nearby_ts[j].coord, &self.local_min.coord);
                                }
                            },
                        }
                        rtip0_pes.a_min = self.para.a0 * (i as f64);
                        rtip0_pes.a_ts = rtip0_pes.a_min * self.para.scale_ts_a0;
                        (pot_rtip, force_rtip) = rtip0_pes.get_energy_force(&s);
                        f_rtip = (&force_rtip * &force_rtip).sum().sqrt();
                    }
                    else
                    {
                        rtip0_pes.sigma_ts = vec![10.0; self.nearby_ts.len()];
                        rtip0_pes.a_min = 0.0;
                        rtip0_pes.a_ts = 0.0;
                        pot_rtip = 0.0;
                        force_rtip = Array2::zeros(s.coord.raw_dim());
                        f_rtip = 0.0;
                    }
                },

                // If self.local_min.atom_add_pot specifies the atoms to be add RTI potential, add the RTI potential on the fragments
                Some(atom_add_pot) =>
                {
                    for j in 0..atom_add_pot.len()
                    {
                        s_fragment.coord[[j, 0]] = s.coord[[atom_add_pot[j], 0]];
                        s_fragment.coord[[j, 1]] = s.coord[[atom_add_pot[j], 1]];
                        s_fragment.coord[[j, 2]] = s.coord[[atom_add_pot[j], 2]];
                    }
                    rtip0_pes.sigma_min = rti_dist(&local_min.coord, &s_fragment.coord);
                    if add_rtip
                    {
                        match self.para.scale_ts_sigma
                        {
                            None =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    rtip0_pes.sigma_ts[j] = rti_dist(&nearby_ts[j].coord, &s_fragment.coord);
                                }
                            },
                            Some(scale_ts_sigma) =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    rtip0_pes.sigma_ts[j] = (0.5 * scale_ts_sigma) * rti_dist(&nearby_ts[j].coord, &local_min.coord);
                                }
                            },
                        }
                        rtip0_pes.a_min = self.para.a0 * (i as f64);
                        rtip0_pes.a_ts = rtip0_pes.a_min * self.para.scale_ts_a0;
                        (pot_rtip, force_rtip_fragment) = rtip0_pes.get_energy_force(&s_fragment);
                        f_rtip = (&force_rtip_fragment * &force_rtip_fragment).sum().sqrt();
                        for j in 0..atom_add_pot.len()
                        {
                            force_rtip[[atom_add_pot[j], 0]] = force_rtip_fragment[[j, 0]];
                            force_rtip[[atom_add_pot[j], 1]] = force_rtip_fragment[[j, 1]];
                            force_rtip[[atom_add_pot[j], 2]] = force_rtip_fragment[[j, 2]];
                        }
                    }
                    else
                    {
                        rtip0_pes.sigma_ts = vec![10.0; self.nearby_ts.len()];
                        rtip0_pes.a_min = 0.0;
                        rtip0_pes.a_ts = 0.0;
                        pot_rtip = 0.0;
                        force_rtip = Array2::zeros(s.coord.raw_dim());
                        f_rtip = 0.0;
                    }
                },
            }

            // Output the information in this iterative step
            if rank == ROOT_RANK
            {
                if (i % self.para.print_step) == 0
                {
                    s.write_pdb(&self.str_output_file, false, i);
                }
                let mut rtip_output = File::options().append(true).open(&self.output_file).expect(&error_file("opening", &self.output_file));
                rtip_output.write_all(format!("{:6} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8}\n", i, rtip0_pes.sigma_min, pot_real, pot_rtip, f_real, f_rtip).as_bytes()).expect(&error_file("writing", &self.output_file));
            }

            // Search the structure of minimum potential energy along the given direction
            force_total = &force_real + &force_rtip;
            dcoord = optimization::min_1d_real_bias(real_pes, &rtip0_pes, &s, &force_total, pot_real+pot_rtip, self.para.pot_epsilon*(s.natom as f64));
            s.coord += &dcoord;

            // Judge whether to stop the loop or not
            if pot_real > pot_real_max
            {
                pot_real_max = pot_real;
            }
            if pot_real < pot_real_min
            {
                pot_real_min = pot_real;
            }
            if pot_real < (pot_real_max - self.para.pot_drop)
            {
                add_rtip = false;
            }
            if (add_rtip == false) && (f_real/(s.natom as f64).sqrt() < self.para.f_epsilon)
            {
                break
            }
            if pot_real > (pot_real_min + self.para.pot_climb) || f_rtip > 1000.0
            {
                break
            }
        }
    }
}





impl<'a> IdwmPathSampling for RepulsivePot<'a>
{
    fn idwm_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P)
    {
        // Parallel processing
        let rank = comm.rank();
        let root_process = comm.process_at_rank(ROOT_RANK);

        // Output the local minimum structure to PDB output file, and output the header to IDWM output file
        if rank == ROOT_RANK
        {
            self.local_min.write_pdb(&self.str_output_file, true, 0);
            let mut idwm_output = File::create(&self.output_file).expect(&error_file("creating", &self.output_file));
            idwm_output.write_all(b"  step        idw_dist        pot_real        pot_idwm          f_real          f_idwm\n").expect(&error_file("writing", &self.output_file));
        }



        // Randomly obtain a nearby structure around self.local_min
        let mut s: System = self.local_min.clone();
        match &self.local_min.atom_add_pot
        {
            // If self.local_min.atom_add_pot has none value, disturb all the atoms
            None =>
            {
                let mut dcoord: Array2<f64> = Array2::random(s.coord.raw_dim(), Uniform::new(0.0, 1.0));
                let dcoord_p: &mut [f64] = dcoord.as_slice_mut().expect(&error_as_slice("dcoord"));
                root_process.broadcast_into(dcoord_p);
                dcoord = &dcoord - dcoord.mean_axis(Axis(0)).unwrap();
                let dcoord_norm: f64 = (&dcoord * &dcoord).sum().sqrt();
                s.coord += &(dcoord * (0.1/dcoord_norm));
            },

            // If self.local_min.atom_add_pot specifies the atoms to be added IDW potential, only disturb these atoms
            Some(atom_add_pot) =>
            {
                let mut dcoord: Array2<f64> = Array2::random((atom_add_pot.len(), 3), Uniform::new(0.0, 1.0));
                let dcoord_p: &mut [f64] = dcoord.as_slice_mut().expect(&error_as_slice("dcoord"));
                root_process.broadcast_into(dcoord_p);
                dcoord = &dcoord - dcoord.mean_axis(Axis(0)).unwrap();
                let dcoord_norm: f64 = (&dcoord * &dcoord).sum().sqrt();
                for i in 0..atom_add_pot.len()
                {
                    s.coord[[atom_add_pot[i], 0]] += dcoord[[i, 0]] * (0.1/dcoord_norm);
                    s.coord[[atom_add_pot[i], 1]] += dcoord[[i, 1]] * (0.1/dcoord_norm);
                    s.coord[[atom_add_pot[i], 2]] += dcoord[[i, 2]] * (0.1/dcoord_norm);
                }
            },
        }



        // Define the variables for pathway sampling
        let mut pot_real: f64;
        let mut force_real: Array2<f64>;
        let mut f_real: f64;

        let mut s_fragment: System;             // A fragment of the current structure containing the atoms to be added IDW potential
        let mut idwm0_pes = match &self.local_min.atom_add_pot
        {
            // If self.local_min.atom_add_pot has none value, create the IDWM PES from the original local minimum and nearby TS
            None =>
            {
                // Initialize the fragents
                s_fragment = System
                {
                    natom: 0,
                    coord: Array2::zeros((0, 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };

                // Obtain the weighted distance matrix for local minimum and nearby TS
                let local_min: Array2<f64> = wei_dist_mat(&self.local_min.coord);
                let mut nearby_ts: Vec<Array2<f64>> = vec![Array2::zeros((self.local_min.natom, self.local_min.natom)); self.nearby_ts.len()];
                for i in 0..self.nearby_ts.len()
                {
                    nearby_ts[i] = wei_dist_mat(&self.nearby_ts[i].coord);
                }

                // Create the IDWM PES from the original local minimum and nearby TS
                Idwm0PES
                {
                    local_min,
                    nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; self.nearby_ts.len()],
                }
            },

            // If self.local_min.atom_add_pot specifies the atoms to be added IDWM potential, create fragments for the local minimum and nearby TS
            Some(atom_add_pot) =>
            {
                // Initialize the fragments
                s_fragment = System
                {
                    natom: atom_add_pot.len(),
                    coord: Array2::zeros((atom_add_pot.len(), 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };
                let mut local_min_coord: Array2<f64> = Array2::zeros((atom_add_pot.len(), 3));
                let mut nearby_ts_coord: Vec<Array2<f64>> = vec![local_min_coord.clone(); self.nearby_ts.len()];

                // Copy the atomic coordinates from the original local minimum and nearby TS to the corresponding fragments
                for i in 0..atom_add_pot.len()
                {
                    local_min_coord[[i, 0]] = self.local_min.coord[[atom_add_pot[i], 0]];
                    local_min_coord[[i, 1]] = self.local_min.coord[[atom_add_pot[i], 1]];
                    local_min_coord[[i, 2]] = self.local_min.coord[[atom_add_pot[i], 2]];
                }
                for j in 0..self.nearby_ts.len()
                {
                    for i in 0..atom_add_pot.len()
                    {
                        nearby_ts_coord[j][[i, 0]] = self.nearby_ts[j].coord[[atom_add_pot[i], 0]];
                        nearby_ts_coord[j][[i, 1]] = self.nearby_ts[j].coord[[atom_add_pot[i], 1]];
                        nearby_ts_coord[j][[i, 2]] = self.nearby_ts[j].coord[[atom_add_pot[i], 2]];
                    }
                }

                // Obtain the weighted distance matrix for local minimum and nearby TS
                let local_min: Array2<f64> = wei_dist_mat(&local_min_coord);
                let mut nearby_ts: Vec<Array2<f64>> = vec![Array2::zeros((atom_add_pot.len(), atom_add_pot.len())); self.nearby_ts.len()];
                for i in 0..self.nearby_ts.len()
                {
                    nearby_ts[i] = wei_dist_mat(&nearby_ts_coord[i]);
                }

                // Create the IDWM PES from the fragments
                Idwm0PES
                {
                    local_min,
                    nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; self.nearby_ts.len()],
                }
            },
        };

        let mut pot_idwm: f64;
        let mut force_idwm: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut force_idwm_fragment: Array2<f64>;
        let mut f_idwm: f64;

        let mut force_total: Array2<f64>;
        let mut dcoord: Array2<f64>;
        let mut pot_real_max: f64 = -1000000000000000.0;
        let mut pot_real_min: f64 = 1000000000000000.0;
        let mut add_idwm: bool = true;



        // Perform the pathway sampling iteractively
        for i in 1..(self.para.max_step+1)
        {
            // Calculate the real potential energy and atomic forces
            (pot_real, force_real) = real_pes.get_energy_force(&s);
            s.pot = pot_real;
            f_real = (&force_real * &force_real).sum().sqrt();

            // Calculate the IDWM potential energy and atomic forces
            match &self.local_min.atom_add_pot
            {
                // If self.local_min.atom_add_pot has none value, add the IDWM potential on all the atoms
                None =>
                {
                    idwm0_pes.sigma_min = idw_dist(&idwm0_pes.local_min, &s.coord);
                    if add_idwm
                    {
                        match self.para.scale_ts_sigma
                        {
                            None =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    idwm0_pes.sigma_ts[j] = idw_dist(&idwm0_pes.nearby_ts[j], &s.coord);
                                }
                            },
                            Some(scale_ts_sigma) =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    idwm0_pes.sigma_ts[j] = (0.5 * scale_ts_sigma) * idw_dist1(&idwm0_pes.nearby_ts[j], &idwm0_pes.local_min);
                                }
                            },
                        }
                        idwm0_pes.a_min = self.para.a0 * (i as f64);
                        idwm0_pes.a_ts = idwm0_pes.a_min * self.para.scale_ts_a0;
                        (pot_idwm, force_idwm) = idwm0_pes.get_energy_force(&s);
                        f_idwm = (&force_idwm * &force_idwm).sum().sqrt();
                    }
                    else
                    {
                        idwm0_pes.sigma_ts = vec![10.0; self.nearby_ts.len()];
                        idwm0_pes.a_min = 0.0;
                        idwm0_pes.a_ts = 0.0;
                        pot_idwm = 0.0;
                        force_idwm = Array2::zeros(s.coord.raw_dim());
                        f_idwm = 0.0;
                    }
                },

                // If self.local_min.atom_add_pot specifies the atoms to be added IDWM potential, add the IDWM potential on the fragments
                Some(atom_add_pot) =>
                {
                    for j in 0..atom_add_pot.len()
                    {
                        s_fragment.coord[[j, 0]] = s.coord[[atom_add_pot[j], 0]];
                        s_fragment.coord[[j, 1]] = s.coord[[atom_add_pot[j], 1]];
                        s_fragment.coord[[j, 2]] = s.coord[[atom_add_pot[j], 2]];
                    }
                    idwm0_pes.sigma_min = idw_dist(&idwm0_pes.local_min, &s_fragment.coord);
                    if add_idwm
                    {
                        match self.para.scale_ts_sigma
                        {
                            None =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    idwm0_pes.sigma_ts[j] = idw_dist(&idwm0_pes.nearby_ts[j], &s_fragment.coord);
                                }
                            },
                            Some(scale_ts_sigma) =>
                            {
                                for j in 0..self.nearby_ts.len()
                                {
                                    idwm0_pes.sigma_ts[j] = (0.5 * scale_ts_sigma) * idw_dist1(&idwm0_pes.nearby_ts[j], &idwm0_pes.local_min);
                                }
                            },
                        }
                        idwm0_pes.a_min = self.para.a0 * (i as f64);
                        idwm0_pes.a_ts = idwm0_pes.a_min * self.para.scale_ts_a0;
                        (pot_idwm, force_idwm_fragment) = idwm0_pes.get_energy_force(&s_fragment);
                        f_idwm = (&force_idwm_fragment * &force_idwm_fragment).sum().sqrt();
                        for j in 0..atom_add_pot.len()
                        {
                            force_idwm[[atom_add_pot[j], 0]] = force_idwm_fragment[[j, 0]];
                            force_idwm[[atom_add_pot[j], 1]] = force_idwm_fragment[[j, 1]];
                            force_idwm[[atom_add_pot[j], 2]] = force_idwm_fragment[[j, 2]];
                        }
                    }
                    else
                    {
                        idwm0_pes.sigma_ts = vec![10.0; self.nearby_ts.len()];
                        idwm0_pes.a_min = 0.0;
                        idwm0_pes.a_ts = 0.0;
                        pot_idwm = 0.0;
                        force_idwm = Array2::zeros(s.coord.raw_dim());
                        f_idwm = 0.0;
                    }
                },
            }

            // Output the information in this iterative step
            if rank == ROOT_RANK
            {
                if (i % self.para.print_step) == 0
                {
                    s.write_pdb(&self.str_output_file, false, i);
                }
                let mut idwm_output = File::options().append(true).open(&self.output_file).expect(&error_file("opening", &self.output_file));
                idwm_output.write_all(format!("{:6} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8}\n", i, idwm0_pes.sigma_min, pot_real, pot_idwm, f_real, f_idwm).as_bytes()).expect(&error_file("writing", &self.output_file));
            }

            // Search the structure of minimum potential energy along the given direction
            force_total = &force_real + &force_idwm;
            dcoord = optimization::min_1d_real_bias(real_pes, &idwm0_pes, &s, &force_total, pot_real+pot_idwm, self.para.pot_epsilon*(s.natom as f64));
            s.coord += &dcoord;

            // Judge whether to stop the loop or not
            if pot_real > pot_real_max
            {
                pot_real_max = pot_real;
            }
            if pot_real < pot_real_min
            {
                pot_real_min = pot_real;
            }
            if pot_real < (pot_real_max - self.para.pot_drop)
            {
                add_idwm = false;
            }
            if (add_idwm == false) && (f_real/(s.natom as f64).sqrt() < self.para.f_epsilon)
            {
                break
            }
            if pot_real > (pot_real_min + self.para.pot_climb) || f_idwm > 1000.0
            {
                break
            }
        }
    }
}





impl<'a> RtipPathSampling for AttractivePot<'a>
{
    fn rtip_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P)
    {
        // Parallel processing
        let rank = comm.rank();

        // Output the initial state structure to PDB output file, and output the header to RTIP output file
        if rank == ROOT_RANK
        {
            self.initial_state.write_pdb(&self.str_output_file, true, 0);
            let mut rtip_output = File::create(&self.output_file).expect(&error_file("creating", &self.output_file));
            rtip_output.write_all(b"  step        rti_dist        pot_real        pot_rtip          f_real          f_rtip\n").expect(&error_file("writing", &self.output_file));
        }

        // Begin from the initial state structure
        let mut s: System = self.initial_state.clone();



        // Define the variables for pathway sampling
        let mut pot_real: f64;
        let mut force_real: Array2<f64>;
        let mut f_real: f64;

        let mut final_state: System;            // A fragment of the original final state structure containing the atoms to be added RTI potential
        let mut _nearby_ts: Vec<System>;
        let mut s_fragment: System;             // A fragment of the current structure containing the atoms to be added RTI potential
        let mut rtip0_pes = match &self.final_state.atom_add_pot
        {
            // If self.final_state.atom_add_pot has none value, create the RTIP PES from the original final state structure
            None =>
            {
                // Initialize the fragments
                final_state = System
                {
                    natom: 0,
                    coord: Array2::zeros((0, 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };
                _nearby_ts = vec![final_state.clone(); 0];
                s_fragment = final_state.clone();

                // Create the RTI potential from the original final state
                Rtip0PES
                {
                    local_min: &self.final_state,
                    nearby_ts: &_nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; 0],
                }
            },

            // If self.final_state.atom_add_pot specifies the atoms to be added RTI potential, create fragments for the final state structure
            Some(atom_add_pot) =>
            {
                // Initialize the fragments
                final_state = System
                {
                    natom: atom_add_pot.len(),
                    coord: Array2::zeros((atom_add_pot.len(), 3)),
                    cell: None,
                    atom_type: None,
                    atom_add_pot: None,
                    mutable: None,
                    pot: 0.0,
                };
                _nearby_ts = vec![final_state.clone(); 0];
                s_fragment = final_state.clone();

                // Copy the atomic coordinates from the original final state structure to its fragment
                for i in 0..atom_add_pot.len()
                {
                    final_state.coord[[i, 0]] = self.final_state.coord[[atom_add_pot[i], 0]];
                    final_state.coord[[i, 1]] = self.final_state.coord[[atom_add_pot[i], 1]];
                    final_state.coord[[i, 2]] = self.final_state.coord[[atom_add_pot[i], 2]];
                }

                // Create the RTIP PES from the fragment
                Rtip0PES
                {
                    local_min: &final_state,
                    nearby_ts: &_nearby_ts,
                    a_min: 0.0,
                    a_ts: 0.0,
                    sigma_min: 10.0,
                    sigma_ts: vec![10.0; 0],
                }
            },
        };

        let mut pot_rtip: f64;
        let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut force_rtip_fragment: Array2<f64>;
        let mut f_rtip: f64;

        let mut force_total: Array2<f64>;
        let mut dcoord: Array2<f64>;
        let mut add_rtip: bool = true;



        // Perform the pathway sampling iteractively
        for i in 1..(self.para.max_step+1)
        {
            // Calculate the real potential energy and atomic forces
            (pot_real, force_real) = real_pes.get_energy_force(&s);
            s.pot = pot_real;
            f_real = (&force_real * &force_real).sum().sqrt();

            // Calculate the RTIP potential energy and atomic forces
            match &self.final_state.atom_add_pot
            {
                // If self.final_state.atom_add_pot has none value, add the RTI potential on all the atoms
                None =>
                {
                    rtip0_pes.sigma_min = rti_dist(&self.final_state.coord, &s.coord);
                    if add_rtip
                    {
                        rtip0_pes.a_min = -self.para.a0 * (i as f64);
                        (pot_rtip, force_rtip) = rtip0_pes.get_energy_force(&s);
                        f_rtip = (&force_rtip * &force_rtip).sum().sqrt();
                    }
                    else
                    {
                        rtip0_pes.a_min = 0.0;
                        pot_rtip = 0.0;
                        force_rtip = Array2::zeros(s.coord.raw_dim());
                        f_rtip = 0.0;
                    }
                },

                // If self.final_state.atom_add_pot specifies the atoms to be added RTI potential, add the RTI potential on the fragments
                Some(atom_add_pot) =>
                {
                    for j in 0..atom_add_pot.len()
                    {
                        s_fragment.coord[[j, 0]] = s.coord[[atom_add_pot[j], 0]];
                        s_fragment.coord[[j, 1]] = s.coord[[atom_add_pot[j], 1]];
                        s_fragment.coord[[j, 2]] = s.coord[[atom_add_pot[j], 2]];
                    }
                    rtip0_pes.sigma_min = rti_dist(&final_state.coord, &s_fragment.coord);
                    if add_rtip
                    {
                        rtip0_pes.a_min = -self.para.a0 * (i as f64);
                        (pot_rtip, force_rtip_fragment) = rtip0_pes.get_energy_force(&s_fragment);
                        f_rtip = (&force_rtip_fragment * &force_rtip_fragment).sum().sqrt();
                        for j in 0..atom_add_pot.len()
                        {
                            force_rtip[[atom_add_pot[j], 0]] = force_rtip_fragment[[j, 0]];
                            force_rtip[[atom_add_pot[j], 1]] = force_rtip_fragment[[j, 1]];
                            force_rtip[[atom_add_pot[j], 2]] = force_rtip_fragment[[j, 2]];
                        }
                    }
                    else
                    {
                        rtip0_pes.a_min = 0.0;
                        pot_rtip = 0.0;
                        force_rtip = Array2::zeros(s.coord.raw_dim());
                        f_rtip = 0.0;
                    }
                },
            }

            // Output the information in this iteractive step
            if rank == ROOT_RANK
            {
                if (i % self.para.print_step) == 0
                {
                    s.write_pdb(&self.str_output_file, false, i);
                }
                let mut rtip_output = File::options().append(true).open(&self.output_file).expect(&error_file("opening", &self.output_file));
                rtip_output.write_all(format!("{:6} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8}\n", i, rtip0_pes.sigma_min, pot_real, pot_rtip, f_real, f_rtip).as_bytes()).expect(&error_file("writing", &self.output_file));
            }

            // Search the structure of minimum potential energy along the given direction
            force_total = &force_real + &force_rtip;
            dcoord = optimization::min_1d_real_bias(real_pes, &rtip0_pes, &s, &force_total, pot_real+pot_rtip, self.para.pot_epsilon*(s.natom as f64));
            s.coord += &dcoord;

            // Judge whether to stop the loop or not
            if rtip0_pes.sigma_min < 1.0 || f_rtip > 1000.0
            {
                add_rtip = false;
            }
            if (add_rtip == false) && (f_real/(s.natom as f64).sqrt() < self.para.f_epsilon)
            {
                break
            }
        }
    }
}





impl<'a> RtipPathSampling for SynthesisPot<'a>
{
    fn rtip_path_sampling<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P)
    {
        // Parallel processing
        let rank = comm.rank();

        // Output the initial state structure to PDB output file, and output the header to RTIP output file
        if rank == ROOT_RANK
        {
            self.initial_state.write_pdb(&self.str_output_file, true, 0);
            let mut rtip_output = File::create(&self.output_file).expect(&error_file("creating", &self.output_file));
            rtip_output.write_all(b"  step        rti_dist        pot_real        pot_rtip          f_real          f_rtip\n").expect(&error_file("writing", &self.output_file));
        }

        // Begin from the initial state structure
        let mut s: System = self.initial_state.clone();



        // Define the variables for pathway sampling
        let mut pot_real: f64;
        let mut force_real: Array2<f64>;
        let mut f_real: f64;

        let mut mol_indices: Vec<usize> = Vec::new();
        for i in 0..self.mol_index.len()
        {
            mol_indices.extend_from_slice(&self.mol_index[i]);
        }
        let mut o_all: Array1<f64>;                 // To reserve the geometric center of all the molecules for synthesis
        let mut o_i: Array1<f64>;               // To reserve the geometric center of a special molecule for synthesis
        let mut index: usize;

        let mut final_state: System = System                // Initialization of the final state of the molecules for synthesis, where their geometric centers coincide
        {
            natom: mol_indices.len(),
            coord: Array2::zeros((mol_indices.len(), 3)),
            cell: None,
            atom_type: None,
            atom_add_pot: None,
            mutable: None,
            pot: 0.0,
        };
        let _nearby_ts: Vec<System> = vec![final_state.clone(); 0];
        let mut s_fragment: System =
        if mol_indices.len() == s.coord.len_of(Axis(0))             // If the molecules for synthesis construct the whole system (i.e. without the other environment atoms or molecules)
        {
            s.atom_add_pot = None;
            System
            {
                natom: 0,
                coord: Array2::zeros((0, 3)),
                cell: None,
                atom_type: None,
                atom_add_pot: None,
                mutable: None,
                pot: 0.0,
            }
        }
        else                // If besides the molecules for synthesis, there is other environment atoms or moleculs
        {
            s.atom_add_pot = Some(mol_indices.clone());
            final_state.clone()
        };

        let mut pot_rtip: f64;
        let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut force_rtip_fragment: Array2<f64>;
        let mut f_rtip: f64;

        let mut force_total: Array2<f64>;
        let mut dcoord: Array2<f64>;
        let mut pot_real_max: f64 = -1000000000000000.0;
        let mut pot_real_min: f64 = 1000000000000000.0;
        let mut pot_real_initial: f64 = 0.0;
        let mut add_rtip: bool = true;

        // Perform the pathway sampling iteractively
        for i in 1..(self.para.max_step+1)
        {
            // Calculate the real potential energy and atomic forces
            (pot_real, force_real) = real_pes.get_energy_force(&s);
            s.pot = pot_real;
            f_real = (&force_real * &force_real).sum().sqrt();

            // Calculate the RTIP potential energy and atomic forces
            let rtip0_pes: Rtip0PES =
            if mol_indices.len() == s.coord.len_of(Axis(0))             // If the molecules for synthesis construct the whole system (i.e. without the other environment atoms or molecules)
            {
                // Define the final state of the molecules for synthesis, where their geometric centers coincide
                o_all = s.coord.mean_axis(Axis(0)).expect(&error_none_value("s.coord"));
                for j in 0..self.mol_index.len()                // For each molecule in synthesis
                {
                    // Find the geometric center of Molecule j
                    o_i = Array1::zeros(3);
                    for k in 0..self.mol_index[j].len()                 // for each atom in Molecule j
                    {
                        o_i += &s.coord.slice(s![self.mol_index[j][k], ..]);
                    }
                    o_i /= self.mol_index[j].len() as f64;
                    // Move Molecule j for the final state
                    for k in 0..self.mol_index[j].len()
                    {
                        index = self.mol_index[j][k];
                        final_state.coord[[index, 0]] = s.coord[[index, 0]] - o_i[0] + o_all[0];
                        final_state.coord[[index, 1]] = s.coord[[index, 1]] - o_i[1] + o_all[1];
                        final_state.coord[[index, 2]] = s.coord[[index, 2]] - o_i[2] + o_all[2];
                    }
                }
                // Construct the synthesis potential energy surface
                if add_rtip
                {
                    let rtip0_pes = Rtip0PES
                    {
                        local_min: &final_state,
                        nearby_ts: &_nearby_ts,
                        a_min: -self.para.a0 * (i as f64),
                        a_ts: 0.0,
                        sigma_min: rti_dist(&final_state.coord, &s.coord),
                        sigma_ts: vec![10.0; 0],
                    };
                    (pot_rtip, force_rtip) = rtip0_pes.get_energy_force(&s);
                    f_rtip = (&force_rtip * &force_rtip).sum().sqrt();

                    rtip0_pes
                }
                else
                {
                    let rtip0_pes = Rtip0PES
                    {
                        local_min: &final_state,
                        nearby_ts: &_nearby_ts,
                        a_min: 0.0,
                        a_ts: 0.0,
                        sigma_min: rti_dist(&final_state.coord, &s.coord),
                        sigma_ts: vec![10.0; 0],
                    };
                    pot_rtip = 0.0;
                    force_rtip = Array2::zeros(s.coord.raw_dim());
                    f_rtip = 0.0;

                    rtip0_pes
                }
            }

            else                // If besides the molecules for synthesis, there is other environment atoms or moleculs
            {
                // Extract the atomic coordinates of the molecules for synthesis
                for j in 0..mol_indices.len()
                {
                    s_fragment.coord[[j, 0]] = s.coord[[mol_indices[j], 0]];
                    s_fragment.coord[[j, 1]] = s.coord[[mol_indices[j], 1]];
                    s_fragment.coord[[j, 2]] = s.coord[[mol_indices[j], 2]];
                }
                // Define the final state of the molecules for synthesis, where their geometric centers coincide
                o_all = s_fragment.coord.mean_axis(Axis(0)).expect(&error_none_value("s_fragment.coord"));
                index = 0;
                for j in 0..self.mol_index.len()                // For each molecule in synthesis
                {
                    // Find the geometric center of Molecule j
                    o_i = Array1::zeros(3);
                    for k in 0..self.mol_index[j].len()                 // for each atom in Molecule j
                    {
                        o_i += &s.coord.slice(s![self.mol_index[j][k], ..]);
                    }
                    o_i /= self.mol_index[j].len() as f64;
                    // Move Molecule j for the final state
                    for k in 0..self.mol_index[j].len()
                    {
                        final_state.coord[[index, 0]] = s.coord[[self.mol_index[j][k], 0]] - o_i[0] + o_all[0];
                        final_state.coord[[index, 1]] = s.coord[[self.mol_index[j][k], 1]] - o_i[1] + o_all[1];
                        final_state.coord[[index, 2]] = s.coord[[self.mol_index[j][k], 2]] - o_i[2] + o_all[2];
                        index += 1;
                    }
                }
                // Construct the synthesis potential energy surface
                if add_rtip
                {
                    let rtip0_pes = Rtip0PES
                    {
                        local_min: &final_state,
                        nearby_ts: &_nearby_ts,
                        a_min: -self.para.a0 * (i as f64),
                        a_ts: 0.0,
                        sigma_min: rti_dist(&final_state.coord, &s_fragment.coord),
                        sigma_ts: vec![10.0; 0],
                    };
                    (pot_rtip, force_rtip_fragment) = rtip0_pes.get_energy_force(&s_fragment);
                    f_rtip = (&force_rtip_fragment * &force_rtip_fragment).sum().sqrt();
                    for j in 0..mol_indices.len()
                    {
                        force_rtip[[mol_indices[j], 0]] = force_rtip_fragment[[j, 0]];
                        force_rtip[[mol_indices[j], 1]] = force_rtip_fragment[[j, 1]];
                        force_rtip[[mol_indices[j], 2]] = force_rtip_fragment[[j, 2]];
                    }

                    rtip0_pes
                }
                else
                {
                    let rtip0_pes = Rtip0PES
                    {
                        local_min: &final_state,
                        nearby_ts: &_nearby_ts,
                        a_min: 0.0,
                        a_ts: 0.0,
                        sigma_min: rti_dist(&final_state.coord, &s_fragment.coord),
                        sigma_ts: vec![10.0; 0],
                    };
                    pot_rtip = 0.0;
                    force_rtip = Array2::zeros(s.coord.raw_dim());
                    f_rtip = 0.0;

                    rtip0_pes
                }
            };

            // Output the information in this iteractive step
            if rank == ROOT_RANK
            {
                if (i % self.para.print_step) == 0
                {
                    s.write_pdb(&self.str_output_file, false, i);
                }
                let mut rtip_output = File::options().append(true).open(&self.output_file).expect(&error_file("opening", &self.output_file));
                rtip_output.write_all(format!("{:6} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8}\n", i, rtip0_pes.sigma_min, pot_real, pot_rtip, f_real, f_rtip).as_bytes()).expect(&error_file("writing", &self.output_file));
            }

            // Search the structure of minimum potential energy along the given direction
            force_total = &force_real + &force_rtip;
            dcoord = optimization::min_1d_real_bias(real_pes, &rtip0_pes, &s, &force_total, pot_real+pot_rtip, self.para.pot_epsilon*(s.natom as f64));
            s.coord += &dcoord;

            // Judge whether to stop the loop or not
            if pot_real > pot_real_max
            {
                pot_real_max = pot_real;
            }
            if pot_real < pot_real_min
            {
                pot_real_min = pot_real;
            }
            if i == 1
            {
                pot_real_initial = pot_real;
            }
            if pot_real < (pot_real_max - self.para.pot_drop) && pot_real > pot_real_initial
            {
                add_rtip = false;
            }
            if (add_rtip == false) && (f_real/(s.natom as f64).sqrt() < self.para.f_epsilon)
            {
                break
            }
            if pot_real > (pot_real_min + self.para.pot_climb) || f_rtip > 1000.0
            {
                break
            }
        }
    }
}










