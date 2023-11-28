use std::fs::File;
use std::io::Write;
use crate::common::constants::*;
use crate::common::error::*;
use crate::pes_exploration::potential::*;
use crate::pes_exploration::rtip::*;
use crate::pes_exploration::traits::*;
use crate::pes_exploration::system::System;
use mpi::traits::*;
use ndarray::{Array2, Axis};
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;





impl<'a> RtipNVTMD for RepulsivePot<'a>
{
    fn rtip_nvt_md<C: Communicator, P: PES>(&self, comm: &C, real_pes: &P)
    {
        // Parallel processing
        let rank = comm.rank();
        let root_process = comm.process_at_rank(ROOT_RANK);

        // Output the local minimum structure to PDB output file, and output the header to RTIP output file
        if rank == ROOT_RANK
        {
            self.local_min.write_pdb(&self.str_output_file, true, 0);
            let mut rtip_output = File::create(&self.output_file).expect(&error_file("creating", &self.output_file));
            rtip_output.write_all(b"  step            time        rti_dist            temp             kin        pot_real        pot_rtip          f_real          f_rtip\n").expect(&error_file("writing", &self.output_file));
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



        // Define the variables for NVT MD
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
                    sigma_min: 0.0,
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
                    sigma_min: 0.0,
                    sigma_ts: vec![10.0; self.nearby_ts.len()],
                }
            },
        };

        let mut pot_rtip: f64;
        let mut force_rtip: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut force_rtip_fragment: Array2<f64>;
        let mut f_rtip: f64;

        let mut force_total: Array2<f64>;
        let dt: f64 = self.para.dt;
        let mut t: f64 = 0.0;
        let mut kin: f64;
        let mut temp: f64;                      // K
        let mut lambda: f64;
        let mut vel: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut acc: Array2<f64> = Array2::zeros(s.coord.raw_dim());
        let mut atom_mass: Vec<f64> = vec![0.0; s.natom];
        for i in 0..s.natom
        {
            atom_mass[i] = s.atom_type.as_ref().expect(&error_none_value("s.atom_type"))[i].get_mass();
        }



        // Perform the NVT MD iteractively
        for i in 1..(self.para.max_step+1)
        {
            // First step of leapfrog method
            t += dt;                    // fs
            vel += &(0.5 * dt * FEMTOSECOND_TO_AU * &acc);                  // A.U.
            s.coord += &(dt * FEMTOSECOND_TO_AU * &vel);                  // A.U.

            // Berendsen thermostat
            kin = 0.0;
            for j in 0..s.natom
            {
                kin += atom_mass[j] * ( vel[[j, 0]].powi(2) + vel[[j, 1]].powi(2) + vel[[j, 2]].powi(2) );          // Hartree
            }
            temp = kin * HARTREE_TO_JOULE / ( BOLTZMANN * 3.0 * (s.natom-1) as f64 );               // K
            if temp < 1.0
            {
                temp = 1.0;
            }
            lambda = ( 1.0 + (dt / self.para.tau) * (self.para.temp_bath / temp - 1.0) ).sqrt();

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
                },
            }

            // Calculate the atomic acceleration from the total force
            force_total = &force_real + &force_rtip;
            for j in 0..s.natom
            {
                acc[[j, 0]] = force_total[[j, 0]] / atom_mass[j];
                acc[[j, 1]] = force_total[[j, 1]] / atom_mass[j];
                acc[[j, 2]] = force_total[[j, 2]] / atom_mass[j];
            }

            // Second step of leapfrog method
            vel += &(0.5 * dt * FEMTOSECOND_TO_AU * &acc);              // A.U.
            vel *= lambda;

            // Calculate kinetic energy
            kin = 0.0;
            for j in 0..s.natom
            {
                kin += atom_mass[j] * ( vel[[j, 0]].powi(2) + vel[[j, 1]].powi(2) + vel[[j, 2]].powi(2) );              // Hartree
            }
            temp = kin * HARTREE_TO_JOULE / ( BOLTZMANN * 3.0 * (s.natom-1) as f64 );               // K
            kin *= 0.5;                 // Hartree

            // Output the information in this iterative step
            if rank == ROOT_RANK
            {
                if (i % self.para.print_step) == 0
                {
                    s.write_pdb(&self.str_output_file, false, i);
                }
                let mut rtip_output = File::options().append(true).open(&self.output_file).expect(&error_file("opening", &self.output_file));
                rtip_output.write_all(format!("{:6} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8} {:15.8}\n", i, t, rtip0_pes.sigma_min, temp, kin, pot_real, pot_rtip, f_real, f_rtip).as_bytes()).expect(&error_file("writing", &self.output_file));
            }
        }
    }
}










