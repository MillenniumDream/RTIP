//! Interatomic distance weighted metric
use crate::pes_exploration::system::System;
use crate::pes_exploration::traits::PES;
use ndarray::{Array1, Array2, array, Axis, s};





fn w(x: f64) -> f64
{
    let n: i32 = 5;
    let sigma: f64 = 3.0;

    ( -(x/sigma).powi(n) ).exp() + 1.0
}

fn dw(x: f64) -> f64
{
    let n: i32 = 5;
    let sigma: f64 = 3.0;

    ( -(x/sigma).powi(n) ).exp() * ( -(n as f64) * x.powi(n-1) / sigma.powi(n) )
}





/// Input the atomic coordinates of a system (coord),
/// output the weighted distance matrix.
///
/// # Parameters
/// ```
/// coord: the atomic coordinates of the system (natom * 3)
/// w_mat: the weighted distance matrix (natom * natom)
/// ```
pub fn wei_dist_mat(coord: &Array2<f64>) -> Array2<f64>
{
    let natom: usize = coord.len_of(Axis(0));

    let mut w_mat: Array2<f64> = Array2::zeros((natom, natom));
    for i in 0..(natom-1)
    {
        for j in (i+1)..natom
        {
            w_mat[[i,j]] = w( ( (coord[[i,0]]-coord[[j,0]]).powi(2) + (coord[[i,1]]-coord[[j,1]]).powi(2) + (coord[[i,2]]-coord[[j,2]]).powi(2) ).sqrt() );
        }
    }

    return w_mat
}

/// Input the weighted distance matrix of a reference system (w_mat),
/// and the atomic coordinates of another system (coord),
/// output the interatomic distance weighted (IDW) metric between them.
///
/// # Parameters
/// ```
/// w_mat0: the weighted distance matrix of the reference system (natom * natom)
/// coord: the atomic coordinates of the system (natom * 3)
/// dist: the output IDW distance between the two systems
/// ```
pub fn idw_dist(w_mat0: &Array2<f64>, coord: &Array2<f64>) -> f64
{
    let natom: usize = coord.len_of(Axis(0));

    let mut dist: f64 = 0.0;
    for i in 0..(natom-1)
    {
        for j in (i+1)..natom
        {
            dist += (w( ( (coord[[i,0]]-coord[[j,0]]).powi(2) + (coord[[i,1]]-coord[[j,1]]).powi(2) + (coord[[i,2]]-coord[[j,2]]).powi(2) ).sqrt() ) - w_mat0[[i,j]]).powi(2);
        }
    }

    return dist.sqrt()
}

pub fn idw_dist1(w_mat0: &Array2<f64>, w_mat1: &Array2<f64>) -> f64
{
    let natom: usize = w_mat0.len_of(Axis(0));

    let mut dist: f64 = 0.0;
    for i in 0..(natom-1)
    {
        for j in (i+1)..natom
        {
            dist += (w_mat1[[i,j]] - w_mat0[[i,j]]).powi(2);
        }
    }

    return dist.sqrt()
}

pub fn idw_dist_vec(w_mat0: &Array2<f64>, coord: &Array2<f64>) -> (f64, Array2<f64>)
{
    let natom: usize = coord.len_of(Axis(0));

    let mut vec_mat: Array2< Array1<f64> > = Array2::from_elem((natom, natom), array![0.0, 0.0, 0.0]);
    let mut dist_mat: Array2<f64> = Array2::zeros((natom, natom));
    let mut w_mat: Array2<f64> = Array2::zeros((natom, natom));
    let mut dw_mat: Array2<f64> = Array2::zeros((natom, natom));
    let mut dist: f64 = 0.0;
    for i in 0..(natom-1)
    {
        for j in (i+1)..natom
        {
            vec_mat[[i,j]] = &coord.slice(s![i, ..]) - &coord.slice(s![j, ..]);            // the (i, j) element is a vector from atom j to atom i
            dist_mat[[i,j]] = (vec_mat[[i,j]][0].powi(2) + vec_mat[[i,j]][1].powi(2) + vec_mat[[i,j]][2].powi(2)).sqrt();           // the (i, j) element is the atomic distance between atom i and j
            w_mat[[i,j]] = w(dist_mat[[i,j]]);
            dw_mat[[i,j]] = dw(dist_mat[[i,j]]);
            dist += (w_mat[[i,j]] - w_mat0[[i,j]]).powi(2);
        }
    }
    dist = dist.sqrt();

    let mut vec: Array2<f64> = Array2::zeros((natom, 3));
    for i in 0..natom
    {
        for j in 0..natom
        {
            if i<j
            {
                vec[[i,0]] += (w_mat[[i,j]] - w_mat0[[i,j]]) * dw_mat[[i,j]] * vec_mat[[i,j]][0] / dist_mat[[i,j]];
                vec[[i,1]] += (w_mat[[i,j]] - w_mat0[[i,j]]) * dw_mat[[i,j]] * vec_mat[[i,j]][1] / dist_mat[[i,j]];
                vec[[i,2]] += (w_mat[[i,j]] - w_mat0[[i,j]]) * dw_mat[[i,j]] * vec_mat[[i,j]][2] / dist_mat[[i,j]];
            }
            if j<i
            {
                vec[[i,0]] -= (w_mat[[j,i]] - w_mat0[[j,i]]) * dw_mat[[j,i]] * vec_mat[[j,i]][0] / dist_mat[[j,i]];
                vec[[i,1]] -= (w_mat[[j,i]] - w_mat0[[j,i]]) * dw_mat[[j,i]] * vec_mat[[j,i]][1] / dist_mat[[j,i]];
                vec[[i,2]] -= (w_mat[[j,i]] - w_mat0[[j,i]]) * dw_mat[[j,i]] * vec_mat[[j,i]][2] / dist_mat[[j,i]];
            }
        }
    }
    vec /= dist;

    return (dist, vec)
}





fn idw_pot(w_mat0: &Array2<f64>, coord: &Array2<f64>, a: f64, sigma: f64) -> f64
{
    let idw_dist: f64 = idw_dist(w_mat0, coord);
    let pot_idw: f64 = a * ( -(idw_dist*idw_dist) / (2.0*sigma*sigma) ).exp();

    return pot_idw
}

fn idw_pot_force(w_mat0: &Array2<f64>, coord: &Array2<f64>, a: f64, sigma: f64) -> (f64, Array2<f64>)
{
    let (idw_dist, idw_vec) = idw_dist_vec(w_mat0, coord);
    let pot_idw: f64 = a * ( -(idw_dist*idw_dist) / (2.0*sigma*sigma) ).exp();
    let force_idw: Array2<f64> = idw_vec * pot_idw * idw_dist / (sigma*sigma);

    return (pot_idw, force_idw)
}





pub struct Idwm0PES
{
    pub local_min: Array2<f64>,
    pub nearby_ts: Vec<Array2<f64>>,
    pub a_min: f64,
    pub a_ts: f64,
    pub sigma_min: f64,
    pub sigma_ts: Vec<f64>,
}

impl PES for Idwm0PES
{
    fn get_energy(&self, s: &System) -> f64
    {
        match &s.atom_add_pot
        {
            // If the input structure do not specify the atoms to be added IDWM potential
            None =>
            {
                let mut pot_idwm: f64 = idw_pot(&self.local_min, &s.coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_idwm += idw_pot(&self.nearby_ts[i], &s.coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_idwm
            },

            // If the input structure specifies the atoms to be added IDWM potential
            Some(atom_add_pot) =>
            {
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot.len(), 3));
                for i in 0..atom_add_pot.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot[i], 2]];
                }

                let mut pot_idwm: f64 = idw_pot(&self.local_min, &s_fragment_coord, self.a_min, self.sigma_min);
                for i in 0..self.nearby_ts.len()
                {
                    pot_idwm += idw_pot(&self.nearby_ts[i], &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                }

                return pot_idwm
            },
        }
    }

    fn get_energy_force(&self, s: &System) -> (f64, Array2<f64>)
    {
        match &s.atom_add_pot
        {
            // If the input structure do not specify the atoms to be added IDWM potential
            None =>
            {
                let (mut pot_idwm, mut force_idwm) = idw_pot_force(&self.local_min, &s.coord, self.a_min, self.sigma_min);
                let mut pot_idwm_ts: f64;
                let mut force_idwm_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_idwm_ts, force_idwm_ts) = idw_pot_force(&self.nearby_ts[i], &s.coord, self.a_ts, self.sigma_ts[i]);
                    pot_idwm += pot_idwm_ts;
                    force_idwm += &force_idwm_ts;
                }

                return (pot_idwm, force_idwm)
            },

            // If the input structure specifies the atoms to be added IDWM potential
            Some(atom_add_pot) =>
            {
                let mut s_fragment_coord: Array2<f64> = Array2::zeros((atom_add_pot.len(), 3));
                for i in 0..atom_add_pot.len()
                {
                    s_fragment_coord[[i, 0]] = s.coord[[atom_add_pot[i], 0]];
                    s_fragment_coord[[i, 1]] = s.coord[[atom_add_pot[i], 1]];
                    s_fragment_coord[[i, 2]] = s.coord[[atom_add_pot[i], 2]];
                }

                let (mut pot_idwm, mut force_idwm_fragment) = idw_pot_force(&self.local_min, &s_fragment_coord, self.a_min, self.sigma_min);
                let mut pot_idwm_ts: f64;
                let mut force_idwm_ts: Array2<f64>;
                for i in 0..self.nearby_ts.len()
                {
                    (pot_idwm_ts, force_idwm_ts) = idw_pot_force(&self.nearby_ts[i], &s_fragment_coord, self.a_ts, self.sigma_ts[i]);
                    pot_idwm += pot_idwm_ts;
                    force_idwm_fragment += &force_idwm_ts;
                }

                let mut force_idwm: Array2<f64> = Array2::zeros(s.coord.raw_dim());
                for i in 0..atom_add_pot.len()
                {
                    force_idwm[[atom_add_pot[i], 0]] = force_idwm_fragment[[i, 0]];
                    force_idwm[[atom_add_pot[i], 1]] = force_idwm_fragment[[i, 1]];
                    force_idwm[[atom_add_pot[i], 2]] = force_idwm_fragment[[i, 2]];
                }

                return (pot_idwm, force_idwm)
            },
        }
    }
}










