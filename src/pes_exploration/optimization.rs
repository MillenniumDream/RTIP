//! Containing the functions for local optimization

use crate::common::error::*;
use crate::common::constants::{GOLDEN_RATIO1, GOLDEN_RATIO2};
use crate::pes_exploration::traits::PES;
use crate::pes_exploration::system::System;
use ndarray::Array2;





/// One-dimensional optimization using the golden section method
/// Input a one-dimensional function fun, which is assumed to decrease along +x direction near x = 0.0
/// the functional value f0 at x = 0.0, and the trial step delta_x > 0.0
/// find out the minimum value of the function within (0, +oo).
/// First, determine the searching region by the extrapolation method with the scaling factor = 1.618
/// And then use the golden section method to find out the minimun functional value
///
/// # Parameters
/// ```
/// fun: the input one-dimensional function
/// f0: the functional value at x = 0.0
/// delta_x: the trial step
/// epsilon: the convergence criteria for the functional value
/// ```
pub fn min_1d<T> (fun: &mut T, f0: f64, epsilon: f64) -> (f64, f64)
where T: FnMut(f64) -> f64
{
    // Determine the searching region by the extrapolation method
    let mut delta_x: f64 = 0.01;
    let mut x0: f64 = 0.0;
    let mut x1: f64 = delta_x;
    let mut f1: f64 = fun(x1);

    // If the trial step is too large, decrease the step unless it is smaller than 0.0000000001
    while (f1 > f0) && (delta_x > 0.000000000000001)
    {
        delta_x /= 10.0;
        x1 = delta_x;
        f1 = fun(x1);
    }

    if f1 > f0
    {
        panic!("{}", error_min_1d());       // Assert that the input fun along +x direction is decreasing
    }

    delta_x *= GOLDEN_RATIO2;
    let mut x2: f64 = x1 + delta_x;
    let mut f2: f64 = fun(x2);

    while (f2 < f1) && (x1 < 0.1)
    {
        delta_x *= GOLDEN_RATIO2;
        x0 = x1;
        x1 = x2;
        x2 = x1 + delta_x;
        f1 = f2;
        f2 = fun(x2);
    }



    // Control the step not to be too large
    if (f2 < f1) && (x1 >= 0.1)
    {
        return (x1, f1)
    }



    // Find out the minimum functional value in the region (a, b) by the golden section method
    let mut a: f64 = x0;
    let mut b: f64 = x2;
    let mut u: f64 = x1;
    let mut v: f64 = a + GOLDEN_RATIO1 * (b - a);
    let mut fu: f64 = f1;
    let mut fv: f64 = fun(v);
    
    loop
    {
        if fu > fv
        {
            if (fu - fv) < epsilon
            {
                return (v, fv)
            }
            else
            {
                a = u;
                u = v;
                v = a + GOLDEN_RATIO1 * (b - a);
                fu = fv;
                fv = fun(v);
            }
        }
        else
        {
            if (fv - fu) < epsilon
            {
                return (u, fu)
            }
            else
            {
                b = v;
                v = u;
                u = b - GOLDEN_RATIO1 * (b - a);
                fv = fu;
                fu = fun(u);
            }
        }
    }
}





/// Assume pot_total = pot_real + pot_bias
/// search the minimum potential energy structure along the force direction
///
/// # Parameters
/// ```
/// real_pes: the input real potential energy surface (PES)
/// bias_pes: the input bias PES
/// s: the input current step structure
/// force_total: the input total force (force_real + force_bias), along which to search the minimum potential energy
/// pot_total: the input total potential energy (pot_real + pot_bias) of the current step structure s
/// epsilon: the energy criterion for convergence (Unit: Hartree)
/// delta_coord: the output coordinate displacement, i.e., (s.coord + delta_coord) has the minimum pot_total along the force_total direction
/// ```
pub fn min_1d_real_bias<P1: PES, P2: PES>(real_pes: &P1, bias_pes: &P2, s: &System, force_total: &Array2<f64>, pot_total: f64, epsilon: f64) -> Array2<f64>
{
    let mut s1: System = s.clone();
    let f_total: f64 = (force_total * force_total).sum().sqrt();        // the norm of the total force

    // Difine a function, input the displacement along the force direction, output the total potential energy
    let mut fun = |x: f64| -> f64
    {
        // Obtain a structure along the force direction
        s1.coord = &s.coord + force_total * (x/f_total);

        // Calculate the real potential energy and bias potential energy, and return their summation
        real_pes.get_energy(&s1) + bias_pes.get_energy(&s1)
    };

    let (x_min, _pot_total) = min_1d(&mut fun, pot_total, epsilon);
    let delta_coord: Array2<f64> = force_total * (x_min/f_total);
    delta_coord
}










