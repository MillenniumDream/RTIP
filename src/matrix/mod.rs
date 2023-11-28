//! About the matrix operations.

use crate::common::constants::PI;
use ndarray::{Array1, Array2, array};
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;





pub fn rand_rot() -> Array2<f64>
{
    let alpha: f64 = Array1::random(1, Uniform::new(0.0, 2.0 * PI))[0];
    let beta: f64 = Array1::random(1, Uniform::new(0.0, 2.0 * PI))[0];
    let gamma: f64 = Array1::random(1, Uniform::new(0.0, 2.0 * PI))[0];

    return array![[alpha.cos()*gamma.cos() - beta.cos()*alpha.sin()*gamma.sin(), -beta.cos()*gamma.cos()*alpha.sin() - alpha.cos()*gamma.sin(), alpha.sin()*beta.sin()],
                  [gamma.cos()*alpha.sin() + alpha.cos()*beta.cos()*gamma.sin(), alpha.cos()*beta.cos()*gamma.cos() - alpha.sin()*gamma.sin(), -alpha.cos()*beta.sin()],
                  [beta.sin()*gamma.sin(), gamma.cos()*beta.sin(), beta.cos()]]
}










