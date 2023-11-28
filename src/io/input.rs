//! About the input files.





/// The structure to contain the parameters of the roto-translational invariant potential (RTIP).
///
/// # Fields
/// ```
/// a0: specifying how quickly the Gaussian potential is increasing (Unit: Hartree)
/// max_step: specifying the maximum iterative steps for the Gaussian potential increasing
/// pot_drop: when the real potential energy (i.e. e_cp2k) becomes smaller than (e_max - e_drop), remove the Gaussian potential and perform a local optimization (Unit: Hartree)
/// pot_epsilon: in each steeping searching, when the difference of two adjacent e_total (i.e. e_cp2k+e_rtip) becomes smaller than e_epsilon, stop the searching (Unit: Hartree/atom)
/// f_epsilon: in local optimization without Gaussian potential, when force_rtip becomes smaller than force_epsilon, stop the optimization (Unit: Hartree/(Bohr*atom))
/// ```
#[derive(Clone)]
pub struct Para
{
    // General parameters
    pub a0: f64,
    pub scale_ts_a0: f64,
    pub scale_ts_sigma: Option<f64>,
    pub max_step: usize,
    pub print_step: usize,

    // Pathway sampling parameters
    pub pot_climb: f64,
    pub pot_drop: f64,
    pub pot_epsilon: f64,
    pub f_epsilon: f64,

    // MD parameters
    pub dt: f64,
    pub tau: f64,
    pub temp_bath: f64,
}





impl Para
{
    pub fn new() -> Self
    {
        Para
        {
            a0: 0.002,
            scale_ts_a0: 1.0,
            scale_ts_sigma: Some(0.25),
            max_step: 1500,
            print_step: 1,

            pot_climb: 0.185,
            pot_drop: 0.02,
            pot_epsilon: 0.00005,
            f_epsilon: 0.001,

            dt: 0.5,
            tau: 500.0,
            temp_bath: 1000.0,
        }
    }
}










