use crate::pes_exploration::system::System;
use crate::io::input::Para;





/// The structure containing all the information about the Gaussian repulsive potential,
/// which is introduced to push the system out of the local minimum along an unknown pathway.
///
/// # Fields
/// ```
/// local_min: the local minimum structure
/// nearby_ts: the known nearby transition state (TS) structures
/// para: containing all the parameters for pushing the structure out of the local minimum
/// str_output_file: the structure output file name
/// output_file: the RTIP output file name
/// ```
pub struct RepulsivePot<'a>
{
    pub local_min: System,
    pub nearby_ts: Vec<System>,
    pub para: &'a Para,
    pub str_output_file: String,
    pub output_file: String,
}





/// The structure containing all the information about the Gaussian attractive potential,
/// which is introduced to draft the system towards the final state.
///
/// # Fields
/// ```
/// initial_state: the current system to be drafted
/// final_state: the objective system to be drafted towards
/// para: containing all the parameters for drafting the initial state towards the final state
/// str_output_file: the structure output file name
/// output_file: the RTIP output file name
/// ```
pub struct AttractivePot<'a>
{
    pub initial_state: System,
    pub final_state: System,
    pub para: &'a Para,
    pub str_output_file: String,
    pub output_file: String,
}





/// The structure containing all the information about the synthesis potential,
/// which is introduced to synthesize a product from several molecules.
///
/// # Fields
/// ```
/// initial_state: the initial system containing several separated molecules
/// mol_index: atomic index of the molecules for synthesis
/// para: containing all the parameters for the molecule synthesis
/// str_output_file: the structure output file name
/// output_file the RTIP output file name
/// ```
pub struct SynthesisPot<'a>
{
    pub initial_state: System,
    pub mol_index: Vec<Vec<usize>>,
    pub para: &'a Para,
    pub str_output_file: String,
    pub output_file: String,
}










