//! About the output files.
use std::fs;
use crate::common::constants::ROOT_RANK;
use crate::common::error::*;
use mpi::traits::Communicator;





/// Specify the output files for roto-translational invariant potential (RTIP)
///
/// # Parameters
/// ```
/// comm: the input communicator
/// index: the input index that specifies where to output the files
/// str_output_file: the structure output file name
/// rtip_output_file: the RTIP output file name
/// ```
pub fn output_rtip<C: Communicator>(comm: &C, index: Option<i32>) -> (String, String)
{
    match index
    {
        // If index exists, create a directory and output the files to it
        Some(index) =>
        {
            if comm.rank() == ROOT_RANK
            {
                let dir: String = format!("{}", index);
                fs::create_dir(&dir).expect(&error_dir("creating", &dir));
            }

            let str_output_file: String = format!("{}/rtip.pdb", index);
            let rtip_output_file: String = format!("{}/rtip.out", index);

            (str_output_file, rtip_output_file)
        },

        // If index non-exists, output the files to the current directory
        None =>
        {
            let str_output_file: String = format!("rtip.pdb");
            let rtip_output_file: String = format!("rtip.out");

            (str_output_file, rtip_output_file)
        },
    }
}

/// Specify the output file for CP2K
///
/// # Parameters
/// ```
/// index: the input index that specifies where to output the CP2K file
/// cp2k_output_file: the CP2K output file name
/// ```
pub fn output_cp2k(index: Option<i32>) -> String
{
    match index
    {
        // If index exists, create a directory and output the file to it
        Some(index) =>
        {
            format!("{}/cp2k.out", index)
        },

        // If index non-exists, output the file to the current directory
        None =>
        {
            format!("cp2k.out")
        },
    }
}










