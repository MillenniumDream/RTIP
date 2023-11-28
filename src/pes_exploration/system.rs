//! About the structure, constraint, and status of the research system.





use crate::common::constants::{Element, BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR};
use crate::external::cp2k::*;
use crate::common::error::*;
use std::fs;
use std::fs::File;
use std::io::Write;
use ndarray::{Array, Array2};





/// The basic structure to describe the research system (e.g. molecule and crystal).
///
/// # Fields
/// ```
/// natom: the number of atoms in the research system
/// coord: the coordinates of the atoms in the research system (natom*3 Array, Unit: Bohr)
/// cell: the unit cell of the research system (3*3 Array, Unit: Bohr)
/// atom_type: the types of the atoms in the research system (natom Vec)
/// add_pot: whether to add the biase Gaussian potential on the atoms or not (natom Vec)
/// mutable: whether the atomic positions is mutable or not (natom*3 Array)
/// pot: the potential energy of the system
/// ```
#[derive(Clone)]
pub struct System
{
    pub natom: usize,
    pub coord: Array2<f64>,
    pub cell: Option< Array2<f64> >,
    pub atom_type: Option< Vec<Element> >,
    pub atom_add_pot: Option< Vec<usize> >,
    pub mutable: Option< Array2<bool> >,
    pub pot: f64,
}





impl System
{
    /// Input a CP2K force environment, get the number and coordinates of the particles from CP2K, and construct a System to return
    ///
    /// # Parameters
    /// ```
    /// force_env: the previously built CP2K force environment
    /// ```
    ///
    /// # Examples
    /// ```
    /// let s = System::from_cp2k(force_env);
    /// ```
    pub fn from_cp2k(force_env: ForceEnv) -> Self
    {
        let mut nparticle: i32 = 0;
        cp2k_get_nparticle(force_env, &mut nparticle);
        let natom: usize = nparticle.try_into().expect(&error_type_transformation("nparticle", "i32", "usize"));           // From i32 to usize

        let mut pos: Vec<f64> = vec![0.0; natom*3];
        cp2k_get_positions(force_env, &mut pos, nparticle*3);
        let coord: Array2<f64> = Array::from_shape_vec((natom,3), pos).expect(&error_none_value("coord"));          // From Vec to Array2
        
        let mut pot: f64 = 0.0;
        cp2k_get_potential_energy(force_env, &mut pot);

        System
        {
            natom,
            coord,
            cell: None,
            atom_type: None,
            atom_add_pot: None,
            mutable: None,
            pot,
        }
    }

    /// Input a CP2K force environment, consume the System and set the coordinates of the particles for CP2K
    ///
    /// # Parameters
    /// ```
    /// force_env: the previously built CP2K force environment
    /// ```
    ///
    /// # Examples
    /// ```
    /// s.to_cp2k(force_env)
    /// ```
    pub fn to_cp2k(&self, force_env: ForceEnv)
    {
        let n_el: i32 = (self.natom * 3).try_into().expect(&error_type_transformation("natom", "usize", "i32"));
        let pos: &[f64] = self.coord.as_slice().expect(&error_none_value("pos"));
        cp2k_set_positions(force_env, pos, n_el);
    }



    /// Create a new PDB file (if already existed, truncate it), or open an old PDB file, and write the structure (in Angstrom) to it
    ///
    /// # Parameters
    /// ```
    /// filename: name of the PDB file to be writen
    /// create_new_file: whether to create a new PDB file or not
    /// step: current step of the structure
    /// ```
    ///
    /// # Examples
    /// ```
    /// s.write_pdb("filename.pdb", true, 1);
    /// ```
    pub fn write_pdb(&self, filename: &str, create_new_file: bool, step: usize)
    {
        let mut pdb = match create_new_file
        {
            // If create_new_file == true, create a new file and write the PDB TITLE
            true =>
            {
                let mut pdb = File::create(filename).expect(&error_file("creating", filename));       // Create the PDB file
                pdb.write_all(b"TITLE     PDB file created by RTIP\n").expect(&error_file("writing", filename));       // Write the PDB TITLE
                pdb
            },
            // If create_new_file == false, open an old file and append to it
            false =>
            {
                File::options().append(true).open(filename).expect(&error_file("opening", filename))
            },
        };

        // Write the potential energy
        pdb.write_all(format!("REMARK    , Step = {:8}, E = {:15.8}\n", step, self.pot).as_bytes()).expect(&error_file("writing", filename));
        
        // Write the atomic elements and coordinates
        match &self.atom_type
        {
            // If the elements of the atoms are reserved in System
            Some(atom_type) =>
            {
                for i in 0..self.natom
                {
                    // Warning: The atomic coordinates should be transformd from Bohr to Angstrom
                    pdb.write_all(format!("ATOM  {:>5} {:<4}              {:8.3}{:8.3}{:8.3}  0.00  0.00          {:>2}\n", i+1, format!("{:?}", atom_type[i]), self.coord[[i,0]] * BOHR_TO_ANGSTROM, self.coord[[i,1]] * BOHR_TO_ANGSTROM, self.coord[[i,2]] * BOHR_TO_ANGSTROM, format!("{:?}", atom_type[i]) ).as_bytes()).expect(&error_file("writing", filename));
                }
            },

            // If the elements of the atoms are not reserved in System
            None =>
            {
                for i in 0..self.natom
                {
                    // Warning: The atomic coordinates should be transformed from Bohr to Angstrom
                    pdb.write_all(format!("ATOM  {:>5} {:<4}              {:8.3}{:8.3}{:8.3}  0.00  0.00          {:>2}\n", i+1, "XX", self.coord[[i,0]] * BOHR_TO_ANGSTROM, self.coord[[i,1]] * BOHR_TO_ANGSTROM, self.coord[[i,2]] * BOHR_TO_ANGSTROM, "XX" ).as_bytes()).expect(&error_file("writing", filename));
                }
            },
        }

        pdb.write_all(b"END\n").expect(&error_file("writing", filename));        // Write the PDB END
    }



    /// Read a system from a XYZ file
    ///
    /// # Parameters
    /// ```
    /// filename: name of the XYZ file to read from
    /// ```
    ///
    /// # Examples
    /// ```
    /// let s = System::read_xyz("filename.xyz");
    /// ```
    pub fn read_xyz(filename: &str) -> Self
    {
        let content = fs::read_to_string(filename).expect(&error_file("reading", filename));                // Read the whole file
        let mut line = content.lines();                 // Take the iteractor for the lines of the file
        
        // Read the number of atom from the first line
        let natom: usize = match line.next()
        {
            Some(value) => value.trim().parse().expect(&error_read_xyz(filename)),
            None => panic!("{}", error_read_xyz(filename)),
        };
        line.next();                // Ignore the second line
        
        // Read the element type and Cartesian coordinates of each atom iteractively
        let mut atom_type: Vec<Element> = vec![Element::H; natom];
        let mut coord: Array2<f64> = Array2::zeros((natom, 3));
        for i in 0..natom
        {
            // Read the element type and atomic coordinates
            let type_coord: Vec<&str> = match line.next()
            {
                Some(value) => value.split_whitespace().collect(),
                None => panic!("{}", error_read_xyz(filename)),
            };
            if type_coord.len() < 4
            {
                panic!("{}", error_read_xyz(filename));
            }

            // Restore the element type and atomic coordinates
            atom_type[i] = Element::from_str(type_coord[0]);
            coord[[i, 0]] = type_coord[1].parse().expect(&error_read_xyz(filename));
            coord[[i, 1]] = type_coord[2].parse().expect(&error_read_xyz(filename));
            coord[[i, 2]] = type_coord[3].parse().expect(&error_read_xyz(filename));
        }
        // Warning: The atomic coordinates should be transformed from Angstrom to Bohr
        coord *= ANGSTROM_TO_BOHR;

        System
        {
            natom,
            coord,
            cell: None,
            atom_type: Some(atom_type),
            atom_add_pot: None,
            mutable: None,
            pot: 0.0,
        }
    }

    /// Create a new XYZ file (if already existed, truncate it), or open an old XYZ file, and write the structure (in Angstrom) to it
    ///
    /// # Parameters
    /// ```
    /// filename: name of the XYZ file to be writen
    /// create_new_file: whether to create a new XYZ file or not
    /// step: current step of the structure
    /// ```
    ///
    /// # Examples
    /// ```
    /// s.write_xyz("filename.xyz", true, 1);
    /// ```
    pub fn write_xyz(&self, filename: &str, create_new_file: bool, step: usize)
    {
        let mut xyz = match create_new_file
        {
            // If create_new_file == true, create a new file
            true =>
            {
                File::create(filename).expect(&error_file("creating", filename))
            },
            // If create_new_file == false, open an old file and append to it
            false =>
            {
                File::options().append(true).open(filename).expect(&error_file("opening", filename))
            },
        };

        // Write the number of atoms
        xyz.write_all(format!("{:8}\n", self.natom).as_bytes()).expect(&error_file("writing", filename));

        // Write the potential energy
        xyz.write_all(format!("Step = {:8}, E = {:15.8}\n", step, self.pot).as_bytes()).expect(&error_file("writing", filename));

        // Write the atomic elements and coordinates
        match &self.atom_type
        {
            // If the elements of the atoms are reserved in System
            Some(atom_type) =>
            {
                for i in 0..self.natom
                {
                    // Warning: The atomic coordinates should be transformd from Bohr to Angstrom
                    xyz.write_all(format!("{:>4} {:20.10} {:20.10} {:20.10}\n", format!("{:?}", atom_type[i]), self.coord[[i,0]] * BOHR_TO_ANGSTROM, self.coord[[i,1]] * BOHR_TO_ANGSTROM, self.coord[[i,2]] * BOHR_TO_ANGSTROM).as_bytes()).expect(&error_file("writing", filename));
                }
            },

            // If the elements of the atoms are not reserved in System
            None =>
            {
                for i in 0..self.natom
                {
                    // Warning: The atomic coordinates should be transformed from Bohr to Angstrom
                    xyz.write_all(format!("{:>4} {:20.10} {:20.10} {:20.10}\n", "XX", self.coord[[i,0]] * BOHR_TO_ANGSTROM, self.coord[[i,1]] * BOHR_TO_ANGSTROM, self.coord[[i,2]] * BOHR_TO_ANGSTROM).as_bytes()).expect(&error_file("writing", filename));
                }
            },
        }
    }
}










