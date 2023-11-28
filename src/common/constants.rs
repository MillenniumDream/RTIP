//! Contains mathematical, physical, and chemical constants.
use crate::common::error::*;
use phf::phf_map;
use mpi::topology::Rank;





pub const ROOT_RANK: Rank = 0;





// Mathematical
pub const PI: f64 = 3.141592653589793;
pub const GOLDEN_RATIO1: f64 = 0.618033988749895;
pub const GOLDEN_RATIO2: f64 = 1.618033988749895;





// Physical

// Unit Conversion
pub const C_LIGHT: f64 = 299792458.0;
pub const H_PLANCK: f64 = 6.62606896E-34;
pub const BOLTZMANN: f64 = 1.3806504E-23;
pub const RYBDERG: f64 = 10973731.568527;

pub const BOHR_TO_ANGSTROM: f64 = 0.52917720859;
pub const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_TO_ANGSTROM;

pub const HARTREE_TO_JOULE: f64 = 2.0 * RYBDERG * H_PLANCK * C_LIGHT;
pub const JOULE_TO_HARTREE: f64 = 1.0 / HARTREE_TO_JOULE;

pub const AU_TO_FEMTOSECOND: f64 = 1.0E15 / (4.0 * PI * RYBDERG * C_LIGHT);
pub const FEMTOSECOND_TO_AU: f64 = 1.0 / AU_TO_FEMTOSECOND;





// Chemical

#[derive(Clone)]
#[derive(Debug)]
pub enum Element
{
    H, He,
    Li, Be, B, C, N, O, F, Ne,
    Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,
    Cs, Ba, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg,
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
    Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
}

// 'STR_ELEMENT' is a static structure of type 'phf::Map', initialized by macro 'phf_map'
static STR_ELEMENT: phf::Map<&'static str, Element> = phf_map!
{
    "H" => Element::H,
    "He" => Element::He,

    "Li" => Element::Li,
    "Be" => Element::Be,
    "B" => Element::B,
    "C" => Element::C,
    "N" => Element::N,
    "O" => Element::O,
    "F" => Element::F,
    "Ne" => Element::Ne,

    "Na" => Element::Na,
    "Mg" => Element::Mg,
    "Al" => Element::Al,
    "Si" => Element::Si,
    "P" => Element::P,
    "S" => Element::S,
    "Cl" => Element::Cl,
    "Ar" => Element::Ar,

    "K" => Element::K,
    "Ca" => Element::Ca,
    "Sc" => Element::Sc,
    "Ti" => Element::Ti,
    "V" => Element::V,
    "Cr" => Element::Cr,
    "Mn" => Element::Mn,
    "Fe" => Element::Fe,
    "Co" => Element::Co,
    "Ni" => Element::Ni,
    "Cu" => Element::Cu,
    "Zn" => Element::Zn,
    "Ga" => Element::Ga,
    "Ge" => Element::Ge,
    "As" => Element::As,
    "Se" => Element::Se,
    "Br" => Element::Br,
    "Kr" => Element::Kr,

    "Rb" => Element::Rb,
    "Sr" => Element::Sr,
    "Y" => Element::Y,
    "Zr" => Element::Zr,
    "Nb" => Element::Nb,
    "Mo" => Element::Mo,
    "Tc" => Element::Tc,
    "Ru" => Element::Ru,
    "Rh" => Element::Rh,
    "Pd" => Element::Pd,
    "Ag" => Element::Ag,
    "Cd" => Element::Cd,
    "In" => Element::In,
    "Sn" => Element::Sn,
    "Sb" => Element::Sb,
    "Te" => Element::Te,
    "I" => Element::I,
    "Xe" => Element::Xe,

    "Cs" => Element::Cs,
    "Ba" => Element::Ba,
    "Hf" => Element::Hf,
    "Ta" => Element::Ta,
    "W" => Element::W,
    "Re" => Element::Re,
    "Os" => Element::Os,
    "Ir" => Element::Ir,
    "Pt" => Element::Pt,
    "Au" => Element::Au,
    "Hg" => Element::Hg,
    "Tl" => Element::Tl,
    "Pb" => Element::Pb,
    "Bi" => Element::Bi,
    "Po" => Element::Po,
    "At" => Element::At,
    "Rn" => Element::Rn,

    "Fr" => Element::Fr,
    "Ra" => Element::Ra,
    "Rf" => Element::Rf,
    "Db" => Element::Db,
    "Sg" => Element::Sg,
    "Bh" => Element::Bh,
    "Hs" => Element::Hs,
    "Mt" => Element::Mt,
    "Ds" => Element::Ds,
    "Rg" => Element::Rg,

    "La" => Element::La,
    "Ce" => Element::Ce,
    "Pr" => Element::Pr,
    "Nd" => Element::Nd,
    "Pm" => Element::Pm,
    "Sm" => Element::Sm,
    "Eu" => Element::Eu,
    "Gd" => Element::Gd,
    "Tb" => Element::Tb,
    "Dy" => Element::Dy,
    "Ho" => Element::Ho,
    "Er" => Element::Er,
    "Tm" => Element::Tm,
    "Yb" => Element::Yb,
    "Lu" => Element::Lu,

    "Ac" => Element::Ac,
    "Th" => Element::Th,
    "Pa" => Element::Pa,
    "U" => Element::U,
    "Np" => Element::Np,
    "Pu" => Element::Pu,
    "Am" => Element::Am,
    "Cm" => Element::Cm,
    "Bk" => Element::Bk,
    "Cf" => Element::Cf,
    "Es" => Element::Es,
    "Fm" => Element::Fm,
    "Md" => Element::Md,
    "No" => Element::No,
    "Lr" => Element::Lr,
};

impl Element
{
    pub fn from_str(element: &str) -> Self
    {
        STR_ELEMENT.get(element).cloned().expect(&error_element())
    }

    pub fn get_mass(&self) -> f64
    {
        match self
        {
            Element::H => 1837.362218829611,
            Element::C => 21894.16671795623,
            Element::O => 29165.12201514224,
            Element::N => 25532.65213254827,
            Element::S => 58450.91924794280,
            Element::P => 56461.71406415092,
            _ => panic!("No such element"),
        }
    }
}










