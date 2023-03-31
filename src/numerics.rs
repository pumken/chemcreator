//! # Numerics
//!
//! The `numerics` module provides functions that give names to numbers.

use crate::naming::NamingError;

/// Returns the numeric prefix for group prefixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 20, a [`NamingError::GroupOccurrence`] without a group is returned.
pub fn minor_numeric(value: i32) -> Result<&'static str, NamingError> {
    let out = match value {
        1 => "",
        2 => "di",
        3 => "tri",
        4 => "tetra",
        5 => "penta",
        6 => "hexa",
        7 => "hepta",
        8 => "octa",
        9 => "nona",
        10 => "deca",
        11 => "undeca",
        12 => "dodeca",
        13 => "trideca",
        14 => "tetradeca",
        15 => "pentadeca",
        16 => "hexadeca",
        17 => "heptadeca",
        18 => "octadeca",
        19 => "nonadeca",
        20 => "icosa",
        21 => "henicosa",
        22 => "docosa",
        23 => "tricosa",
        24 => "tetracosa",
        25 => "pentacosa",
        26 => "hexacosa",
        27 => "heptacosa",
        28 => "octacosa",
        29 => "nonacosa",
        30 => "triconta",
        31 => "hentriconta",
        32 => "dotriconta",
        33 => "tritriconta",
        34 => "tetratriconta",
        35 => "pentatriconta",
        36 => "hexatriconta",
        37 => "heptatriconta",
        38 => "octatriconta",
        39 => "nonatriconta",
        40 => "tetraconta",
        41 => "hentetraconta",
        42 => "dotetraconta",
        _ => return Err(NamingError::GroupOccurrence(None, value)),
    };
    Ok(out)
}

/// Returns the numeric prefix for branch prefixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 20, a [`NamingError::GroupOccurrence`] without a group is returned.
pub fn branch_numeric(value: i32) -> Result<&'static str, NamingError> {
    let out = match value {
        1 => "",
        2 => "bis",
        3 => "tris",
        4 => "tetrakis",
        5 => "pentakis",
        6 => "hexakis",
        7 => "heptakis",
        8 => "octakis",
        9 => "nonakis",
        10 => "decakis",
        11 => "undecakis",
        12 => "dodecakis",
        13 => "tridecakis",
        14 => "tetradecakis",
        15 => "pentadecakis",
        16 => "hexadecakis",
        17 => "heptadecakis",
        18 => "octadecakis",
        19 => "nonadecakis",
        20 => "icosakis",
        _ => return Err(NamingError::GroupOccurrence(None, value)),
    };
    Ok(out)
}

/// Returns the numeric prefix for suffixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 84, a [`NamingError::CarbonCount`] is returned.
pub fn major_numeric(value: i32) -> Result<&'static str, NamingError> {
    let out = match value {
        1 => "meth",
        2 => "eth",
        3 => "prop",
        4 => "but",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        11 => "undec",
        12 => "dodec",
        13 => "tridec",
        14 => "tetradec",
        15 => "pentadec",
        16 => "hexadec",
        17 => "heptadec",
        18 => "octadec",
        19 => "nonadec",
        20 => "icos",
        21 => "henicos",
        22 => "docos",
        23 => "tricos",
        24 => "tetracos",
        25 => "pentacos",
        26 => "hexacos",
        27 => "heptacos",
        28 => "octacos",
        29 => "nonacos",
        30 => "triacont",
        31 => "untriacont",
        32 => "dotriacont",
        33 => "tritriacont",
        34 => "tetratriacont",
        35 => "pentatriacont",
        36 => "hexatriacont",
        37 => "heptatriacont",
        38 => "octatriacont",
        39 => "nonatriacont",
        40 => "tetracont",
        41 => "hentetracont",
        42 => "dotetracont",
        43 => "tritetracont",
        44 => "tetratetracont",
        45 => "pentatetracont",
        46 => "hexatetracont",
        47 => "heptatetracont",
        48 => "octatetracont",
        49 => "nonatetracont",
        50 => "pentacont",
        51 => "henpentacont",
        52 => "dopentacont",
        53 => "tripentacont",
        54 => "tetrapentacont",
        55 => "pentapentacont",
        56 => "hexapentacont",
        57 => "heptapentacont",
        58 => "octapentacont",
        59 => "nonapentacont",
        60 => "hexacont",
        61 => "henhexacont",
        62 => "dohexacont",
        63 => "trihexacont",
        64 => "tetrahexacont",
        65 => "pentahexacont",
        66 => "hexahexacont",
        67 => "heptahexacont",
        68 => "octahexacont",
        69 => "nonahexacont",
        70 => "heptacont",
        71 => "henheptacont",
        72 => "doheptacont",
        73 => "triheptacont",
        74 => "tetraheptacont",
        75 => "pentaheptacont",
        76 => "hexaheptacont",
        77 => "heptaheptacont",
        78 => "octaheptacont",
        79 => "nonaheptacont",
        80 => "octacont",
        81 => "henoctacont",
        82 => "dooctacont",
        83 => "trioctacont",
        84 => "tetraoctacont",
        _ => return Err(NamingError::CarbonCount(value)),
    };
    Ok(out)
}
