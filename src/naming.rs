//! # Naming
//!
//! The `naming` module contains functions that makes final validations and finds the name of
//! an organic molecule as a [`Branch`].

use crate::chain::primary_chain;
use crate::groups::Fallible;
use crate::groups::InvalidGraphError::Other;
use crate::molecule::Group::Alkane;
use crate::molecule::{Atom, Branch, Cell, Group, Substituent};
use crate::spatial::GridState;
use crate::{chain, groups, validation};
use chain::get_all_chains;
use groups::link_groups;
use std::fmt::{Display, Formatter};
use thiserror::Error;
use validation::{check_structure, check_valence};

/// Determines the name of the molecule on the given `graph`.
///
/// ## Errors
///
/// If the molecule on the given `graph` is discontinuous, cyclic, or contains invalid bonding,
/// an [`InvalidGraphError`] will be returned.
pub fn name_molecule(graph: &GridState) -> Fallible<String> {
    let cells = graph
        .find_all(|cell| cell.is_atom())
        .iter()
        .map(|&cell| {
            if let Cell::Atom(it) = cell {
                it.to_owned()
            } else {
                panic!("is_atom check failed in name_molecule")
            }
        })
        .collect::<Vec<Atom>>();

    // Initial checks
    if graph.is_empty() {
        return Ok("".to_string());
    }
    check_structure(graph)?;
    check_valence(cells, graph)?;

    // Preliminary chain
    let all_chains = get_all_chains(graph)?;
    let chain = primary_chain(graph, all_chains, None)?;

    // Group-linked branch
    let mut branch = link_groups(graph, chain, None)?;
    branch = branch.index_corrected()?;

    process_name(branch).map_err(|e| Other(e.to_string()))
}

/// Generates the name of the organic molecule by the given `branch`.
pub(crate) fn process_name(branch: Branch) -> Result<String, NamingError> {
    let len = branch.chain.len() as i32;
    let collection = SubFragmentCollection::new(branch);

    let name = format!(
        "{}{}{}{}",
        prefix(collection.secondary_group_fragments())?,
        major_numeric(len)?,
        bonding(collection.unsaturated_group_fragments())?,
        suffix(collection.primary_group_fragment())?,
    );

    if let Some(sub) = substitute_names(&name) {
        Ok(sub.to_string())
    } else {
        Ok(name)
    }
}

/// Returns a corrected name for rules that are not worth it to add in.
pub(crate) fn substitute_names(name: &str) -> Option<&str> {
    let out = match name {
        "1-fluoromethanoyl fluoride" => "methanoyl difluoride",
        "1-chloromethanoyl chloride" => "methanoyl dichloride",
        "1-bromomethanoyl bromide" => "methanoyl dibromide",
        "1-iodomethanoyl iodide" => "methanoyl diiodide",
        _ => return None,
    };
    Some(out)
}

/// Returns the prefix sequence string of the given `fragments`.
///
/// ## Errors
///
/// If the occurrences of one of the `fragments` exceeds the limit of [`minor_numeric`], a
/// [`NamingError::GroupOccurrence`] is returned.
pub(crate) fn prefix(mut fragments: Vec<SubFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.subst.to_string());
    let mut out = vec![];

    for fragment in fragments.into_iter() {
        let locants = if matches!(&fragment.subst, Substituent::Branch(branch) if has_number(&branch.to_string()))
        {
            NamingError::add_substituent(complex_branch_locants(fragment.locants), &fragment.subst)?
        } else {
            NamingError::add_substituent(locants(fragment.locants), &fragment.subst)?
        };
        out.push(format!("{}{}", locants, fragment.subst))
    }

    Ok(out.join("-"))
}

fn has_number(s: &str) -> bool {
    for c in s.chars() {
        if c.is_ascii_digit() {
            return true;
        }
    }
    false
}

fn bonding(mut fragments: Vec<SubFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.subst.to_string());
    if fragments.is_empty() {
        Ok("an".to_string())
    } else {
        let mut sequence = vec![];

        for fragment in fragments.into_iter() {
            let locants = NamingError::add_substituent(locants(fragment.locants), &fragment.subst)?;

            sequence.push(format!("{}{}", locants, fragment.subst));
        }

        let mut out = sequence.join("-");
        out.insert(0, '-');
        Ok(out)
    }
}

fn suffix(fragment: SubFragment) -> Result<String, NamingError> {
    if let Substituent::Group(group) = fragment.subst {
        let locants =
            NamingError::add_substituent(locants(fragment.locants.clone()), &fragment.subst)?;

        let suffix = match group {
            Group::Carboxyl if fragment.locants.len() == 2 => return Ok("edioic acid".to_string()),
            Group::Carboxyl => return Ok("oic acid".to_string()),
            Group::Aldehyde if fragment.locants.len() == 2 => return Ok("edial".to_string()),
            Group::Aldehyde => return Ok("al".to_string()),
            Group::Carbonyl => "one",
            Group::Hydroxyl => "ol",
            Group::Nitrile if fragment.locants.len() == 2 => return Ok("edinitrile".to_string()),
            Group::Nitrile => return Ok("onitrile".to_string()),
            Group::AcidHalide(it) if fragment.locants.len() == 2 => {
                return Ok(format!("edioyl di{it}"));
            }
            Group::AcidHalide(it) => return Ok(format!("oyl {it}")),
            Group::Amide if fragment.locants.len() == 2 => return Ok("ediamide".to_string()),
            Group::Amide => return Ok("amide".to_string()),
            Group::Amine => "amine",
            _ => return Ok("e".to_string()),
        };

        Ok(format!("-{locants}{suffix}"))
    } else {
        panic!("branch provided to suffix function")
    }
}

/// Returns a [`String`] containing a locant prefix with a minor numeric prefix appended.
///
/// ## Errors
///
/// If the number of `locations` given exceeds the limit of [`minor_numeric`], a
/// [`NamingError::GroupOccurrence`] without a group is returned.
pub fn locants(mut locations: Vec<i32>) -> Result<String, NamingError> {
    locations.sort();

    let numbers = locations
        .iter()
        .map(|&it| (it + 1).to_string())
        .collect::<Vec<String>>()
        .join(",");

    let out = format!("{numbers}-{}", minor_numeric(locations.len() as i32)?,);
    Ok(out)
}

pub fn complex_branch_locants(mut locations: Vec<i32>) -> Result<String, NamingError> {
    locations.sort();

    let numbers = locations
        .iter()
        .map(|&it| (it + 1).to_string())
        .collect::<Vec<String>>()
        .join(",");

    let out = format!("{numbers}-{}", branch_numeric(locations.len() as i32)?,);
    Ok(out)
}

/// Returns the numeric prefix for group prefixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 10, a [`NamingError::GroupOccurrence`] without a group is returned.
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
        _ => return Err(NamingError::GroupOccurrence(None, value)),
    };
    Ok(out)
}

/// Returns the numeric prefix for branch prefixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 10, a [`NamingError::GroupOccurrence`] without a group is returned.
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
        _ => return Err(NamingError::GroupOccurrence(None, value)),
    };
    Ok(out)
}

/// Returns the numeric prefix for suffixes for the given `value`.
///
/// ## Errors
///
/// If `value` exceeds 32, a [`NamingError::CarbonCount`] is returned.
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
        32 => "duotriacont",
        _ => return Err(NamingError::CarbonCount(value)),
    };
    Ok(out)
}

/// A collection of [`SubFragment`]s representing the [`Substituents`] attached to a [`Branch`]
/// and the locations of their occurrences.
#[derive(Clone, Debug, PartialEq)]
pub struct SubFragmentCollection {
    pub collection: Vec<SubFragment>,
}

impl SubFragmentCollection {
    pub fn new(branch: Branch) -> SubFragmentCollection {
        let mut out = SubFragmentCollection { collection: vec![] };

        branch
            .groups
            .into_iter()
            .enumerate()
            .for_each(|(index, link)| {
                link.into_iter().for_each(|it| {
                    out.push(it, index as i32);
                })
            });

        out
    }

    /// Adds an occurrence of the given `substituent` to the collection. If this is not the first
    /// occurrence of the `substituent`, the location is simply added to its corresponding
    /// [`SubFragment`].
    fn push(&mut self, substituent: Substituent, index: i32) {
        let item = self
            .collection
            .iter_mut()
            .find(|it| it.subst == substituent);

        match item {
            Some(it) => it.locants.push(index),
            None => self
                .collection
                .push(SubFragment::new(vec![index], substituent)),
        }
    }

    /// Returns the primary [`Group`] of the collection, i.e., that which should be used as a suffix.
    /// If no primary groups exist, [`Alkane`] is returned.
    pub fn primary_group(&self) -> Group {
        let primary = self
            .collection
            .iter()
            .map(|fragment| fragment.subst.to_owned())
            .filter(|group| group.seniority().is_some())
            .max_by_key(|group| group.seniority());

        match primary {
            Some(Substituent::Group(it)) => it,
            _ => Alkane,
        }
    }

    /// Returns the [`SubFragment`] for the primary group.
    pub fn primary_group_fragment(&self) -> SubFragment {
        self.collection
            .iter()
            .find(|&fragment| fragment.subst == Substituent::Group(self.primary_group()))
            .map_or_else(SubFragment::default, SubFragment::to_owned)
    }

    /// Returns all non-primary groups not including
    pub fn secondary_group_fragments(&self) -> Vec<SubFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.subst != Substituent::Group(self.primary_group()))
            .filter(|&fragment| !fragment.subst.is_cxc_group())
            .map(SubFragment::to_owned)
            .collect::<Vec<SubFragment>>()
    }

    pub fn unsaturated_group_fragments(&self) -> Vec<SubFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.subst.is_cxc_group())
            .map(SubFragment::to_owned)
            .collect::<Vec<SubFragment>>()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SubFragment {
    pub locants: Vec<i32>,
    pub subst: Substituent,
}

impl SubFragment {
    pub fn new(locants: Vec<i32>, subst: Substituent) -> SubFragment {
        SubFragment { locants, subst }
    }
}

impl Default for SubFragment {
    fn default() -> Self {
        SubFragment {
            locants: vec![0],
            subst: Substituent::Group(Alkane),
        }
    }
}

impl Display for SubFragment {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            locants(self.locants.clone()).unwrap(),
            self.subst
        )
    }
}

#[derive(Debug, Error, PartialEq)]
pub enum NamingError {
    #[error("A branch was found with too many carbons ({}).", .0)]
    CarbonCount(i32),
    #[error("Found too many occurrences of {} ({}).", NamingError::format(.0), .1)]
    GroupOccurrence(Option<Substituent>, i32),
}

impl NamingError {
    fn format(subst: &Option<Substituent>) -> String {
        let unwrapped = subst.as_ref().expect("substituent should be provided");
        match unwrapped {
            Substituent::Branch(it) => format!("the {} branch", it),
            Substituent::Group(it) => match it {
                Group::Hydrogen | Alkane => panic!("process could not handle: {it}"),
                Group::AcidHalide(_) | Group::Aldehyde | Group::Carboxyl | Group::Nitrile => {
                    panic!("excess of 2-max group: {it}")
                }
                Group::Alkene => "double bonds",
                Group::Alkyne => "triple bonds",
                Group::Bromo => "bromine",
                Group::Chloro => "chlorine",
                Group::Fluoro => "fluorine",
                Group::Iodo => "iodine",
                Group::Hydroxyl => "hydroxyl",
                Group::Carbonyl => "carbonyl",
                Group::Ester => "ester",
                Group::Ether => "ether",
                Group::Amide => "amide",
                Group::Amine => "amine",
            }
            .to_string(),
        }
    }

    fn add_substituent(
        e: Result<String, NamingError>,
        subst: &Substituent,
    ) -> Result<String, NamingError> {
        e.map_err(|e| {
            if let NamingError::GroupOccurrence(_, count) = e {
                NamingError::GroupOccurrence(Some(subst.clone()), count)
            } else {
                panic!("found unexpected error from locants")
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Group::{Bromo, Carbonyl, Chloro, Hydroxyl, Iodo};
    use std::str::FromStr;

    #[test]
    fn prefix_creates_correct_strings() {
        let branch =
            Branch::from_str("0: bromo, iodo; 1: oxo, hydroxy; 2: bromo, hydroxy").unwrap();
        let collection = SubFragmentCollection::new(branch);
        let str = prefix(collection.secondary_group_fragments()).unwrap();

        assert_eq!(str, "1,3-dibromo-2,3-dihydroxy-1-iodo")
    }

    #[test]
    fn collect_groups_aggregates_correctly() {
        let branch =
            Branch::from_str("0: bromo, iodo; 1: oxo, hydroxy; 2: bromo, hydroxy").unwrap();
        let groups = SubFragmentCollection::new(branch);
        let expected = vec![
            SubFragment::new(vec![0, 2], Substituent::Group(Bromo)),
            SubFragment::new(vec![0], Substituent::Group(Iodo)),
            SubFragment::new(vec![1], Substituent::Group(Carbonyl)),
            SubFragment::new(vec![1, 2], Substituent::Group(Hydroxyl)),
        ];

        assert_eq!(groups.collection, expected);
    }

    #[test]
    fn primary_group_doesnt_ignore_halogens() {
        let collection = SubFragmentCollection {
            collection: vec![
                SubFragment::new(vec![0], Substituent::Group(Bromo)),
                SubFragment::new(vec![0], Substituent::Group(Chloro)),
            ],
        };
        assert_eq!(collection.primary_group(), Alkane);
    }

    #[test]
    fn locants_formats_correctly() {
        let str = locants(vec![0, 1, 2]).unwrap();

        assert_eq!(str, "1,2,3-tri");
    }

    #[test]
    fn locants_no_prefix_for_1() {
        let str = locants(vec![2]).unwrap();

        assert_eq!(str, "3-")
    }
}
