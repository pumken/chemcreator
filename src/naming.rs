//! # Naming
//!
//! The `naming` module contains functions that makes final validations and finds the name of
//! an organic molecule as a [`Branch`].

use crate::groups::Fallible;
use crate::groups::InvalidGraphError::Other;
use crate::molecule::Group::Alkane;
use crate::molecule::{Branch, Cell, Group, Substituent};
use crate::spatial::GridState;
use crate::{chain, groups, validation};
use std::fmt::{Display, Formatter};
use thiserror::Error;

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
                it
            } else {
                panic!("is_atom check failed in name_molecule")
            }
        })
        .collect();

    // Initial checks
    if graph.is_empty() {
        return Ok("".to_string());
    }
    validation::check_structure(graph)?;
    validation::check_valence(cells, graph)?;

    // Preliminary chain
    let all_chains = chain::get_all_chains(graph)?;
    let chain = chain::longest_chain(all_chains)?;

    // Group-linked branch
    let mut branch = groups::link_groups(graph, chain, None)?;
    branch = branch.index_corrected();

    process_name(branch).map_err(|e| Other(e.to_string()))
}

pub(crate) fn process_name(branch: Branch) -> Result<String, NamingError> {
    let len = branch.chain.len() as i32;
    let collection = GroupCollection::new(branch);

    Ok(format!(
        "{}{}{}{}",
        prefix(collection.secondary_group_fragments())?,
        major_numeric(len)?,
        bonding(collection.chain_group_fragments())?,
        suffix(collection.primary_group_fragment())?,
    ))
}

fn prefix(mut fragments: Vec<SubFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.subst.to_string());
    let out = fragments
        .into_iter()
        .map(|fragment| format!("{}{}", locants(fragment.locants).unwrap(), fragment.subst))
        .collect::<Vec<String>>()
        .join("-");
    Ok(out)
}

fn bonding(mut fragments: Vec<SubFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.subst.to_string());
    if fragments.is_empty() {
        Ok("an".to_string())
    } else {
        let mut out = fragments
            .into_iter()
            .map(|fragment| format!("{}{}", locants(fragment.locants).unwrap(), fragment.subst))
            .collect::<Vec<String>>()
            .join("-");
        out.insert(0, '-');
        Ok(out)
    }
}

fn suffix(fragment: SubFragment) -> Result<String, NamingError> {
    if let Substituent::Group(group) = fragment.subst {
        let locations = locants(fragment.locants.clone())?;
        let suffix = match group {
            Group::Carboxyl if fragment.locants.len() == 2 => return Ok("edioic acid".to_string()),
            Group::Carboxyl => return Ok("oic acid".to_string()),
            Group::Carbonyl if fragment.locants == vec![0, fragment.locants.len() as i32 - 1] => {
                return Ok("edial".to_string())
            }
            Group::Carbonyl if fragment.locants == vec![0] => return Ok("al".to_string()),
            Group::Carbonyl => "one",
            Group::Hydroxyl => "ol",
            Group::Nitrile => "onitrile",
            _ => return Ok("e".to_string()),
        };

        Ok(format!("-{locations}{suffix}"))
    } else {
        panic!("branch provided to suffix function")
    }
}

/// Returns a [`String`] containing a locant prefix with a minor numeric prefix appended.
fn locants(mut locations: Vec<i32>) -> Result<String, NamingError> {
    locations.sort();
    let numbers = locations
        .iter()
        .map(|&it| (it + 1).to_string())
        .collect::<Vec<String>>()
        .join(",");
    let out = format!("{numbers}-{}", minor_numeric(locations.len() as i32)?,);
    Ok(out)
}

fn minor_numeric(value: i32) -> Result<&'static str, NamingError> {
    let out = match value {
        1 => "",
        2 => "di",
        3 => "tri",
        4 => "tetr",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        _ => return Err(NamingError::GroupOccurrence(None, value)),
    };
    Ok(out)
}

fn major_numeric(value: i32) -> Result<&'static str, NamingError> {
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

#[derive(Clone, Debug, PartialEq)]
pub struct GroupCollection {
    pub collection: Vec<SubFragment>,
}

impl GroupCollection {
    pub fn new(branch: Branch) -> GroupCollection {
        let mut out = GroupCollection { collection: vec![] };

        branch
            .groups
            .into_iter()
            .enumerate()
            .for_each(|(index, link)| {
                link.into_iter().for_each(|it| {
                    out.push_fragment(it, index as i32);
                })
            });

        out
    }

    fn push_fragment(&mut self, subst: Substituent, index: i32) {
        let item = self.collection.iter_mut().find(|it| it.subst == subst);

        match item {
            Some(it) => it.locants.push(index),
            None => self.collection.push(SubFragment::new(vec![index], subst)),
        }
    }

    pub fn primary_group(&self) -> Group {
        let primary = self
            .collection
            .iter()
            .map(|fragment| fragment.subst.to_owned())
            .filter(|group| group.priority().is_some())
            .max_by_key(|group| group.priority());

        match primary {
            Some(Substituent::Group(it)) => it,
            _ => Alkane,
        }
    }

    pub fn primary_group_fragment(&self) -> SubFragment {
        self.collection
            .iter()
            .find(|&fragment| fragment.subst == Substituent::Group(self.primary_group()))
            .map_or_else(SubFragment::default, SubFragment::to_owned)
    }

    pub fn secondary_group_fragments(&self) -> Vec<SubFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.subst != Substituent::Group(self.primary_group()))
            .filter(|&fragment| !fragment.subst.is_chain_group())
            .map(SubFragment::to_owned)
            .collect::<Vec<SubFragment>>()
    }

    pub fn chain_group_fragments(&self) -> Vec<SubFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.subst.is_chain_group())
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
    #[error("Found too many occurrences of the {:?} group ({}).", .0, .1)]
    GroupOccurrence(Option<Group>, i32),
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
        let collection = GroupCollection::new(branch);
        let str = prefix(collection.secondary_group_fragments()).unwrap();

        assert_eq!(str, "1,3-dibromo-2,3-dihydroxy-1-iodo")
    }

    #[test]
    fn collect_groups_aggregates_correctly() {
        let branch =
            Branch::from_str("0: bromo, iodo; 1: oxo, hydroxy; 2: bromo, hydroxy").unwrap();
        let groups = GroupCollection::new(branch);
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
        let collection = GroupCollection {
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
