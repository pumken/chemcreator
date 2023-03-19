//! # Naming
//!
//! The `naming` module contains functions that makes final validations and finds the name of
//! an organic molecule as a [`Branch`].

use crate::groups::Fallible;
use crate::groups::InvalidGraphError::Other;
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
    let branch = groups::link_groups(graph, chain)?;
    // group_indexed_chain.check_chain_index()

    process_name(&branch).map_err(|e| Other(e.to_string()))
}

fn process_name(branch: &Branch) -> Result<String, NamingError> {
    let collection = GroupCollection::new(branch);

    Ok(format!(
        "{}{}{}{}",
        prefix(collection.secondary_group_fragments())?,
        major_numeric(branch.chain.len() as i32)?,
        bonding(collection.chain_group_fragments())?,
        suffix(collection.primary_group_fragment())?,
    ))
}

fn prefix(mut fragments: Vec<GroupFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.group.to_string());
    let out = fragments
        .into_iter()
        .map(|fragment| format!("{}{}", locants(fragment.locants).unwrap(), fragment.group))
        .collect::<Vec<String>>()
        .join("-");
    Ok(out)
}

fn bonding(mut fragments: Vec<GroupFragment>) -> Result<String, NamingError> {
    fragments.sort_by_key(|fragment| fragment.group.to_string());
    if fragments.is_empty() {
        Ok("an".to_string())
    } else {
        let mut out = fragments
            .into_iter()
            .map(|fragment| format!("{}{}", locants(fragment.locants).unwrap(), fragment.group))
            .collect::<Vec<String>>()
            .join("-");
        out.insert(0, '-');
        Ok(out)
    }
}

fn suffix(fragment: GroupFragment) -> Result<String, NamingError> {
    let locations = locants(fragment.locants)?;
    let suffix = match fragment.group {
        Group::Carboxyl => return Ok("oic acid".to_string()),
        Group::Carbonyl => "one",
        Group::Hydroxyl => "ol",
        _ => return Ok("e".to_string()),
    };

    Ok(format!("-{locations}{suffix}"))
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
        4 => "tetra",
        5 => "penta",
        6 => "hexa",
        7 => "hepta",
        8 => "octa",
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
        _ => return Err(NamingError::CarbonCount(value)),
    };
    Ok(out)
}

#[derive(Clone, Debug, PartialEq)]
pub struct GroupCollection {
    pub collection: Vec<GroupFragment>,
}

impl GroupCollection {
    pub fn new(branch: &Branch) -> GroupCollection {
        let mut out = GroupCollection { collection: vec![] };

        branch.groups.iter().enumerate().for_each(|(index, link)| {
            link.iter().for_each(|it| {
                if let Substituent::Group(group) = it {
                    out.push_group(group.to_owned(), index as i32);
                }
            })
        });

        out
    }

    fn push_group(&mut self, group: Group, index: i32) {
        let item = self.collection.iter_mut().find(|it| it.group == group);

        match item {
            Some(it) => it.locants.push(index),
            None => self.collection.push(GroupFragment::new(vec![index], group)),
        }
    }

    pub fn primary_group(&self) -> Option<Group> {
        self.collection
            .iter()
            .map(|fragment| fragment.group)
            .max_by_key(|&group| group.priority())
    }

    pub fn primary_group_fragment(&self) -> GroupFragment {
        self.collection
            .iter()
            .find(|&fragment| fragment.group == self.primary_group().unwrap())
            .map_or_else(GroupFragment::default, GroupFragment::to_owned)
    }

    pub fn secondary_group_fragments(&self) -> Vec<GroupFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.group != self.primary_group().unwrap())
            .filter(|&fragment| !fragment.group.is_chain_group())
            .map(GroupFragment::to_owned)
            .collect::<Vec<GroupFragment>>()
    }

    pub fn chain_group_fragments(&self) -> Vec<GroupFragment> {
        self.collection
            .iter()
            .filter(|&fragment| fragment.group.is_chain_group())
            .map(GroupFragment::to_owned)
            .collect::<Vec<GroupFragment>>()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GroupFragment {
    locants: Vec<i32>,
    group: Group,
}

impl GroupFragment {
    pub fn new(locants: Vec<i32>, group: Group) -> GroupFragment {
        GroupFragment { locants, group }
    }
}

impl Default for GroupFragment {
    fn default() -> Self {
        GroupFragment {
            locants: vec![0],
            group: Group::Alkane,
        }
    }
}

impl Display for GroupFragment {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            locants(self.locants.clone()).unwrap(),
            self.group
        )
    }
}

#[derive(Debug, Error, PartialEq)]
enum NamingError {
    #[error("A branch was found with too many carbons ({}).", .0)]
    CarbonCount(i32),
    #[error("Found too many occurrences of the {:?} group ({}).", .0, .1)]
    GroupOccurrence(Option<Group>, i32),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Group::{Bromo, Carbonyl, Hydroxyl, Iodo};
    use std::str::FromStr;

    #[test]
    fn prefix_creates_correct_strings() {
        let branch =
            Branch::from_str("0: bromo, iodo; 1: oxo, hydroxyl; 2: bromo, hydroxyl").unwrap();
        let collection = GroupCollection::new(&branch);
        let str = prefix(collection.secondary_group_fragments()).unwrap();

        assert_eq!(str, "1,3-dibromo-2,3-dihydroxyl-1-iodo")
    }

    #[test]
    fn collect_groups_aggregates_correctly() {
        let branch =
            Branch::from_str("0: bromo, iodo; 1: oxo, hydroxyl; 2: bromo, hydroxyl").unwrap();
        let groups = GroupCollection::new(&branch);
        let expected = vec![
            GroupFragment::new(vec![0, 2], Bromo),
            GroupFragment::new(vec![0], Iodo),
            GroupFragment::new(vec![1], Carbonyl),
            GroupFragment::new(vec![1, 2], Hydroxyl),
        ];

        assert_eq!(groups.collection, expected);
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
