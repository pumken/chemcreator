//! # Naming
//!
//! The `naming` module contains functions that makes final validations and finds the name of
//! an organic molecule as a [`Branch`].

use crate::groups::Fallible;
use crate::molecule::{Branch, Cell, Group, Substituent};
use crate::spatial::GridState;
use crate::{chain, groups, validation};
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

    process_name(&branch)
}

fn process_name(branch: &Branch) -> Fallible<String> {
    let max_priority = branch
        .groups
        .iter()
        .flatten()
        .filter(|&subst| matches!(subst, Substituent::Group(_)))
        .map(|subst| {
            if let Substituent::Group(it) = subst {
                it
            } else {
                panic!("Non-group found after filter in process_name")
            }
        })
        .max_by_key(|&group| group.priority());

    let name = match max_priority {
        None => Ok("".to_string()),
        Some(Group::Carboxyl) => carboxylic_acid(branch),
        Some(Group::Carbonyl) => aldehyde_ketone(branch),
        Some(Group::Hydroxyl) => alcohol(branch),
        _ => Ok("".to_string()),
    };
    Ok(format!(
        "{}{}",
        prefix(branch, max_priority.unwrap()).unwrap(),
        name.unwrap()
    ))
}

fn carboxylic_acid(branch: &Branch) -> Result<String, NamingError> {
    let numeric = major_numeric_prefix(branch.chain.len() as i32)?;
    Ok(format!("{numeric}anoic acid"))
}

fn aldehyde_ketone(branch: &Branch) -> Result<String, NamingError> {
    let numeric = major_numeric_prefix(branch.chain.len() as i32)?;
    Ok(format!("{numeric}anone"))
}

fn alcohol(branch: &Branch) -> Result<String, NamingError> {
    let numeric = major_numeric_prefix(branch.chain.len() as i32)?;
    Ok(format!("{numeric}anol"))
}

fn prefix(branch: &Branch, main_group: &Group) -> Result<String, NamingError> {
    let mut prefixes: Vec<(Vec<usize>, Group)> = vec![];
    prefixes.sort_by(|a, b| a.1.to_string().cmp(&b.1.to_string()));
    let out = prefixes
        .into_iter()
        .map(|tuple| group_numeric(tuple.0, tuple.1))
        .collect::<Result<Vec<String>, NamingError>>()?
        .join("-");
    Ok(out)
}

fn group_numeric(mut locations: Vec<usize>, group: Group) -> Result<String, NamingError> {
    locations.sort();
    let numbers = locations
        .iter()
        .map(|&it| it.to_string())
        .collect::<Vec<String>>()
        .join(",");
    let out = format!(
        "{numbers}-{}{}",
        minor_numeric_prefix(&group, locations.len() as i32)?,
        group
    );
    Ok(out)
}

fn minor_numeric_prefix(group: &Group, value: i32) -> Result<&'static str, NamingError> {
    let out = match value {
        1 => "",
        2 => "di",
        3 => "tri",
        4 => "tetra",
        5 => "penta",
        6 => "hexa",
        7 => "hepta",
        8 => "octa",
        _ => return Err(NamingError::GroupOccurrence(group.to_owned(), value)),
    };
    Ok(out)
}

fn major_numeric_prefix(value: i32) -> Result<&'static str, NamingError> {
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

#[derive(Debug, Error, PartialEq)]
enum NamingError {
    #[error("A branch was found with too many carbons ({}).", .0)]
    CarbonCount(i32),
    #[error("Found too many occurrences of the {} group ({}).", .0, .1)]
    GroupOccurrence(Group, i32),
}

#[cfg(test)]
mod tests {
    use super::*;
}
