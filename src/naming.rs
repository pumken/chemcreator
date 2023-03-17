//! # Naming
//!
//! The `naming` module contains functions that makes final validations and finds the name of
//! an organic molecule as a [`Branch`].

use crate::groups::{Fallible, InvalidGraphError};
use crate::molecule::{Branch, Cell};
use crate::spatial::GridState;
use crate::{chain, groups, validation};

fn process_name(branch: &Branch) -> Fallible<String> {
    Ok("".to_string())
}

/// Determines the name of the molecule on the given `graph`.
///
/// ## Errors
///
/// If the molecule on the given `graph` is discontinuous, cyclic, or contains invalid bonding,
/// an [`InvalidGraphError`] will be returned.
pub fn name_molecule(graph: &GridState) -> Fallible<String> {
    let cells = graph.find_all(|cell| cell.is_atom())
        .iter()
        .map(|&cell| if let Cell::Atom(it) = cell {
            it
        } else {
            panic!("is_atom check failed in name_molecule")
        })
        .collect();

    // Initial checks
    if graph.is_empty() { return Ok("".to_string()); }
    validation::check_structure(graph)?;
    validation::check_valence(cells, graph)?;

    // Preliminary chain
    let all_chains = chain::get_all_chains(graph)?;
    let chain = chain::longest_chain(all_chains)?;

    // Group-linked branch
    let branch = groups::link_groups(graph, chain)?;
    // group_indexed_chain.check_chain_index()


    let name: Fallible<String> = process_name(&branch);
    match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    }
}
