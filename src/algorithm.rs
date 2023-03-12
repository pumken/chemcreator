//! # Algorithm
//!
//! The `algorithm` module contains the functions needed to find the name of an arbitrary
//! organic molecule.

use ruscii::spatial::Vec2;
use thiserror::Error;
use crate::grid::{GridState, Pointer};
use crate::molecule::{Cell, Element, Group};

enum Substituent {
    Branch(Branch),
    Group(Group),
}

struct Branch {
    cells: Vec<Cell>,
    chain: Vec<Substituent>,
}

pub fn name_molecule(graph: &GridState) -> Result<String, InvalidGraphError> {
    check_valence(graph)?;
    let chain = continue_chain(graph);
    // let group_indexed_chain = link_groups();
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
    return match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    }
}

/// A recursive function that traverses sequential carbon bonds and returns a [`Vec`] of carbon
/// chains as [`Vec`]s of [`Cell::Atom`]s.
pub(crate) fn continue_chain(graph: &GridState) {
    let endpoints = endpoint_carbons(graph);

}

pub(crate) fn endpoint_carbons(graph: &GridState) -> Result<Vec<&Cell>, InvalidGraphError> {
    let all_carbons = graph.find_all(|cell| {
        match cell {
            Cell::Atom(atom) => match atom.element {
                Element::C => true,
                _ => false
            }
            _ => false
        }
    });
    let mut out = vec![];

    for carbon in all_carbons {
        let ptr = Pointer::new(carbon, graph);
        if ptr.bonded_carbons()? <= 1 {
            out.push(carbon);
        }
    }
    Ok(out)
}

/// Retrieves all the atoms in the [`GridState`] and returns `Ok` if all of their valence shells
/// are filled, i.e., that the sum of bond orders across all of their bonds is equivalent to their
/// [`Element::bond_number`].
///
/// ## Errors
///
/// Returns an [`InvalidGraphError::OverfilledValence`] or [`InvalidGraphError::UnfilledValence`]
/// for the first cell for which its valence shell is not correctly filled.
fn check_valence(graph: &GridState) -> Result<(), InvalidGraphError> {
    let cells = graph.find_all(|cell| cell.is_atom());

    for cell in cells {
        if let Cell::Atom(atom) = cell {
            let ptr = Pointer::new(cell, graph);
            let bond_count = match ptr.bond_count() {
                Ok(it) => it,
                Err(it) => panic!("{}", it)
            };
            if bond_count > atom.element.bond_number() {
                return Err(InvalidGraphError::OverfilledValence(atom.pos));
            } else if bond_count < atom.element.bond_number() {
                return Err(InvalidGraphError::UnfilledValence(atom.pos));
            }
        }
    }
    Ok(())
}

#[derive(Error, Debug)]
pub enum InvalidGraphError {
    #[error("Molecule is not continuous.")]
    Discontinuity,
    #[error("Molecule is cyclic.")]
    Cycle,
    #[error("Cell at ({}, {}) is missing bonds.", .0.x, .0.y)]
    UnfilledValence(Vec2),
    #[error("Cell at ({}, {}) has too many bonds.", .0.x, .0.y)]
    OverfilledValence(Vec2),
    #[error("Bond at ({}, {}) is incomplete.", .0.x, .0.y)]
    IncompleteBond(Vec2),
    #[error("This combination of groups is not supported.")]
    UnsupportedGroups,
}