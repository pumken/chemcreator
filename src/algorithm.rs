//! # Algorithm
//!
//! The `algorithm` module contains the functions needed to find the name of an arbitrary
//! organic molecule.

use ruscii::spatial::Vec2;
use thiserror::Error;
use crate::algorithm::InvalidGraphError::{Cycle, Discontinuity};
use crate::grid::{GridState, Pointer};
use crate::molecule::{Atom, Cell, Element, Group};
use crate::nested_vec;

enum Substituent {
    Branch(Branch),
    Group(Group),
}

struct Branch {
    cells: Vec<Cell>,
    chain: Vec<Substituent>,
}

/// Determines the name of the molecule on the given `graph`.
///
/// ## Errors
///
/// If the molecule on the given `graph` is discontinuous, cyclic, or contains invalid bonding,
/// an [`InvalidGraphError`] will be returned.
pub fn name_molecule(graph: &GridState) -> Result<String, InvalidGraphError> {
    let cells = graph.find_all(|cell| cell.is_atom())
        .iter()
        .map(|&cell| if let Cell::Atom(it) = cell {
            it
        } else {
            panic!("is_atom check failed in name_molecule")
        })
        .collect();

    if graph.is_empty() {
        return Ok("".to_string());
    }
    check_structure(graph)?;
    check_valence(cells, graph)?;
    let chain = continue_chain(graph);
    // let group_indexed_chain = link_groups();
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
    return match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    };
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

/// Checks if the molecule on the [`GridState`] is contiguous.
///
/// ## Errors
///
/// If [`GridState`] contains an invalid molecule, [`Discontinuity`] or [`Cycle`] is returned.
/// An empty [`GridState`] does not return an error.
fn check_structure(graph: &GridState) -> Result<(), InvalidGraphError> {
    let starting_cell = match graph.find(|cell| match cell {
        Cell::None(_) => false,
        _ => true,
    }) {
        Some(it) => it,
        None => return Ok(()),
    };
    let filled_pos_directions = graph.filled_cells()
        .iter()
        .map(|&cell| cell.pos())
        .collect::<Vec<Vec2>>();
    let connectivity = get_connected_cells(starting_cell.pos(), graph)?;

    for pos in filled_pos_directions {
        if match graph.get(pos).unwrap() {
            Cell::Atom(_) | Cell::Bond(_) => true,
            Cell::None(_) => false,
        } != connectivity[pos.x as usize][pos.y as usize] {
            return Err(Discontinuity)
        }
    }

    Ok(())
}

/// Returns all [`Vec2`]s that are connected to the given `pos`.
///
/// ## Panics
///
/// This function panics if the given `pos` is not a valid point on the given `graph` or if the
/// given `pos` on the `graph` is a [`Cell::None`].
///
/// ## Errors
///
/// If this function traverses the molecule and finds that it is not simply connected, [`Cycle`]
/// will be returned.
fn get_connected_cells(pos: Vec2, graph: &GridState) -> Result<Vec<Vec<bool>>, InvalidGraphError> {
    match graph.get(pos).expect("pos should be a valid point on the graph.") {
        Cell::None(_) => panic!("Passed empty cell ({}, {}) to get_connected_cells", pos.x, pos.y),
        _ => {}
    }

    let mut searched_points = nested_vec![graph.size.x; graph.size.y; false];

    fn accumulate_components(
        pos: Vec2,
        previous_pos: Option<Vec2>,
        accumulator: &mut Vec<Vec<bool>>,
        graph: &GridState
    ) -> Result<(), InvalidGraphError> {
        accumulator[pos.x as usize][pos.y as usize] = true;
        let ptr = Pointer { graph, pos };
        for cell in ptr.connected() {
            match previous_pos {
                Some(it) => if cell.pos() == it {
                    continue
                },
                None => {}
            }
            match cell {
                Cell::None(_) => {}
                _ => if !accumulator[cell.pos().x as usize][cell.pos().y as usize] {
                    accumulate_components(cell.pos(), Some(pos), accumulator, graph)?
                } else {
                    return Err(Cycle)
                }
            }
        }
        Ok(())
    }

    accumulate_components(pos, None, &mut searched_points, graph)?;
    Ok(searched_points)
}

/// Returns `Ok` if all of the valence shells of the given [`Cell::Atom`]s
/// are filled, i.e., that the sum of bond orders across all of their bonds is equivalent to their
/// [`Element::bond_number`].
///
/// ## Errors
///
/// Returns an [`InvalidGraphError::OverfilledValence`] or [`InvalidGraphError::UnfilledValence`]
/// for the first cell for which its valence shell is not correctly filled.
fn check_valence(atoms: Vec<&Atom>, graph: &GridState) -> Result<(), InvalidGraphError> {
    for atom in atoms {
        let ptr = Pointer { graph, pos: atom.pos };
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