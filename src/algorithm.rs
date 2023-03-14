//! # Algorithm
//!
//! The `algorithm` module contains the functions needed to find the name of an arbitrary
//! organic molecule.

use std::cmp::Ordering;
use ruscii::spatial::Vec2;
use thiserror::Error;
use crate::algorithm::InvalidGraphError::{Cycle, Discontinuity, Other};
use crate::spatial::GridState;
use crate::molecule::{Atom, Cell, Element, Group};
use crate::nested_vec;
use crate::pointer::Pointer;

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

    // Initial checks
    if graph.is_empty() { return Ok("".to_string()); }
    check_structure(graph)?;
    check_valence(cells, graph)?;

    // Preliminary chain
    let all_chains = get_all_chains(graph)?;
    let chain = longest_chain(all_chains);
    // let group_indexed_chain = link_groups();
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
    match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    }
}

pub(crate) fn debug_chain(graph: &GridState) -> Result<Vec<Atom>, InvalidGraphError> {
    let all_chains = get_all_chains(graph)?;
    Ok(longest_chain(all_chains)?)
}

/// Gets the longest of the given [`Vec`] of chains, assuming that it is non-empty.
///
/// ## Errors
///
/// If there are no `chains`, this function will return an `Err`.
pub(crate) fn longest_chain(chains: Vec<Vec<Atom>>) -> Result<Vec<Atom>, InvalidGraphError> {
    Ok(match chains.iter()
        .max_by(|&a, &b| a.len().cmp(&b.len())) {
        None => return Err(Other("No carbon chain found.".to_string())),  // FIXME this is returned when a bond is placed at the edge
        Some(it) => it
    }.to_owned())
}

pub(crate) fn get_all_chains(graph: &GridState) -> Result<Vec<Vec<Atom>>, InvalidGraphError> {
    let endpoints = endpoint_carbons(graph)?.iter()
        .map(|&cell| match cell {
            Cell::Atom(it) => it.to_owned(),
            _ => panic!("endpoint_carbons returned non-atom cell")
        })
        .collect::<Vec<Atom>>();
    // Nested Vec hell
    let mut out: Vec<Vec<Atom>> = vec![];

    for endpoint in endpoints {
        out.extend(endpoint_head_chains(endpoint.to_owned(), graph)?);
    }
    Ok(out)
}

/// Takes a given `endpoint` and returns all continuous carbon chains starting at it.
///
/// ## Errors
///
/// Returns [`InvalidGraphError`] if any invalid structures are found while traversing the graph.
fn endpoint_head_chains(endpoint: Atom, graph: &GridState) -> Result<Vec<Vec<Atom>>, InvalidGraphError> {
    let mut accumulator = vec![vec![]];

    accumulate_carbons(
        endpoint.pos,
        None,
        0usize,
        &mut accumulator,
        graph,
    )?;

    Ok(accumulator)
}

/// Adds the given `pos` to the branch in the `accumulator` with the given `branch_index`. After,
/// this function is called on the unvisited carbon neighbors of the atom at the given `pos`.
///
/// ## Panics
///
/// This function panics if the given `pos` is not a [`Cell::Atom`].
///
/// ## Errors
///
/// If the given `graph` contains an invalid molecule, [`Discontinuity`] or [`Cycle`] is returned.
fn accumulate_carbons(
    pos: Vec2,
    previous_pos: Option<Vec2>,
    branch_index: usize,
    accumulator: &mut Vec<Vec<Atom>>,
    graph: &GridState,
) -> Result<(), InvalidGraphError> {
    let next_carbons = next_carbons(pos, previous_pos, graph)?;

    accumulator[branch_index].push(match graph.get(pos) {
        Ok(Cell::Atom(it)) => it.to_owned(),
        _ => panic!("Non-atom or invalid cell passed to accumulate_carbons")
    });

    let new_branches = create_branches(
        accumulator,
        branch_index,
        match next_carbons.len() {
            0 => return Ok(()),
            it => it - 1
        }
    );
    for (i, carbon) in next_carbons.iter().enumerate() {
        accumulate_carbons(
            carbon.pos,
            Some(pos),
            new_branches[i],
            accumulator,
            graph
        )?
    }

    Ok(())
}

/// Creates `count` clones of the branch at the given `branch_index` in the `accumulator`. The
/// returned [`Vec`] contains the indexes of the branch and its copies.
fn create_branches(accumulator: &mut Vec<Vec<Atom>>, branch_index: usize, count: usize) -> Vec<usize> {
    let mut out = vec![branch_index];

    for _ in 0usize..count {
        accumulator.push(accumulator[branch_index].clone());
        out.push(accumulator.len() - 1)
    };

    out
}

/// Gets the next bonded carbons. Does not include the cell at the `previous_pos` if there is one.
///
/// ## Errors
///
/// If one of the bonds to the current cell is found to be dangling, an
/// [`IncompleteBond`] will be returned.
fn next_carbons(pos: Vec2, previous_pos: Option<Vec2>, graph: &GridState) -> Result<Vec<Atom>, InvalidGraphError> {
    let ptr = Pointer { graph, pos };
    let mut out = ptr.bonded_carbons()?;

    if let Some(it) = previous_pos {
        out.retain(|atom| atom.pos != it);
    }
    Ok(out)
}

/// Returns references to all cells containing endpoint carbon atoms, i.e. those that have exactly
/// one or no carbon neighbors.
///
/// ## Errors
///
/// If one of the bonds to the current cell is found to be dangling, an
/// [`IncompleteBond`] will be returned.
pub(crate) fn endpoint_carbons(graph: &GridState) -> Result<Vec<&Cell>, InvalidGraphError> {
    let all_carbons = graph.find_all(|cell| {
        match cell {
            Cell::Atom(atom) => matches!(atom.element, Element::C),
            _ => false
        }
    });
    let mut out = vec![];

    for carbon in all_carbons {
        let ptr = Pointer::new(carbon, graph);
        if ptr.bonded_carbon_count()? <= 1 {
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
    let starting_cell = match graph.find(|cell| !matches!(cell, Cell::None(_))) {
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
            return Err(Discontinuity);
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
    if let Cell::None(_) = graph.get(pos).expect("pos should be a valid point on the graph.") {
        panic!("Passed empty cell ({}, {}) to get_connected_cells", pos.x, pos.y)
    }

    let mut searched_points = nested_vec![graph.size.x; graph.size.y; false];

    fn accumulate_components(
        pos: Vec2,
        previous_pos: Option<Vec2>,
        accumulator: &mut Vec<Vec<bool>>,
        graph: &GridState,
    ) -> Result<(), InvalidGraphError> {
        accumulator[pos.x as usize][pos.y as usize] = true;
        let ptr = Pointer { graph, pos };
        for cell in ptr.connected() {
            if let Some(it) = previous_pos {
                if cell.pos() == it {
                    continue;
                }
            }
            match cell {
                Cell::None(_) => {}
                _ => if !accumulator[cell.pos().x as usize][cell.pos().y as usize] {
                    accumulate_components(cell.pos(), Some(pos), accumulator, graph)?
                } else {
                    return Err(Cycle);
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
        match bond_count.cmp(&atom.element.bond_number()) {
            Ordering::Less => return Err(InvalidGraphError::UnfilledValence(atom.pos)),
            Ordering::Greater => return Err(InvalidGraphError::OverfilledValence(atom.pos)),
            _ => {}
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
    #[error("{}", .0)]
    Other(String)
}
