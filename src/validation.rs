use ruscii::spatial::Vec2;
use std::cmp::Ordering;
use crate::chain;
use crate::algorithm::{Fallible, InvalidGraphError};
use crate::algorithm::InvalidGraphError::Discontinuity;
use crate::molecule::{Atom, Cell};
use crate::pointer::Pointer;
use crate::spatial::GridState;

/// Checks if the molecule on the [`GridState`] is contiguous.
///
/// ## Errors
///
/// If [`GridState`] contains an invalid molecule, [`Discontinuity`] or [`Cycle`] is returned.
/// An empty [`GridState`] does not return an error.
pub fn check_structure(graph: &GridState) -> Fallible<()> {
    let starting_cell = match graph.find(|cell| !matches!(cell, Cell::None(_))) {
        Some(it) => it,
        None => return Ok(()),
    };
    let filled_pos_directions = graph.filled_cells()
        .iter()
        .map(|&cell| cell.pos())
        .collect::<Vec<Vec2>>();
    let connectivity = chain::get_connected_cells(starting_cell.pos(), graph)?;

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

/// Returns `Ok` if all of the valence shells of the given [`Cell::Atom`]s
/// are filled, i.e., that the sum of bond orders across all of their bonds is equivalent to their
/// [`Element::bond_number`].
///
/// ## Errors
///
/// Returns an [`InvalidGraphError::OverfilledValence`] or [`InvalidGraphError::UnfilledValence`]
/// for the first cell for which its valence shell is not correctly filled.
pub fn check_valence(atoms: Vec<&Atom>, graph: &GridState) -> Fallible<()> {
    for atom in atoms {  // this function eventually could be removed and incorporated into bonded() ?
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
