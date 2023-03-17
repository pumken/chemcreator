//! # Validation
//!
//! The `validation` module provides two functions that checks the [`GridState`] for three invalid
//! structures: discontinuities, cycles, and atoms with unfilled or overfilled valence shells.

use crate::chain;
use crate::groups::InvalidGraphError::Discontinuity;
use crate::groups::{Fallible, InvalidGraphError};
use crate::molecule::{Atom, Cell};
use crate::pointer::Pointer;
use crate::spatial::GridState;
use ruscii::spatial::Vec2;
use std::cmp::Ordering;

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
    let filled_pos_directions = graph
        .filled_cells()
        .iter()
        .map(|&cell| cell.pos())
        .collect::<Vec<Vec2>>();
    let connectivity = chain::get_connected_cells(starting_cell.pos(), graph)?;

    for pos in filled_pos_directions {
        if match graph.get(pos).unwrap() {
            Cell::Atom(_) | Cell::Bond(_) => true,
            Cell::None(_) => false,
        } != connectivity[pos.x as usize][pos.y as usize]
        {
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
    for atom in atoms {
        // this function eventually could be removed and incorporated into bonded() ?
        let ptr = Pointer {
            graph,
            pos: atom.pos,
        };
        let bond_count = match ptr.bond_count() {
            Ok(it) => it,
            Err(it) => panic!("{}", it),
        };
        match bond_count.cmp(&atom.element.bond_number()) {
            Ordering::Less => return Err(InvalidGraphError::UnfilledValence(atom.pos)),
            Ordering::Greater => return Err(InvalidGraphError::OverfilledValence(atom.pos)),
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::{Single, Triple};
    use crate::molecule::Element::{C, H, O};
    use crate::test_utils::unwrap_atom;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn check_structure_validates() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [1, 0; A(H)],
            [0, 1; A(H)],
            [2, 1; A(H)],
            [1, 2; A(H)]
        );
        let ok = check_structure(&graph);

        assert!(matches!(ok, Ok(_)));
    }

    #[test]
    fn check_structure_returns_discontinuity() {
        let graph = graph_with!(3, 1,
            [0, 0; A(C)],
            [2, 0; A(C)]
        );

        let err = check_structure(&graph);

        assert!(matches!(err, Err(InvalidGraphError::Discontinuity)));
    }

    #[test]
    fn check_structure_returns_cycle() {
        let graph = graph_with!(3, 3,
            [0, 0; A(C)],
            [2, 0; A(C)],
            [0, 2; A(C)],
            [2, 2; A(C)],
            [1, 0; B(Single)],
            [0, 1; B(Single)],
            [2, 1; B(Single)],
            [1, 2; B(Single)]
        );
        let err = check_structure(&graph);

        assert!(matches!(err, Err(InvalidGraphError::Cycle)));
    }

    #[test]
    fn check_valence_validates() {
        let graph = graph_with!(7, 3,
            [1, 1; A(C)],
            [3, 1; A(O)],
            [5, 1; A(H)],
            [4, 1; B(Single)],
            [2, 1; B(Single)],
            [0, 1; B(Triple)]
        );
        let input_atoms = vec![
            graph.get(Vec2::xy(1, 1)).unwrap(),
            graph.get(Vec2::xy(3, 1)).unwrap(),
            graph.get(Vec2::xy(5, 1)).unwrap(),
        ]
        .iter()
        .map(|&cell| unwrap_atom(cell))
        .collect::<Vec<Atom>>();
        let references = input_atoms.iter().collect::<Vec<&Atom>>();

        assert!(matches!(check_valence(references, &graph), Ok(_)));
    }

    #[test]
    fn check_valence_returns_unfilled() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [0, 1; B(Single)],
            [2, 1; B(Single)],
            [1, 2; B(Single)]
        );
        let input_atom = unwrap_atom(graph.get(Vec2::xy(1, 1)).unwrap());
        let err = check_valence(vec![&input_atom], &graph);

        assert_eq!(err, Err(InvalidGraphError::UnfilledValence(Vec2::xy(1, 1))));
    }

    #[test]
    fn check_valence_returns_overfilled() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [0, 1; B(Triple)],
            [2, 1; B(Triple)],
            [1, 2; B(Triple)]
        );
        let input_atom = unwrap_atom(graph.get(Vec2::xy(1, 1)).unwrap());
        let err = check_valence(vec![&input_atom], &graph);

        assert_eq!(
            err,
            Err(InvalidGraphError::OverfilledValence(Vec2::xy(1, 1)))
        );
    }
}
