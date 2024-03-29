//! # Validation
//!
//! The `validation` module provides two functions that checks the [`GridState`] for three invalid
//! structures: discontinuities, cycles, and atoms with unfilled or overfilled valence shells.

use crate::chain::get_connected_cells;
use crate::groups::Fallible;
use crate::groups::InvalidGraphError::{self, Discontinuity, Other};
use crate::molecule::Group::{Alkene, Alkyne};
use crate::molecule::{Atom, Branch, Cell, Substituent};
use crate::naming::{prefix, SubFragment, SubFragmentCollection};
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
    let connectivity = get_connected_cells(starting_cell.pos(), graph)?;
    let filled_pos_directions: Vec<Vec2> = graph
        .filled_cells()
        .iter()
        .map(|&cell| cell.pos())
        .collect();

    for pos in filled_pos_directions {
        if graph.get(pos).unwrap().is_not_empty() != connectivity[pos] {
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
pub fn check_valence(atoms: Vec<Atom>, graph: &GridState) -> Fallible<()> {
    for atom in atoms {
        let ptr = Pointer::new(graph, atom.pos);
        let bond_count = match ptr.bond_count() {
            Ok(it) => it,
            Err(it) => panic!("{it}"),
        };

        match bond_count.cmp(&atom.element.bond_number()) {
            Ordering::Less => {
                return Err(InvalidGraphError::UnfilledValence(atom.pos, atom.element));
            }
            Ordering::Greater => {
                return Err(InvalidGraphError::OverfilledValence(atom.pos, atom.element));
            }
            _ => {}
        }
    }

    Ok(())
}

impl Branch {
    /// Checks if chain indexes are in the correct direction.
    pub fn index_corrected(self) -> Fallible<Branch> {
        let original = self.clone();
        let reversed = self.reversed();
        let original_collection = SubFragmentCollection::new(original.clone());
        let reversed_collection = SubFragmentCollection::new(reversed.clone());

        match cmp(
            original_collection.primary_group_fragment(),
            reversed_collection.primary_group_fragment(),
        ) {
            Ordering::Less => return Ok(original),
            Ordering::Greater => return Ok(reversed),
            Ordering::Equal => {}
        }

        let original_multiple_bonds = original_collection.unsaturated_group_fragments();
        let reversed_multiple_bonds = reversed_collection.unsaturated_group_fragments();

        let original_alkenes = original_multiple_bonds
            .iter()
            .find(|&sub| matches!(sub.subst, Substituent::Group(Alkene)));
        let reversed_alkenes = reversed_multiple_bonds
            .iter()
            .find(|&sub| matches!(sub.subst, Substituent::Group(Alkene)));

        if let (Some(org), Some(rev)) = (original_alkenes, reversed_alkenes) {
            match cmp(org.to_owned(), rev.to_owned()) {
                Ordering::Less => return Ok(original),
                Ordering::Greater => return Ok(reversed),
                Ordering::Equal => {}
            }
        }

        let original_alkynes = original_multiple_bonds
            .iter()
            .find(|&sub| matches!(sub.subst, Substituent::Group(Alkyne)));
        let reversed_alkynes = reversed_multiple_bonds
            .iter()
            .find(|&sub| matches!(sub.subst, Substituent::Group(Alkyne)));

        if let (Some(org), Some(rev)) = (original_alkynes, reversed_alkynes) {
            match cmp(org.to_owned(), rev.to_owned()) {
                Ordering::Less => return Ok(original),
                Ordering::Greater => return Ok(reversed),
                Ordering::Equal => {}
            }
        }

        let original_prefixes_str = prefix(original_collection.secondary_group_fragments())
            .map_err(|e| Other(e.to_string()))?;
        let reversed_prefixes_str = prefix(reversed_collection.secondary_group_fragments())
            .map_err(|e| Other(e.to_string()))?;

        match original_prefixes_str.cmp(&reversed_prefixes_str) {
            Ordering::Equal | Ordering::Less => Ok(original),
            Ordering::Greater => Ok(reversed),
        }
    }
}

fn cmp(first: SubFragment, second: SubFragment) -> Ordering {
    let mut first_locants = first.locants;
    first_locants.sort();
    let mut second_locants = second.locants;
    second_locants.sort();

    for (index, locant) in first_locants.iter().enumerate() {
        match locant.cmp(&second_locants[index]) {
            Ordering::Equal => continue,
            Ordering::Less => return Ordering::Less,
            Ordering::Greater => return Ordering::Greater,
        }
    }

    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::{Single, Triple};
    use crate::molecule::Element::{C, H, O};
    use crate::molecule::Group::Alkane;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn check_structure_validates() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [1, 0; A(H)],
            [0, 1; A(H)],
            [2, 1; A(H)],
            [1, 2; A(H)],
        );
        let ok = check_structure(&graph);

        assert!(matches!(ok, Ok(_)));
    }

    #[test]
    fn check_structure_returns_discontinuity() {
        let graph = graph_with!(3, 1,
            [0, 0; A(C)],
            [2, 0; A(C)],
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
            [1, 2; B(Single)],
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
            [0, 1; B(Triple)],
        );
        let input_atoms = vec![
            graph.get(Vec2::xy(1, 1)).unwrap(),
            graph.get(Vec2::xy(3, 1)).unwrap(),
            graph.get(Vec2::xy(5, 1)).unwrap(),
        ]
        .iter()
        .map(|&cell| cell.unwrap_atom())
        .collect::<Vec<Atom>>();

        assert!(matches!(check_valence(input_atoms, &graph), Ok(_)));
    }

    #[test]
    fn check_valence_returns_unfilled() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [0, 1; B(Single)],
            [2, 1; B(Single)],
            [1, 2; B(Single)],
        );
        let input_atom = graph.get(Vec2::xy(1, 1)).unwrap().unwrap_atom();
        let err = check_valence(vec![input_atom], &graph);

        assert_eq!(
            err,
            Err(InvalidGraphError::UnfilledValence(Vec2::xy(1, 1), C))
        );
    }

    #[test]
    fn check_valence_returns_overfilled() {
        let graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [0, 1; B(Triple)],
            [2, 1; B(Triple)],
            [1, 2; B(Triple)],
        );
        let input_atom = graph.get(Vec2::xy(1, 1)).unwrap().unwrap_atom();
        let err = check_valence(vec![input_atom], &graph);

        assert_eq!(
            err,
            Err(InvalidGraphError::OverfilledValence(Vec2::xy(1, 1), C))
        );
    }

    #[test]
    fn cmp_locants() {
        let a = SubFragment::new(vec![1, 1, 3], Substituent::Group(Alkane));
        let b = SubFragment::new(vec![2, 3, 4], Substituent::Group(Alkane));
        let c = SubFragment::new(vec![1, 2, 2], Substituent::Group(Alkane));
        let d = SubFragment::new(vec![1, 2, 2], Substituent::Group(Alkane));

        assert_eq!(cmp(a.clone(), b.clone()), Ordering::Less);
        assert_eq!(cmp(b, c.clone()), Ordering::Greater);
        assert_eq!(cmp(a, c.clone()), Ordering::Less);
        assert_eq!(cmp(c, d), Ordering::Equal);
    }
}
