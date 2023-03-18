//! # Groups
//!
//! The `groups` module provides functionality for identifying functional groups on a branch.

use crate::chain;
use crate::groups::InvalidGraphError::{Other, UnrecognizedGroup};
use crate::molecule::Group::{Carbonyl, Carboxyl, Hydroxyl};
use crate::molecule::{Atom, BondOrder, Branch, Element, Group, GroupNode, Substituent};
use crate::pointer::Pointer;
use crate::spatial::{FromVec2, GridState};
use ruscii::spatial::{Direction, Vec2};
use thiserror::Error;

/// Generates a [`Branch`] from the given `chain` containing all functional groups attached
/// to it.
pub(crate) fn link_groups(graph: &GridState, chain: Vec<Atom>) -> Fallible<Branch> {
    let mut branch = Branch::new(chain);

    fn accumulate_groups(
        graph: &GridState,
        accumulator: &mut Branch,
        index: usize,
    ) -> Fallible<()> {
        if index >= accumulator.chain.len() {
            return Ok(());
        }

        let group_nodes = group_directions(graph, accumulator, index)?
            .iter()
            .map(|&dir| group_node_tree(graph, accumulator.chain[index].pos, dir))
            .collect::<Fallible<Vec<GroupNode>>>()?;

        let groups = convert_nodes(group_nodes)?;
        accumulator.groups.push(groups);

        accumulate_groups(graph, accumulator, index + 1)
    }

    accumulate_groups(graph, &mut branch, 0usize)?;
    Ok(branch)
}

pub(crate) fn debug_branches(graph: &GridState) -> Fallible<Branch> {
    let all_chains = chain::get_all_chains(graph)?;
    let chain = chain::longest_chain(all_chains)?;
    link_groups(graph, chain)
}

/// Converts and combines the given `group_nodes` into [`Substituent`]s.
///
/// ## Errors
///
/// If a given [`GroupNode`] is not recognized, [`UnrecognizedGroup`] will be returned.
fn convert_nodes(group_nodes: Vec<GroupNode>) -> Fallible<Vec<Substituent>> {
    let groups = group_nodes
        .into_iter()
        .map(identify_single_bond_group)
        .collect::<Fallible<Vec<Group>>>()?;

    let new_groups = group_patterns(groups);

    Ok(new_groups)
}

/// Returns the [`Group`] corresponding to the structure of the given `node`.
///
/// ## Errors
///
/// Returns [`UnrecognizedGroup`] if the structure is not valid.
fn identify_single_bond_group(node: GroupNode) -> Fallible<Group> {
    let string = node.to_string();
    let out = match string.as_str() {
        "1O(1H)" => Hydroxyl,
        "2O" => Carbonyl,
        _ => return Err(UnrecognizedGroup),
    };

    Ok(out)
}

/// Recognizes compound groups from the given `groups` and converts all into [`Substituent`]s.
fn group_patterns(mut groups: Vec<Group>) -> Vec<Substituent> {
    let mut out = vec![];

    loop {
        if groups.contains(&Carbonyl) && groups.contains(&Hydroxyl) {
            groups.retain(|it| it != &Carbonyl && it != &Hydroxyl);
            out.push(Substituent::Group(Carboxyl));
            continue;
        }
        break;
    }

    let mut rest = groups
        .into_iter()
        .map(Substituent::Group)
        .collect::<Vec<Substituent>>();
    out.append(&mut rest);

    out
}

/// Constructs a [`GroupNode`] from the given `pos` in the given `direction`.
///
/// ## Panics
///
/// If the given `pos` does not point to a valid [`Cell::Atom`], or [`Direction::None`] is passed,
/// this function will panic.
///
/// ## Errors
///
/// If any invalid structures were found in the bonded group, an [`InvalidGraphError`] will be
/// returned.
pub(crate) fn group_node_tree(
    graph: &GridState,
    pos: Vec2,
    direction: Direction,
) -> Fallible<GroupNode> {
    let ptr = Pointer::new(graph, pos);
    let bond = ptr.bond_order(direction).unwrap();
    let atom = ptr.traverse_bond(direction)?;
    let mut next = vec![];

    for direction in next_directions(graph, atom.pos, pos)? {
        next.push(group_node_tree(graph, atom.pos, direction)?)
    }

    Ok(GroupNode {
        bond,
        atom: atom.element,
        next,
    })
}

/// Returns a [`Vec`] of [`Direction`] from the [`Atom`] at the given `pos` to bonded atoms
/// not including the one at the `previous_pos`.
///
/// ## Panics
///
/// This function assumes that the `previous_pos` is orthogonal to the given `pos`, that `pos`
/// points to a valid [`Cell::Atom`], and that there are no dangling bonds. If any of these
/// contracts are broken, this function will panic.
fn next_directions(graph: &GridState, pos: Vec2, previous_pos: Vec2) -> Fallible<Vec<Direction>> {
    let ptr = Pointer::new(graph, pos);
    let out = ptr
        .bonded()?
        .into_iter()
        .map(|atom| {
            Direction::from_points(pos, atom.pos)
                .map_err(|_| Other("An unexpected error occurred (groups/next_directions)."))
        })
        .filter(|dir| {
            if let Ok(it) = dir {
                it != &Direction::from_points(pos, previous_pos).unwrap()
            } else {
                true
            }
        })
        .collect::<Fallible<Vec<Direction>>>()?;

    Ok(out)
}

/// Returns a [`Vec`] of [`Direction`]s from the [`Atom`] at the given `index` to functional
/// groups.
///
/// ## Errors
///
/// If one of the bonds to the current cell is found to be dangling, an
/// [`IncompleteBond`] will be returned.
fn group_directions(
    graph: &GridState,
    accumulator: &Branch,
    index: usize,
) -> Fallible<Vec<Direction>> {
    let ptr = Pointer::new(graph, accumulator.chain[index].pos);
    let directions = ptr
        .connected_directions()
        .into_iter()
        .filter(|&direction| {
            let first_element = ptr.traverse_bond(direction).unwrap().element;
            let single_bond = matches!(ptr.bond_order(direction).unwrap(), BondOrder::Single);
            let hydrocarbon =
                matches!(first_element, Element::C) || matches!(first_element, Element::H);
            !single_bond || !hydrocarbon
        })
        .collect::<Vec<Direction>>();

    Ok(directions)
}

pub(crate) type Fallible<T> = Result<T, InvalidGraphError>;

/// The group of invalid structures that can appear on the [`GridState`].
#[derive(Error, Debug, PartialEq)]
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
    #[error("Unrecognized group.")]
    UnrecognizedGroup,
    #[error("{}", .0)]
    Other(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::{Double, Single};
    use crate::molecule::Element::{C, H, O};
    use crate::molecule::Group::{Bromo, Ether};
    use crate::test_utils::GW::{A, B};

    #[test]
    fn link_groups_recognizes_groups() {
        let graph = graph_with!(6, 3,
            [0, 1; A(H)],
            [1, 0; A(H)], [1, 1; A(C)], [1, 2; A(H)],
            [2, 1; B(Single)],
            [3, 0; A(H)], [3, 1; A(C)], [3, 2; A(H)],
            [4, 1; A(O)],
            [5, 1; A(H)]
        );
        let chain = vec![
            Atom {
                element: C,
                pos: Vec2::xy(1, 1),
            },
            Atom {
                element: C,
                pos: Vec2::xy(3, 1),
            },
        ];
        let branch = link_groups(&graph, chain.clone()).unwrap();
        let expected = Branch {
            chain,
            groups: vec![vec![], vec![Substituent::Group(Hydroxyl)]],
        };

        assert_eq!(branch, expected);
    }

    #[test]
    fn identify_single_bond_group_converts_group_node() {
        let hydroxyl_node = GroupNode {
            bond: Single,
            atom: O,
            next: vec![GroupNode {
                bond: Single,
                atom: H,
                next: vec![],
            }],
        };
        let carbonyl_node = GroupNode {
            bond: Double,
            atom: O,
            next: vec![],
        };
        let hydroxyl = identify_single_bond_group(hydroxyl_node).unwrap();
        let carbonyl = identify_single_bond_group(carbonyl_node).unwrap();

        assert_eq!(hydroxyl, Hydroxyl);
        assert_eq!(carbonyl, Carbonyl);
    }

    #[test]
    fn group_patterns_retains_groups() {
        let groups = vec![Carbonyl, Bromo, Ether];
        let expected = groups
            .iter()
            .map(|group| Substituent::Group(group.to_owned()))
            .collect::<Vec<Substituent>>();

        assert_eq!(group_patterns(groups), expected);
    }

    #[test]
    fn group_patterns_combines_to_carboxyl() {
        let groups = vec![Carbonyl, Hydroxyl];

        assert_eq!(
            group_patterns(groups)[0],
            vec![Substituent::Group(Carboxyl)][0]
        );
    }

    #[test]
    fn graph_node_tree() {
        let graph = graph_with!(1, 5,
            [0, 0; A(C)],
            [0, 1; B(Single)],
            [0, 2; A(O)],
            [0, 3; B(Single)],
            [0, 4; A(H)]
        );
        let a = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up).unwrap();
        let b = GroupNode {
            bond: Single,
            atom: O,
            next: vec![GroupNode {
                bond: Single,
                atom: H,
                next: vec![],
            }],
        };

        assert_eq!(a, b);
    }

    #[test]
    fn graph_node_tree_recognizes_implicit_bond() {
        let graph = graph_with!(1, 3,
            [0, 0; A(C)],
            [0, 1; A(O)],
            [0, 2; A(H)]
        );
        let a = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up).unwrap();
        let b = GroupNode {
            bond: Single,
            atom: O,
            next: vec![GroupNode {
                bond: Single,
                atom: H,
                next: vec![],
            }],
        };

        assert_eq!(a, b);
    }

    #[test]
    fn next_directions_omits_previous() {
        let graph = graph_with!(3, 3,
            [0, 1; A(C)],
            [1, 0; A(H)], [1, 1; A(C)],
            [2, 1; A(C)]
        );
        let directions = next_directions(&graph, Vec2::xy(1, 1), Vec2::xy(0, 1)).unwrap();
        let expected = vec![Direction::Down, Direction::Right];

        assert_eq!(directions, expected);
    }

    #[test]
    fn next_directions_behaves_at_boundaries() {
        let graph = graph_with!(3, 2,
            [0, 1; A(C)],
            [1, 0; A(H)],
            [1, 1; A(C)],
            [2, 1; A(C)]
        );
        let directions = next_directions(&graph, Vec2::xy(1, 1), Vec2::xy(0, 1)).unwrap();
        let expected = vec![Direction::Down, Direction::Right];

        assert_eq!(directions, expected);
    }

    #[test]
    fn group_directions_recognizes_groups() {
        let graph = graph_with!(5, 5,
            [0, 2; A(C)],
            [1, 2; B(Single)],
            [2, 0; A(O)], [2, 1; B(Double)], [2, 2; A(C)], [2, 3; A(O)], [2, 4; A(H)]
        );
        let branch = Branch {
            chain: vec![Atom {
                element: C,
                pos: Vec2::xy(2, 2),
            }],
            groups: vec![],
        };
        let directions = group_directions(&graph, &branch, 0usize).unwrap();

        assert_eq!(directions, vec![Direction::Up, Direction::Down]);
    }
}
