//! # Groups
//!
//! The `groups` module provides functionality for identifying functional groups on a branch.

use ruscii::spatial::{Direction, Vec2};
use thiserror::Error;
use crate::groups::InvalidGraphError::{Other, UnrecognizedGroup};
use crate::spatial::{FromVec2, GridState, ToVec2};
use crate::molecule::{Atom, BondOrder, Branch, Cell, Element, Group, GroupNode, Substituent};
use crate::molecule::Group::{Carbonyl, Carboxyl, Hydroxyl};
use crate::{chain, validation};
use crate::pointer::Pointer;

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
    let branch = link_groups(graph, chain)?;
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
    match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    }
}

pub(crate) fn link_groups(graph: &GridState, chain: Vec<Atom>) -> Fallible<Branch> {
    let mut branch = Branch::new(chain);

    accumulate_groups(graph, &mut branch, 0usize)?;
    Ok(branch)
}

pub(crate) fn debug_branches(graph: &GridState) -> Fallible<Branch> {
    let all_chains = chain::get_all_chains(graph)?;
    let chain = chain::longest_chain(all_chains)?;
    link_groups(graph, chain)
}

fn accumulate_groups(
    graph: &GridState,
    accumulator: &mut Branch,
    index: usize,
) -> Fallible<()> {
    let group_nodes = group_directions(graph, accumulator, index)?
        .iter()
        .map(|&dir| group_node_tree(graph, accumulator.chain[index].pos, dir))
        .collect::<Fallible<Vec<GroupNode>>>()?;

    Ok(())
}

fn convert_nodes(group_nodes: Vec<GroupNode>) -> Fallible<Vec<Substituent>> {
    let mut out = vec![];
    let mut groups = vec![];

    for node in group_nodes {
        groups.push(identify_single_bond_group(node)?)
    }

    let new_groups = group_patterns(groups);

    Ok(out)
}

/// Returns the [`Group`] corresponding to the structure of the given `node`.
///
/// ## Errors
///
/// Returns [`UnrecognizedGroup`] if the structure is not valid.
fn identify_single_bond_group(node: GroupNode) -> Fallible<Group> {
    let string = node.to_string();
    let id = string.as_str();

    let out = match id {
        "1O(1H)" => Hydroxyl,
        "2O" => Carbonyl,
        _ => return Err(UnrecognizedGroup)
    };
    Ok(out)
}

fn group_patterns(mut groups: Vec<Group>) -> Vec<Substituent> {
    let mut out = vec![];

    while !groups.is_empty() {
        if groups.contains(&Carbonyl) && groups.contains(&Hydroxyl) {
            groups.retain(|it| it != &Carbonyl || it != &Hydroxyl);
            out.push(Substituent::Group(Carboxyl));
            continue
        }
        out.push(Substituent::Group(groups[0].clone()));
        groups.remove(0);
    }
    out
}

pub(crate) fn group_node_tree(graph: &GridState, pos: Vec2, direction: Direction) -> Fallible<GroupNode> {
    let ptr = Pointer { graph, pos };
    let bond = ptr.bond_order(direction).unwrap();
    let atom = ptr.traverse_bond(direction).unwrap();
    let mut next = vec![];

    for direction in next_directions(graph, atom.pos, pos)? {
        next.push(group_node_tree(graph, atom.pos, direction)?)
    }

    Ok(GroupNode { bond, atom: atom.element, next })
}

/// Returns a [`Vec`] of [`Direction`] from the [`Atom`] at the given `pos` to bonded atoms
/// not including the one at the `previous_pos`.
///
/// ## Panics
///
/// This function assumes that the `previous_pos` is orthogonal to the given `pos`, that `pos`
/// points to a valid [`Cell::Atom`], and that there are no dangling bonds. If any of these
/// contracts are broken, this function will panic.
fn next_directions(
    graph: &GridState,
    pos: Vec2,
    previous_pos: Vec2,
) -> Fallible<Vec<Direction>> {
    let ptr = Pointer { graph, pos };
    let bonded = ptr.bonded()?;
    let mut out = vec![];

    for atom in bonded {
        let result = match Direction::from_points(pos, atom.pos) {
            Ok(it) => it,
            Err(_) => return Err(Other("An unexpected error occurred (G151)."))
        };
        out.push(result)
    }
    out.retain(|&dir| dir != Direction::from_points(pos, previous_pos).unwrap());

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
    accumulator: &mut Branch,
    index: usize,
) -> Fallible<Vec<Direction>> {
    let ptr = Pointer { graph, pos: accumulator.chain[index].pos };
    let cells = ptr.connected();
    let mut out = vec![];

    for cell in cells {
        let traversal_ptr = Pointer { graph, pos: cell.pos() };
        let direction = Direction::from_points(ptr.pos, traversal_ptr.pos)
            .expect("Connected cells should be orthogonal");

        if !matches!(graph.get(ptr.pos + direction.to_vec2()), Ok(Cell::Bond(_))) {
            continue
        }

        match traversal_ptr.bond_order(direction) {
            Some(BondOrder::Double) | Some(BondOrder::Triple) => {
                out.push(direction);
                continue
            }
            _ => {}
        }

        match traversal_ptr.traverse_bond(direction)?.element {
            Element::C | Element::H => {}
            _ => out.push(direction),
        }
    }

    Ok(out)
}

pub(crate) type Fallible<T> = Result<T, InvalidGraphError>;

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
    use crate::molecule::Element::{C, H, O};
    use crate::molecule::BondOrder::Single;
    use crate::graph_with;
    use crate::test_utils::GW::{A, B};
    use super::*;

    #[test]
    fn graph_node_tree() {
        let graph = graph_with!(1, 5,
            [0, 0; A(C)],
            [0, 1; B(Single)],
            [0, 2; A(O)],
            [0, 3; B(Single)],
            [0, 4; A(H)]
        );
        let a = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up)
            .unwrap();
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
        let a = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up)
            .unwrap();
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
        let directions = next_directions(&graph, Vec2::xy(1, 1), Vec2::xy(0, 1))
            .unwrap();
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
        let directions = next_directions(&graph, Vec2::xy(1, 1), Vec2::xy(0, 1))
            .unwrap();
        let expected = vec![Direction::Down, Direction::Right];

        assert_eq!(directions, expected);
    }
}
