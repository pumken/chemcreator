//! # Groups
//!
//! The `groups` module provides functionality for identifying functional groups on a branch.

use crate::chain::{endpoint_head_chains, get_all_chains, parent_chain};
use crate::compound;
use crate::groups::InvalidGraphError::{Other, UnrecognizedGroup};
use crate::molecule::Group::{
    self, AcidHalide, Aldehyde, Alkene, Alkyne, Amide, Amine, Bromo, Carbonyl, Carboxyl, Chloro,
    Fluoro, Hydrogen, Hydroxyl, Iodo, Nitrile,
};
use crate::molecule::Halogen::{Bromine, Chlorine, Fluorine, Iodine};
use crate::molecule::{Atom, BondOrder, Branch, Element, GroupNode, Substituent};
use crate::pointer::Pointer;
use crate::spatial::{GridState, InvertDirection};
use ruscii::spatial::{Direction, Vec2};
use thiserror::Error;

/// Generates a [`Branch`] from the given `chain` containing all functional groups attached
/// to it, except for the group starting at the given `parent`.
pub(crate) fn link_groups(
    graph: &GridState,
    chain: Vec<Atom>,
    parent: Option<Atom>,
) -> Fallible<Branch> {
    let mut branch = Branch::new(chain);

    fn accumulate_groups(
        graph: &GridState,
        accumulator: &mut Branch,
        index: usize,
        parent: Option<Atom>,
    ) -> Fallible<()> {
        if index >= accumulator.chain.len() {
            return Ok(());
        }

        let directions = group_directions(graph, accumulator, index, parent.clone())?;

        let group_nodes = directions
            .0
            .into_iter()
            .map(|dir| group_node_tree(graph, accumulator.chain[index].pos, dir))
            .collect::<Fallible<Vec<GroupNode>>>()?;

        let mut chain_nodes: Vec<Substituent> = directions
            .1
            .into_iter()
            .map(|dir| branch_node_tree(graph, accumulator.chain[index].pos, dir))
            .collect::<Fallible<Vec<Vec<Atom>>>>()?
            .into_iter()
            .map(|chain| {
                link_groups(
                    graph,
                    chain,
                    Some(Atom {
                        element: Element::C,
                        pos: accumulator.chain[index].pos,
                    }),
                )
            })
            .collect::<Fallible<Vec<Branch>>>()?
            .into_iter()
            .map(Substituent::Branch)
            .collect();

        let mut groups = convert_nodes(group_nodes)?;
        groups.append(&mut chain_nodes);
        accumulator.groups.push(groups);

        accumulate_groups(graph, accumulator, index + 1, parent)
    }

    accumulate_groups(graph, &mut branch, 0usize, parent)?;
    Ok(branch)
}

pub(crate) fn debug_branches(graph: &GridState) -> Fallible<Branch> {
    let all_chains = get_all_chains(graph)?;
    let chain = parent_chain(graph, all_chains, None)?;
    link_groups(graph, chain, None)
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
        "1H" => Hydrogen,
        "1B" => Bromo,
        "1L" => Chloro,
        "1F" => Fluoro,
        "1I" => Iodo,
        "1O(1H)" => Hydroxyl,
        "1N(1H)(1H)" => Amine,
        "2O" => Carbonyl,
        "2C" => Alkene,
        "3C" => Alkyne,
        "3N" => Nitrile,
        _ => return Err(UnrecognizedGroup),
    };

    Ok(out)
}

/// Recognizes compound groups from the given `groups` and converts all into [`Substituent`]s.
fn group_patterns(mut groups: Vec<Group>) -> Vec<Substituent> {
    let mut out = vec![];

    loop {
        compound!(groups, out,
            [Carbonyl, Hydroxyl => Carboxyl],
            [Carbonyl, Fluoro => AcidHalide(Fluorine)],
            [Carbonyl, Chloro => AcidHalide(Chlorine)],
            [Carbonyl, Bromo => AcidHalide(Bromine)],
            [Carbonyl, Iodo => AcidHalide(Iodine)],
            [Carbonyl, Amine => Amide],
            [Carbonyl, Hydrogen => Aldehyde],
        );
        break;
    }
    groups.retain(|it| it != &Hydrogen);

    let mut rest: Vec<Substituent> = groups.into_iter().map(Substituent::Group).collect();
    out.append(&mut rest);

    out
}

/// Takes the given `groups`, matches them against the cases given and returns the resultant
/// groups to `out`.
#[macro_export]
macro_rules! compound {
    ($groups:expr, $out:expr, $([$first:expr, $second:expr => $comp:expr],)*) => {
        $(
        if $groups.contains(&$first) && $groups.contains(&$second) {
            if let Some(index) = $groups.iter().position(|&x| x == $first) {
                $groups.remove(index);
            }
            if let Some(index) = $groups.iter().position(|&x| x == $second) {
                $groups.remove(index);
            }
            $out.push(Substituent::Group($comp));
            continue;
        }
        )*
    };
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

    if atom.element != Element::C {
        for direction in next_directions(graph, atom.pos, pos)? {
            next.push(group_node_tree(graph, atom.pos, direction)?)
        }
    }

    Ok(GroupNode {
        bond,
        atom: atom.element,
        next,
    })
}

pub(crate) fn branch_node_tree(
    graph: &GridState,
    pos: Vec2,
    direction: Direction,
) -> Fallible<Vec<Atom>> {
    let atom = Pointer::new(graph, pos).traverse_bond(direction)?;
    let chains = endpoint_head_chains(atom, graph, Some(pos))?;
    parent_chain(
        graph,
        chains,
        Some(
            graph
                .get(pos)
                .map_err(|_| Other("An unexpected error occurred.".to_string()))?
                .unwrap_atom(),
        ),
    )
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
            Direction::try_from(atom.pos - pos).inv().map_err(|_| {
                Other("An unexpected error occurred (groups/next_directions).".to_string())
            })
        })
        .filter(|dir| {
            if let Ok(it) = dir {
                it != &Direction::try_from(previous_pos - pos).inv().unwrap()
            } else {
                true
            }
        })
        .collect::<Fallible<Vec<Direction>>>()?;

    Ok(out)
}

/// Returns a [`Vec`] of [`Direction`]s from the [`Atom`] at the given `index` to functional
/// groups and side chains, respectively.
///
/// ## Errors
///
/// If one of the bonds to the current cell is found to be dangling, an
/// [`IncompleteBond`] will be returned.
fn group_directions(
    graph: &GridState,
    accumulator: &Branch,
    index: usize,
    parent: Option<Atom>,
) -> Fallible<(Vec<Direction>, Vec<Direction>)> {
    let ptr = Pointer::new(graph, accumulator.chain[index].pos);
    let directions = ptr
        .connected_directions()
        .into_iter()
        .filter(|&direction| {
            let opposite_atom = ptr.traverse_bond(direction).unwrap();
            let single_bond = matches!(ptr.bond_order(direction).unwrap(), BondOrder::Single);
            let in_chain = accumulator.chain.contains(&opposite_atom);
            let parent = Some(opposite_atom) == parent;
            !(parent || in_chain && single_bond)
        })
        .filter(|&direction| {
            if index > 0 {
                direction
                    != Direction::try_from(
                        accumulator.chain[index - 1].pos - accumulator.chain[index].pos,
                    )
                    .inv()
                    .expect("Consecutive indexes should return orthogonal points.")
            } else {
                true
            }
        })
        .partition(|&direction| {
            let opposite_atom = ptr.traverse_bond(direction).unwrap();
            let carbon = matches!(opposite_atom.element, Element::C);
            let in_chain = accumulator.chain.contains(&opposite_atom);
            in_chain || !carbon
        });

    Ok(directions)
}

/// A type that is only available when the [`GridState`] is valid.
pub(crate) type Fallible<T> = Result<T, InvalidGraphError>;

/// The group of invalid structures that can appear on the [`GridState`].
#[derive(Error, Debug, PartialEq)]
pub enum InvalidGraphError {
    #[error("Molecule is not continuous.")]
    Discontinuity,
    #[error("Molecule is cyclic.")]
    Cycle,
    #[error("Bond at {0} has no inconsistent order.")]
    InconsistentBond(Vec2),
    #[error("Atom at {0} is missing bonds ({1} requires {}).", .1.bond_number())]
    UnfilledValence(Vec2, Element),
    #[error("Atom at {0} has too many bonds ({1} requires {}).", .1.bond_number())]
    OverfilledValence(Vec2, Element),
    #[error("Bond at {0} is incomplete.")]
    IncompleteBond(Vec2),
    #[error("Unrecognized group.")]
    UnrecognizedGroup,
    #[error("{0}")]
    Other(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::{Double, Single};
    use crate::molecule::Element::{Br, Cl, C, F, H, I, O};
    use crate::molecule::Group::Ether;
    use crate::spatial::EnumAll;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn link_groups_recognizes_groups() {
        let graph = graph_with!(6, 3,
            [0, 1; A(H)],
            [1, 0; A(H)], [1, 1; A(C)], [1, 2; A(Cl)],
            [2, 1; B(Single)],
            [3, 0; A(H)], [3, 1; A(C)], [3, 2; A(H)],
            [4, 1; A(O)],
            [5, 1; A(H)],
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
        let branch = link_groups(&graph, chain.clone(), None).unwrap();
        let expected = Branch {
            chain,
            groups: vec![
                vec![Substituent::Group(Chloro)],
                vec![Substituent::Group(Hydroxyl)],
            ],
            parent_alpha: None,
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
        let groups = vec![Carbonyl, Ether];
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
    fn group_node_tree_parses_structure() {
        let graph = graph_with!(1, 5,
            [0, 0; A(C)],
            [0, 1; B(Single)],
            [0, 2; A(O)],
            [0, 3; B(Single)],
            [0, 4; A(H)],
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
    fn group_node_tree_recognizes_implicit_bond() {
        let graph = graph_with!(1, 3,
            [0, 0; A(C)],
            [0, 1; A(O)],
            [0, 2; A(H)],
        );
        let node = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up).unwrap();
        let expected = GroupNode {
            bond: Single,
            atom: O,
            next: vec![GroupNode {
                bond: Single,
                atom: H,
                next: vec![],
            }],
        };

        assert_eq!(node, expected);
    }

    #[test]
    fn group_node_tree_recognizes_carbon_group() {
        let graph = graph_with!(3, 3,
            [0, 0; A(H)], [0, 1; A(C)], [0, 2; A(H)],
            [1, 1; B(Double)],
            [2, 0; A(H)], [2, 1; A(C)], [2, 2; A(H)],
        );
        let node = group_node_tree(&graph, Vec2::xy(0, 1), Direction::Right).unwrap();
        let expected = GroupNode {
            bond: Double,
            atom: C,
            next: vec![],
        };

        assert_eq!(node, expected);
    }

    #[test]
    fn next_directions_omits_previous() {
        let graph = graph_with!(3, 3,
            [0, 1; A(C)],
            [1, 0; A(H)], [1, 1; A(C)],
            [2, 1; A(C)],
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
            [2, 1; A(C)],
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
            [2, 0; A(C)], [2, 1; B(Double)], [2, 2; A(C)], [2, 3; A(O)], [2, 4; A(H)],
        );
        let branch = Branch {
            chain: vec![
                Atom {
                    element: C,
                    pos: Vec2::xy(0, 2),
                },
                Atom {
                    element: C,
                    pos: Vec2::xy(2, 2),
                },
            ],
            groups: vec![],
            parent_alpha: None,
        };
        let directions = group_directions(&graph, &branch, 1usize, None).unwrap();

        assert_eq!(directions.0, vec![Direction::Up]);
        assert_eq!(directions.1, vec![Direction::Down]);
    }

    #[test]
    fn group_directions_recognizes_halogens() {
        let graph = graph_with!(3, 3,
            [0, 1; A(Cl)],
            [1, 0; A(Br)],
            [1, 1; A(C)],
            [1, 2; A(I)],
            [2, 1; A(F)],
        );
        let branch = Branch {
            chain: vec![Atom {
                element: C,
                pos: Vec2::xy(1, 1),
            }],
            groups: vec![],
            parent_alpha: None,
        };
        let directions = group_directions(&graph, &branch, 0usize, None).unwrap();

        assert_eq!(directions.0, Direction::all());
    }
}
