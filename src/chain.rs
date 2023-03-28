//! # Chain
//!
//! The `chain` module contains functions that allow for identifying the parent chain of an
//! organic molecule.

use crate::groups::InvalidGraphError::{Cycle, Other};
use crate::groups::{link_groups, Fallible};
use crate::molecule::Group::{Alkene, Alkyne};
use crate::molecule::{Atom, Branch, Cell, Element, Substituent};
use crate::nested_vec;
use crate::pointer::Pointer;
use crate::spatial::GridState;
use ruscii::spatial::Vec2;

/// Gets the parent chain of the given [`Vec`] of chains according to
/// [IUPAC rules](https://en.wikipedia.org/wiki/IUPAC_nomenclature_of_organic_chemistry) on
/// Wikipedia, accessed 26 Mar 2023. The parent chain is, in order of precedence, that which:
/// 1. has the greatest number of occurrences of the primary functional group.
/// 2. has the greatest number of multiple bonds.
/// 3. is the longest.
/// 3. has the greatest number of occurrences of prefix substituents.
/// 4. has the greatest number of single bonds.
pub fn parent_chain(
    graph: &GridState,
    chains: Vec<Vec<Atom>>,
    parent: Option<Atom>,
) -> Fallible<Vec<Atom>> {
    if chains.is_empty() {
        return Err(Other("No carbon chain found.".to_string()));
    }

    let branches = chains
        .into_iter()
        .map(|chain| link_groups(graph, chain, parent.clone()))
        .collect::<Fallible<Vec<Branch>>>()?;

    let mut unique_branches = vec![];

    for branch in branches.into_iter() {
        if !unique_branches.contains(&branch.reversed()) {
            unique_branches.push(branch);
        }
    }

    if unique_branches.is_empty() {
        panic!("unique_branches is empty!")
    }

    let primary_group = &unique_branches
        .iter()
        .flat_map(|it| it.groups.to_owned())
        .flatten()
        .filter_map(|it| {
            if let Substituent::Group(group) = it {
                if group.seniority().is_some() {
                    Some(group)
                } else {
                    None
                }
            } else {
                None
            }
        })
        .max_by_key(|it| it.seniority());

    let mut filtration = unique_branches.to_owned();

    // 1
    if let Some(it) = &primary_group {
        let most_primary_groups = chain_max_by(unique_branches, |branch| {
            branch
                .groups
                .iter()
                .flatten()
                .filter(|&subst| matches!(subst, Substituent::Group(group) if group == it))
                .count()
        });
        if let Ok(it) = most_primary_groups {
            return Ok(it);
        }
        filtration = most_primary_groups.unwrap_err();
    }

    // 2
    let most_multiple_bonds = chain_max_by(filtration, |branch| {
        branch
            .groups
            .iter()
            .flatten()
            .filter(|&subst| {
                matches!(
                    subst,
                    Substituent::Group(Alkene) | Substituent::Group(Alkyne)
                )
            })
            .count()
    });
    if let Ok(it) = most_multiple_bonds {
        return Ok(it);
    }
    let filtration = most_multiple_bonds.unwrap_err();

    // 3
    let longest = chain_max_by(filtration, |branch| branch.chain.len());
    if let Ok(it) = longest {
        return Ok(it);
    }
    let filtration = longest.unwrap_err();

    // 4
    let most_prefix_subst = chain_max_by(filtration, |branch| {
        branch
            .groups
            .iter()
            .flatten()
            .filter(|&subst| {
                if let Some(it) = primary_group {
                    matches!(subst, Substituent::Group(group) if group != it)
                        || matches!(subst, Substituent::Branch(_))
                } else {
                    true
                }
            })
            .count()
    });
    if let Ok(it) = most_prefix_subst {
        return Ok(it);
    }
    let filtration = most_prefix_subst.unwrap_err();

    // TODO 5

    Ok(filtration[0].chain.to_owned())
}

/// Returns the maximum of the given `branches` based on the values yielded by given `f` or filters
/// them.
///
/// ## Errors
///
/// If a maximum cannot be attained (e.g., if there are no branches), an [`Err`] containing all
/// original [`Branch`]es is returned. If there are multiple branches matching the maximum, an
/// [`Err`] containing matching [`Branch`]es is returned.
pub(crate) fn chain_max_by<F>(branches: Vec<Branch>, f: F) -> Result<Vec<Atom>, Vec<Branch>>
where
    F: Fn(&Branch) -> usize,
{
    let max: usize = branches.iter().map(&f).max().ok_or(branches.to_owned())?;

    let primary_by_max = branches
        .into_iter()
        .filter(|it| f(it) == max)
        .collect::<Vec<Branch>>();

    if primary_by_max.len() == 1 {
        Ok(primary_by_max[0].to_owned().chain)
    } else {
        Err(primary_by_max)
    }
}

/// Returns a [`Vec`] of chains between every endpoint on the graph, thereby yielding every
/// possible candidate for the primary chain.
///
/// ## Errors
///
/// Returns [`InvalidGraphError`] if any invalid structures are found while traversing the graph.
pub(crate) fn get_all_chains(graph: &GridState) -> Fallible<Vec<Vec<Atom>>> {
    let endpoints = endpoint_carbons(graph)?
        .iter()
        .map(|&cell| match cell {
            Cell::Atom(it) => it.to_owned(),
            _ => panic!("endpoint_carbons returned non-atom cell"),
        })
        .collect::<Vec<Atom>>();
    // Nested Vec hell
    let mut out: Vec<Vec<Atom>> = vec![];

    for endpoint in endpoints {
        out.extend(endpoint_head_chains(endpoint.to_owned(), graph, None)?);
    }
    Ok(out)
}

/// Takes a given `endpoint` and returns all continuous carbon chains starting at it.
///
/// ## Errors
///
/// Returns [`InvalidGraphError`] if any invalid structures are found while traversing the graph.
pub(crate) fn endpoint_head_chains(
    endpoint: Atom,
    graph: &GridState,
    previous_pos: Option<Vec2>,
) -> Fallible<Vec<Vec<Atom>>> {
    let mut accumulator = vec![vec![]];

    accumulate_carbons(endpoint.pos, previous_pos, 0usize, &mut accumulator, graph)?;

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
) -> Fallible<()> {
    let next_carbons = next_carbons(pos, previous_pos, graph)?;

    accumulator[branch_index].push(match graph.get(pos) {
        Ok(Cell::Atom(it)) => it.to_owned(),
        _ => panic!("Non-atom or invalid cell passed to accumulate_carbons"),
    });

    let new_branches = create_branches(
        accumulator,
        branch_index,
        match next_carbons.len() {
            0 => return Ok(()),
            it => it - 1,
        },
    );
    for (i, carbon) in next_carbons.iter().enumerate() {
        accumulate_carbons(carbon.pos, Some(pos), new_branches[i], accumulator, graph)?
    }

    Ok(())
}

/// Creates `count` clones of the branch at the given `branch_index` in the `accumulator`. The
/// returned [`Vec`] contains the indexes of the branch and its copies.
fn create_branches(
    accumulator: &mut Vec<Vec<Atom>>,
    branch_index: usize,
    count: usize,
) -> Vec<usize> {
    let mut out = vec![branch_index];

    for _ in 0usize..count {
        accumulator.push(accumulator[branch_index].clone());
        out.push(accumulator.len() - 1)
    }

    out
}

/// Gets the next bonded carbons. Does not include the cell at the `previous_pos` if there is one.
///
/// ## Errors
///
/// If one of the bonds to the current cell is found to be dangling, an
/// [`IncompleteBond`] will be returned.
fn next_carbons(pos: Vec2, previous_pos: Option<Vec2>, graph: &GridState) -> Fallible<Vec<Atom>> {
    let ptr = Pointer::new(graph, pos);
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
pub(crate) fn endpoint_carbons(graph: &GridState) -> Fallible<Vec<&Cell>> {
    let all_carbons = graph.find_all(|cell| match cell {
        Cell::Atom(atom) => matches!(atom.element, Element::C),
        _ => false,
    });
    let mut out = vec![];

    for carbon in all_carbons {
        let ptr = Pointer::new(graph, carbon.pos());
        if ptr.bonded_carbon_count()? <= 1 {
            out.push(carbon);
        }
    }
    Ok(out)
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
pub(crate) fn get_connected_cells(pos: Vec2, graph: &GridState) -> Fallible<Vec<Vec<bool>>> {
    if let Cell::None(_) = graph
        .get(pos)
        .expect("pos should be a valid point on the graph")
    {
        panic!(
            "Passed empty cell ({}, {}) to get_connected_cells",
            pos.x, pos.y
        )
    }

    let mut searched_points = nested_vec![graph.size.x; graph.size.y; false];

    fn accumulate_components(
        pos: Vec2,
        previous_pos: Option<Vec2>,
        accumulator: &mut Vec<Vec<bool>>,
        graph: &GridState,
    ) -> Fallible<()> {
        accumulator[pos.x as usize][pos.y as usize] = true;
        let ptr = Pointer::new(graph, pos);
        for cell in ptr.connected() {
            if let Some(it) = previous_pos {
                if cell.pos() == it {
                    continue;
                }
            }
            match cell {
                Cell::None(_) => {}
                _ => {
                    if !accumulator[cell.pos().x as usize][cell.pos().y as usize] {
                        accumulate_components(cell.pos(), Some(pos), accumulator, graph)?
                    } else {
                        return Err(Cycle);
                    }
                }
            }
        }
        Ok(())
    }

    accumulate_components(pos, None, &mut searched_points, graph)?;
    Ok(searched_points)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::Single;
    use crate::molecule::Element::{C, H};
    use crate::test_utils::GW::{A, B};

    // Also checks accumulate_carbons
    #[test]
    fn endpoint_head_chains_finds_all_atoms() {
        let graph = graph_with!(5, 3,
            [0, 0; A(C)],
            [1, 0; A(C)],
            [2, 0; A(C)], [2, 1; A(C)], [2, 2; A(C)],
            [3, 0; A(C)],
            [4, 0; A(C)],
        );
        let accumulator = endpoint_head_chains(
            Atom {
                element: C,
                pos: Vec2::zero(),
            },
            &graph,
            None,
        )
        .unwrap();

        assert_eq!(accumulator.len(), 2);
        assert_eq!(accumulator[0].len(), 5);
        assert_eq!(accumulator[1].len(), 5);
        assert_ne!(accumulator[0], accumulator[1]);
    }

    #[test]
    fn endpoint_head_chains_ignores_previous() {
        let graph = graph_with!(5, 3,
            [0, 0; A(C)],
            [1, 0; A(C)],
            [2, 0; A(C)], [2, 1; A(C)], [2, 2; A(C)],
            [3, 0; A(C)],
            [4, 0; A(C)],
        );
        let accumulator = endpoint_head_chains(
            Atom {
                element: C,
                pos: Vec2::xy(1, 0),
            },
            &graph,
            Some(Vec2::zero()),
        )
        .unwrap();

        assert_eq!(accumulator.len(), 2);
        assert_eq!(accumulator[0].len(), 4);
        assert_eq!(accumulator[1].len(), 4);
        assert_ne!(accumulator[0], accumulator[1]);
        assert!(!accumulator
            .into_iter()
            .flatten()
            .collect::<Vec<Atom>>()
            .contains({
                if let Cell::Atom(it) = graph.get(Vec2::zero()).unwrap() {
                    it
                } else {
                    panic!("")
                }
            }));
    }

    #[test]
    fn create_branches_copies_correctly() {
        let atom = Atom {
            element: C,
            pos: Vec2::zero(),
        };
        let mut accumulator = vec![vec![atom.clone()], vec![]];
        create_branches(&mut accumulator, 0, 2);

        assert_eq!(accumulator.len(), 4);
        assert_eq!(accumulator[1], vec![]);
        assert_eq!(accumulator[2], vec![atom.clone()]);
        assert_eq!(accumulator[3], vec![atom]);
    }

    #[test]
    fn next_carbons_omits_non_carbons() {
        let graph = graph_with!(3, 3,
            [0, 1; A(H)],
            [1, 0; A(H)],
            [1, 1; A(C)],
            [1, 2; A(C)],
            [2, 1; A(C)],
        );
        let atoms = next_carbons(Vec2::xy(1, 1), None, &graph).unwrap();
        let expected = vec![
            graph.get(Vec2::xy(1, 2)).unwrap(),
            graph.get(Vec2::xy(2, 1)).unwrap(),
        ]
        .iter()
        .map(|&cell| cell.unwrap_atom())
        .collect::<Vec<Atom>>();

        assert_eq!(atoms, expected);
    }

    #[test]
    fn next_carbons_omits_previous() {
        let graph = graph_with!(3, 3,
            [0, 1; A(H)],
            [1, 0; A(H)],
            [1, 1; A(C)],
            [1, 2; A(C)],
            [2, 1; A(C)],
        );
        let atom = next_carbons(Vec2::xy(1, 1), Some(Vec2::xy(1, 2)), &graph).unwrap();
        let expected = vec![graph.get(Vec2::xy(2, 1)).unwrap().unwrap_atom()];

        assert_eq!(atom, expected);
    }

    #[test]
    fn endpoint_carbons_returns_one_carbon_neighbor() {
        let graph = graph_with!(7, 3,
            [0, 1; A(H)],
            [1, 0; A(H)], [1, 1; A(C)], [1, 2; A(H)],
            [2, 1; B(Single)],
            [3, 0; A(H)], [3, 1; A(C)], [3, 2; A(H)],
            [4, 1; B(Single)],
            [5, 0; A(H)], [5, 1; A(C)], [5, 2; A(H)],
            [6, 1; A(H)],
        );
        let cells = endpoint_carbons(&graph).unwrap();
        let expected = vec![
            graph.get(Vec2::xy(1, 1)).unwrap(),
            graph.get(Vec2::xy(5, 1)).unwrap(),
        ];

        assert_eq!(cells, expected);
    }

    #[test]
    fn endpoint_carbons_returns_zero_carbon_neighbor() {
        let graph = graph_with!(3, 3,
            [0, 1; A(H)],
            [1, 0; A(H)], [1, 1; A(C)], [1, 2; A(H)],
            [2, 1; A(H)],
        );
        let cell = endpoint_carbons(&graph).unwrap();
        let expected = vec![graph.get(Vec2::xy(1, 1)).unwrap()];

        assert_eq!(cell, expected);
    }

    #[test]
    fn get_connected_cells_only_returns_connected() {
        let graph = graph_with!(3, 3,
            [0, 0; A(C)],
            [1, 0; A(C)],
            [2, 0; A(C)],
            [0, 2; A(C)],
        );
        let network = get_connected_cells(Vec2::xy(0, 0), &graph).unwrap();
        let expected = vec![
            vec![true, false, false],
            vec![true, false, false],
            vec![true, false, false],
        ];

        assert_eq!(network, expected);
    }
}
