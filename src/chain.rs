use ruscii::spatial::Vec2;
use crate::algorithm::Fallible;
use crate::algorithm::InvalidGraphError::{Cycle, Other};
use crate::molecule::{Atom, Cell, Element};
use crate::nested_vec;
use crate::pointer::Pointer;
use crate::spatial::GridState;

pub(crate) fn debug_chain(graph: &GridState) -> Fallible<Vec<Atom>> {
    let all_chains = get_all_chains(graph)?;
    longest_chain(all_chains)
}

/// Gets the longest of the given [`Vec`] of chains, assuming that it is non-empty.
///
/// ## Errors
///
/// If there are no `chains`, this function will return an `Err`.
pub(crate) fn longest_chain(chains: Vec<Vec<Atom>>) -> Fallible<Vec<Atom>> {
    Ok(match chains.iter()
        .max_by(|&a, &b| a.len().cmp(&b.len())) {
        None => return Err(Other("No carbon chain found.".to_string())),  // FIXME this is returned when a bond is placed at the edge
        Some(it) => it
    }.to_owned())
}

pub(crate) fn get_all_chains(graph: &GridState) -> Fallible<Vec<Vec<Atom>>> {
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
fn endpoint_head_chains(endpoint: Atom, graph: &GridState) -> Fallible<Vec<Vec<Atom>>> {
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
) -> Fallible<()> {
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
fn next_carbons(pos: Vec2, previous_pos: Option<Vec2>, graph: &GridState) -> Fallible<Vec<Atom>> {
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
pub(crate) fn endpoint_carbons(graph: &GridState) -> Fallible<Vec<&Cell>> {
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
    if let Cell::None(_) = graph.get(pos).expect("pos should be a valid point on the graph.") {
        panic!("Passed empty cell ({}, {}) to get_connected_cells", pos.x, pos.y)
    }

    let mut searched_points = nested_vec![graph.size.x; graph.size.y; false];

    fn accumulate_components(
        pos: Vec2,
        previous_pos: Option<Vec2>,
        accumulator: &mut Vec<Vec<bool>>,
        graph: &GridState,
    ) -> Fallible<()> {
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
