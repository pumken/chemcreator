use thiserror::Error;
use crate::grid::{GridState, Pointer};
use crate::molecule::{Cell, Element, Group};

enum Substituent {
    Branch(Branch),
    Group(Group)
}

struct Branch {
    cells: Vec<Cell>,
    chain: Vec<Substituent>
}

fn name_molecule(graph: &GridState) -> String {
    let chain = link_in_chain(graph);
    // let group_indexed_chain = link_groups();
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
}

fn link_in_chain(graph: &GridState) {
    let endpoints = endpoint_carbons(graph);
}

fn endpoint_carbons(graph: &GridState) -> Vec<&Cell> {
    graph.find_all(|cell| {
        if cell.is_not_atom() { return false; }
        let ptr = Pointer::new(cell, graph);
        ptr.bonded()
            .expect("is_not_atom() check should not fail")
            .iter()
            .filter(|atom| if let Cell::Atom(it) = atom {
                it == Element::C
            } else {
                false
            })
            .collect()
    })
}

#[derive(Error, Debug)]
pub(crate) enum InvalidGraphError<'a> {
    #[error("Molecule is not continuous.")]
    Discontinuity,
    #[error("Molecule is cyclic.")]
    Cycle,
    #[error("Cell at ({}, {}) is missing bonds.", .0.pos().x, .0.pos().y)]
    UnfilledValence(&'a Cell),
    #[error("Cell at ({}, {}) has too many bonds.", .0.pos().x, .0.pos().y)]
    OverfilledValence(&'a Cell),
    #[error("Bond at ({}, {}) is incomplete.", .0.pos().x, .0.pos().y)]
    IncompleteBond(&'a Cell),
    #[error("This combination of groups is not supported.")]
    UnsupportedGroups,
}