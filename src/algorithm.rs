use ruscii::spatial::Vec2;
use thiserror::Error;
use crate::grid::{GridState, Pointer};
use crate::molecule::{Cell, Element, Group};

enum Substituent {
    Branch(Branch),
    Group(Group),
}

struct Branch {
    cells: Vec<Cell>,
    chain: Vec<Substituent>,
}

pub fn name_molecule(graph: &GridState) -> Result<String, InvalidGraphError> {
    let chain = link_in_chain(graph);
    // let group_indexed_chain = link_groups();
    // group_indexed_chain.check_chain_index()
    // group_indexed_chain.name();
    return match graph.simple_counter() {
        Ok(it) => Ok(it),
        Err(_) => Err(InvalidGraphError::UnsupportedGroups)
    }
}

fn link_in_chain(graph: &GridState) {
    let endpoints = endpoint_carbons(graph);
}

pub(crate) fn endpoint_carbons(graph: &GridState) -> Result<Vec<&Cell>, InvalidGraphError> {
    let all_carbons = graph.find_all(|cell| {
        match cell {
            Cell::Atom(atom) => match atom.element {
                Element::C => true,
                _ => false
            }
            _ => false
        }
    });
    let mut out = vec![];

    for carbon in all_carbons {
        let ptr = Pointer::new(carbon, graph);
        if ptr.bonded_carbons()? <= 1 {
            out.push(carbon);
        }
    }
    Ok(out)
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
}