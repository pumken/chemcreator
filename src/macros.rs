//! # Macros
//!
//! Not to be confused with Rust's `macro_rules!` declarations, the `macros` module contains
//! common actions that should be automatically performed for the user when they make an input.

use std::ops::{Index, IndexMut};
use crate::molecule::{Cell, ComponentType};
use crate::molecule::Element::{C, H};
use crate::pointer::Pointer;
use crate::spatial::{EnumAll, GridState};
use ruscii::spatial::{Direction, Vec2};

#[derive(Clone, Debug, PartialEq)]
struct CellBlock<'a> {
    pub cells: Vec<CellRow<'a>>,
    pub direction: Direction,
    pub graph: &'a GridState,
}

impl<'a> CellBlock<'a> {
    pub fn new(graph: &GridState) -> CellBlock {
        CellBlock {
            cells: vec![],
            direction: Direction::Up,
            graph,
        }
    }
}

impl<'a> Index<usize> for CellBlock<'a> {
    type Output = CellRow<'a>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index]
    }
}

impl<'a> IndexMut<usize> for CellBlock<'a> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.cells[index]
    }
}

#[derive(Clone, Debug, PartialEq)]
struct CellRow<'a> {
    cells: Vec<Vec2>,
    graph: &'a GridState,
}

impl<'a> CellRow<'a> {
    pub fn new(graph: &'a GridState, cells: Vec<Vec2>) -> CellRow {
        CellRow { cells, graph }
    }

    pub fn borrow(&self, index: usize) -> Result<&'a Cell, String> {
        self.graph.get(self.cells[index])
    }
}

macro_rules! block {
    ($graph:expr, $([$(($xi:expr, $yi:expr)),*],)*) => {
        {
            let mut block = CellBlock::new($graph);
            $(
            {
                let temp = vec![$(Vec2::xy($xi, $yi),)*];
                block.cells.push(CellRow::new($graph, temp));
            }
            )*
            block
        }
    };
}

pub fn invoke_macro(graph: &mut GridState, new: ComponentType, previous: ComponentType) {
    match new {
        ComponentType::Element(C) => hydrogen_extension(graph),
        ComponentType::None => {} // just to appease Clippy
        _ => {}
    }
}

pub fn hydrogen_extension(graph: &mut GridState) {
    let bond_count = Pointer::new(graph, graph.cursor).bond_count().unwrap();
    let mut bonds_needed = graph.current_cell().unwrap().unwrap_atom().element.bond_number() - bond_count;

    for direction in Direction::all() {
        if bonds_needed <= 0 {
            return;
        }
        if hydrogen_arm(graph, direction) {
            bonds_needed -= 1;
        }
    }
}

fn hydrogen_arm(graph: &mut GridState, direction: Direction) -> bool {
    let selection = block!(graph,
        [(-1, 0), (0, 1), (1, 0)],
    );
    let mut ptr = Pointer::new(graph, graph.cursor);
    let first_neighbor = ptr.move_ptr(direction) && !ptr.borrow().unwrap().is_empty();
    let pos = ptr.pos;

    if first_neighbor {
        graph.put(pos, ComponentType::Element(H));
        true
    } else {
        false
    }
}
