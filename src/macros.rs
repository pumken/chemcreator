//! # Macros
//!
//! Not to be confused with Rust's `macro_rules!` declarations, the `macros` module contains
//! common actions that should be automatically performed for the user when they make an input.

use std::ops::{Index, IndexMut};
use crate::molecule::{Cell, ComponentType};
use crate::molecule::Element::{C, H, O};
use crate::pointer::Pointer;
use crate::spatial::{EnumAll, GridState};
use ruscii::spatial::{Direction, Vec2};
use crate::molecule::BondOrder::{Double, Single};

#[derive(Clone, Debug, PartialEq)]
struct CellBlock<'a> {
    pub cells: Vec<Vec<Vec2>>,
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

    pub fn borrow(&self, group: usize, index: usize) -> Result<&'a Cell, String> {
        let relative = self.cells[group][index];
        let rotated = match self.direction {
            Direction::Up => relative,
            Direction::Down => -relative,
            Direction::Right => Vec2::xy(relative.y, -relative.x),
            Direction::Left => Vec2::xy(-relative.y, relative.x),
            Direction::None => panic!("Direction::None passed to borrow")
        };
        let absolute = rotated + self.graph.cursor;

        self.graph.get(absolute)
    }
}

macro_rules! block {
    ($graph:expr, $([$(($xi:expr, $yi:expr)),*],)*) => {
        {
            let mut block = CellBlock::new($graph);
            $(
            {
                let temp = vec![$(Vec2::xy($xi, $yi),)*];
                block.cells.push(temp);
            }
            )*
            block
        }
    };
}

pub fn invoke_macro(graph: &mut GridState, new: ComponentType, previous: ComponentType) {
    match new {
        ComponentType::Element(C) => {
            for direction in Direction::all() {
                let mut block = block!(graph,
                    [(0, 1), (0, 2)],
                );
                block.direction = direction;

                let first = match block.borrow(0, 0) {
                    Ok(it) => it,
                    Err(_) => continue,
                };
                let first_pos = first.pos();
                let second = match block.borrow(0, 1) {
                    Ok(it) => it,
                    Err(_) => continue,
                };

                let condition = second.is_atom() && second.unwrap_atom().element == C;

                if condition {
                    graph.put(first_pos, ComponentType::Order(Single));
                }
            }
            hydrogen_extension(graph)
        }
        ComponentType::Element(O) => {
            for direction in Direction::all() {
                let mut block = block!(graph,
                    [(0, 1), (0, 2)],
                );
                block.direction = direction;

                let first = match block.borrow(0, 0) {
                    Ok(it) => it,
                    Err(_) => continue,
                };
                let first_pos = first.pos();
                let second = match block.borrow(0, 1) {
                    Ok(it) => it,
                    Err(_) => continue,
                };

                let condition = second.is_atom() && second.unwrap_atom().element == C;

                if condition {
                    graph.put(first_pos, ComponentType::Order(Double));
                }
            }
        }
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
