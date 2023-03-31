//! # Macros
//!
//! Not to be confused with Rust's `macro_rules!` declarations, the `macros` module contains
//! common actions that should be automatically performed for the user when they make an input.

use crate::molecule::BondOrder::{Double, Single};
use crate::molecule::Element::{C, H, O};
use crate::molecule::{BondOrder, BondOrientation, Cell, ComponentType};
use crate::pointer::Pointer;
use crate::spatial::{EnumAll, GridState, ToVec2};
use ruscii::spatial::{Direction, Vec2};
use std::cmp::Ordering;

use crate::molecule::BondOrientation::Horiz;

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
            Direction::None => panic!("Direction::None passed to borrow"),
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

pub fn invoke_macro(graph: &mut GridState, new: ComponentType, _previous: ComponentType) {
    match new {
        ComponentType::Element(C) => {
            let _ = carbon_extension(graph) || methane_creation(graph) || carbon_extension_alt(graph);
        }
        ComponentType::Element(O) => {
            let _ = carbonyl_extension(graph) || hydroxyl_extension(graph);
        }
        ComponentType::Order(_) => {
            cxc_bond_correction(graph)
        }
        _ => {}
    }
}

fn cxc_bond_correction(graph: &mut GridState) {
    let dirs = if graph.current_cell().unwrap().unwrap_bond().orient == Horiz {
        [Direction::Left, Direction::Right]
    } else {
        [Direction::Up, Direction::Down]
    };
    let ptr = Pointer::new(graph, graph.cursor);
    let first = match ptr.traverse_bond(dirs[0]) {
        Ok(it) => it,
        Err(_) => return,
    };
    let second = match ptr.traverse_bond(dirs[1]) {
        Ok(it) => it,
        Err(_) => return,
    };

    if first.element == C && second.element == C {
        hydrogen_correction(graph, first.pos);
        hydrogen_correction(graph, second.pos)
    }
}

fn carbon_extension(graph: &mut GridState) -> bool {
    let mut block = block!(graph, [(0, 1), (0, 2)],);

    for direction in Direction::all() {
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

        let atom_condition = second.is_atom() && second.unwrap_atom().element == C;
        let bond_condition = second.is_bond()
            && BondOrientation::from(direction) == second.unwrap_bond().orient
            && !first.is_atom();

        if atom_condition || bond_condition {
            graph.put(first_pos, ComponentType::Order(Single));
            hydrogen_correction(graph, graph.cursor);
            return true;
        }
    }

    false
}

fn carbon_extension_alt(graph: &mut GridState) -> bool {
    let mut block = block!(graph, [(0, -1), (0, 1)],);

    for direction in Direction::all() {
        block.direction = direction;

        let first = match block.borrow(0, 0) {
            Ok(it) => it,
            Err(_) => continue,
        };
        let second = match block.borrow(0, 1) {
            Ok(it) => it,
            Err(_) => continue,
        };
        let second_pos = second.pos();

        let condition = first.is_atom() && first.unwrap_atom().element == C &&
            (second.is_empty() || (second.is_atom() && second.unwrap_atom().element == C));

        if condition {
            graph.put(second_pos, ComponentType::Element(C));
            graph.put(graph.cursor, ComponentType::Order(Single));
            hydrogen_correction(graph, second_pos);
            return true;
        }
    }

    false
}

fn methane_creation(graph: &mut GridState) -> bool {
    let ptr = Pointer::new(graph, graph.cursor);

    if ptr.connected_directions().is_empty() {
        hydrogen_correction(graph, graph.cursor);
        return true;
    }

    false
}

fn hydroxyl_extension(graph: &mut GridState) -> bool {
    let mut block = block!(graph, [(0, 1), (0, -1)],);

    for direction in Direction::all() {
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
        let second_pos = second.pos();

        let condition = first.is_atom() && first.unwrap_atom().element == C;

        if condition {
            graph.put(second_pos, ComponentType::Element(H));
            hydrogen_correction(graph, first_pos);
            return true;
        }
    }

    false
}

fn carbonyl_extension(graph: &mut GridState) -> bool {
    let mut block = block!(graph, [(0, 1), (0, 2)],);

    for direction in Direction::all() {
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
        let second_pos = second.pos();

        let condition = second.is_atom() && second.unwrap_atom().element == C;

        if condition {
            graph.put(first_pos, ComponentType::Order(Double));
            hydrogen_correction(graph, second_pos);
            return true;
        }
    }

    false
}

fn hydrogen_fill(graph: &mut GridState, pos: Vec2) -> bool {
    let ptr = Pointer::new(graph, pos);
    let first_neighbor = match ptr.borrow() {
        Ok(it) => it.is_empty(),
        Err(_) => return false,
    };

    if first_neighbor {
        graph.put(pos, ComponentType::Element(H));
        true
    } else {
        false
    }
}

pub fn hydrogen_correction(graph: &mut GridState, pos: Vec2) {
    let bond_count = Pointer::new(graph, pos).bond_count().unwrap();
    let mut bonds_needed = graph.get(pos).unwrap().unwrap_atom().element.bond_number() - bond_count;

    for direction in Direction::all() {
        let adjusted = pos + direction.to_vec2();

        match bonds_needed.cmp(&0) {
            Ordering::Equal => return,
            Ordering::Less => {
                if hydrogen_remove(graph, adjusted) {
                    bonds_needed += 1;
                }
            }
            Ordering::Greater => {
                if hydrogen_fill(graph, adjusted) {
                    bonds_needed -= 1;
                }
            }
        }
    }
}

fn hydrogen_remove(graph: &mut GridState, pos: Vec2) -> bool {
    let cell = graph.get(pos).unwrap();

    if cell.is_atom() && cell.unwrap_atom().element == H {
        graph.put(pos, ComponentType::None);
        true
    } else {
        false
    }
}

pub fn bond_conversion(
    graph: &mut GridState,
    pos: Vec2,
    order: BondOrder,
    orient: BondOrientation,
) {
    let directions = if orient == Horiz {
        [Direction::Left, Direction::Right]
    } else {
        [Direction::Up, Direction::Down]
    };

    'outer: for direction in directions {
        let mut pos = pos;

        loop {
            pos += direction.to_vec2();
            let current_cell = match graph.get_mut(pos) {
                Ok(it) => it,
                Err(_) => continue 'outer,
            };

            match current_cell {
                Cell::Bond(it) if it.orient == orient => it.order = order,
                _ => continue 'outer,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn hydrogen_correction_leaves_full_bonds() {
        let mut graph = graph_with!(3, 3,
            [1, 1; A(C)],
            [0, 1; A(H)],
            [1, 0; A(H)],
            [1, 2; B(Double)],
            [2, 1; B(Single)],
        );
        hydrogen_correction(&mut graph, Vec2::xy(1, 1));
        let ptr = Pointer::new(&graph, Vec2::xy(1, 1));

        assert_eq!(ptr.bond_count(), Ok(4));
    }
}
