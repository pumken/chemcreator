//! # Macros
//!
//! Not to be confused with Rust's `macro_rules!` declarations, the `macros` module contains
//! common actions that should be automatically performed for the user when they make an input.

use crate::molecule::BondOrder::Single;
use crate::molecule::Element::{C, H};
use crate::molecule::ComponentType;
use crate::pointer::Pointer;
use crate::spatial::{EnumAll, GridState, ToVec2};
use ruscii::spatial::Direction;

pub fn invoke_macro(graph: &mut GridState) {
    match graph
        .current_cell()
        .expect("cell should be within bounds")
        .comp()
    {
        ComponentType::Element(C) => fill_hydrogen(graph),
        ComponentType::None => {} // just to appease Clippy
        _ => {}
    }
}

pub fn fill_hydrogen(graph: &mut GridState) {
    for direction in Direction::all() {
        hydrogen_arm(graph, direction);
    }
}

fn hydrogen_arm(graph: &mut GridState, direction: Direction) {
    let mut ptr = Pointer::new(graph, graph.cursor);
    let first_neighbor = ptr.move_ptr(direction) && !ptr.borrow().unwrap().is_empty();
    let second_neighbor = ptr.move_ptr(direction) && !ptr.borrow().unwrap().is_empty();
    let pos = ptr.pos;

    if first_neighbor && second_neighbor {
        graph.put(pos, ComponentType::Element(H));
        graph.put(pos - direction.to_vec2(), ComponentType::Order(Single));
    } else if first_neighbor && !second_neighbor {
        graph.put(pos, ComponentType::Element(H));
    }
}
