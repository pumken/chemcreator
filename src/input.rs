//! # Input
//!
//! The `input` module contains functions that interpret user input.

use crate::groups::debug_branches;
use crate::macros::invoke_macro;
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Element::{Cl, C, H, I, N, O};
use crate::molecule::{ComponentType, Element};
use crate::naming::name_molecule;
use crate::spatial::GridState;
use crate::{AppState, Mode};
use ruscii::app::State;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::Direction;
use Element::{Br, F};

pub(crate) fn input_insert_mode(app_state: &State, state: &mut AppState, graph: &mut GridState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::B) => {
                state.key = "B";
                update(state, graph, ComponentType::Element(Br));
            }
            KeyEvent::Pressed(Key::C) => {
                state.key = "C";
                update(state, graph, ComponentType::Element(C));
            }
            KeyEvent::Pressed(Key::F) => {
                state.key = "F";
                update(state, graph, ComponentType::Element(F));
            }
            KeyEvent::Pressed(Key::H) => {
                state.key = "H";
                update(state, graph, ComponentType::Element(H));
            }
            KeyEvent::Pressed(Key::I) => {
                state.key = "I";
                update(state, graph, ComponentType::Element(I));
            }
            KeyEvent::Pressed(Key::L) => {
                state.key = "L";
                update(state, graph, ComponentType::Element(Cl));
            }
            KeyEvent::Pressed(Key::N) => {
                state.key = "N";
                update(state, graph, ComponentType::Element(N));
            }
            KeyEvent::Pressed(Key::O) => {
                state.key = "O";
                update(state, graph, ComponentType::Element(O));
            }
            KeyEvent::Pressed(Key::F5) => {
                graph.clear_all();
                state.key = "F5";
                update(state, graph, ComponentType::None);
            }
            KeyEvent::Pressed(Key::F7) => {
                state.macros_enabled = !state.macros_enabled;
                state.key = "F7";
            }
            KeyEvent::Pressed(Key::F12) => {
                state.debug = match debug_branches(graph) {
                    Ok(it) => it.to_string(),
                    Err(it) => it.to_string(),
                };
                state.key = "F12";
            }
            KeyEvent::Pressed(Key::Num1) => {
                state.key = "1";
                update(state, graph, ComponentType::Order(Single));
            }
            KeyEvent::Pressed(Key::Num2) => {
                state.key = "2";
                update(state, graph, ComponentType::Order(Double));
            }
            KeyEvent::Pressed(Key::Num3) => {
                state.key = "3";
                update(state, graph, ComponentType::Order(Triple));
            }
            KeyEvent::Pressed(Key::Backspace) => {
                state.key = "Backspace";
                update(state, graph, ComponentType::None);
            }
            KeyEvent::Pressed(Key::Right) => {
                graph.move_cursor(Direction::Right);
                state.key = "→";
            }
            KeyEvent::Pressed(Key::Left) => {
                graph.move_cursor(Direction::Left);
                state.key = "←";
            }
            KeyEvent::Pressed(Key::Up) => {
                graph.move_cursor(Direction::Up);
                state.key = "↑";
            }
            KeyEvent::Pressed(Key::Down) => {
                graph.move_cursor(Direction::Down);
                state.key = "↓";
            }
            KeyEvent::Pressed(Key::Esc) => state.mode = Mode::Normal,
            _ => {
                state.key = match key_event {
                    KeyEvent::Pressed(_) => "Pressed",
                    KeyEvent::Released(_) => "",
                }
            }
        }
    }
}

pub(crate) fn input_view_mode(app_state: &State, state: &mut AppState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::Esc) => app_state.stop(),
            KeyEvent::Pressed(Key::F8) => {
                state.mode = Mode::Insert;
                state.name = "".to_string()
            }
            _ => (),
        }
    }
}

//noinspection RsBorrowChecker
pub(crate) fn update(state: &mut AppState, graph: &mut GridState, comp: ComponentType) {
    let previous = graph
        .current_cell()
        .expect("cell should be within bounds")
        .comp();

    match comp {
        ComponentType::Element(it) => graph.put_atom(it),
        ComponentType::Order(it) => graph.put_bond(it),
        ComponentType::None => graph.clear_cell(),
    }

    if state.macros_enabled {
        invoke_macro(graph, comp, previous);
    }

    (state.name, state.err) = match name_molecule(graph) {
        Ok(it) => (it, "".to_string()),
        Err(it) => ("unidentified".to_string(), it.to_string()),
    };
}
