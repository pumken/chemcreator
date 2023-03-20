//! # Input
//!
//! The `input` module contains functions that interpret user input.

use std::sync::mpsc;
use std::thread;
use std::time::Duration;
use async_std::task;
use crate::groups::debug_branches;
use crate::macros::invoke_macro;
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Element;
use crate::molecule::Element::{C, H, N, O};
use crate::naming::name_molecule;
use crate::spatial::GridState;
use crate::{AppState, Mode};
use ruscii::app::State;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::Direction;

pub(crate) fn input_insert_mode(app_state: &State, state: &mut AppState, graph: &mut GridState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::B) => {
                graph.put_atom(Element::Br);
                state.key = "B";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::C) => {
                graph.put_atom(C);
                state.key = "C";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::F) => {
                graph.put_atom(Element::F);
                state.key = "F";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::H) => {
                graph.put_atom(H);
                state.key = "H";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::I) => {
                graph.put_atom(Element::I);
                state.key = "I";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::L) => {
                graph.put_atom(Element::Cl);
                state.key = "L";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::N) => {
                graph.put_atom(N);
                state.key = "N";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::O) => {
                graph.put_atom(O);
                state.key = "O";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::F5) => {
                graph.clear_all();
                state.key = "F5";
                update(state, graph);
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
                graph.put_bond(Single);
                state.key = "1";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::Num2) => {
                graph.put_bond(Double);
                state.key = "2";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::Num3) => {
                graph.put_bond(Triple);
                state.key = "3";
                update(state, graph);
            }
            KeyEvent::Pressed(Key::Backspace) => {
                graph.clear_cell();
                state.key = "Backspace";
                update(state, graph);
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
pub(crate) fn update(state: &mut AppState, graph: &mut GridState) {
    if state.macros_enabled {
        invoke_macro(graph);
    }

    let (tx, rx) = mpsc::channel();
    let (txr, rxr) = mpsc::channel();
    state.rx = Some(rx);
    thread::spawn(move || {
        let gs: GridState = rxr.recv().unwrap();
        let result = name_molecule(&gs);
        tx.send(result).unwrap();
    });
    txr.send(graph.clone()).unwrap()
}
