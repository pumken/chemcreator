//! # Input
//!
//! The `input` module contains functions that interpret user input.

use ruscii::app::State;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::Direction;
use crate::algorithm::debug_chain;
use crate::{AppState, Mode};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Element::{C, H, O};
use crate::spatial::GridState;

pub(crate) fn input_insert_mode(app_state: &State, state: &mut AppState, graph: &mut GridState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::C) => {
                graph.put_atom(C);
                state.key = "C";
            }
            KeyEvent::Pressed(Key::H) => {
                graph.put_atom(H);
                state.key = "H";
            }
            KeyEvent::Pressed(Key::O) => {
                graph.put_atom(O);
                state.key = "O";
            }
            KeyEvent::Pressed(Key::F5) => {
                state.debug = match debug_chain(&graph) {
                    Ok(it) => {
                        it.iter()
                            .fold("".to_string(), |a, b| {
                                format!("{a} ({}, {}),", b.pos.x, b.pos.y)
                            })
                    }
                    Err(it) => it.to_string()
                };
            }
            KeyEvent::Pressed(Key::Num1) => {
                graph.put_bond(Single);
                state.key = "1";
            }
            KeyEvent::Pressed(Key::Num2) => {
                graph.put_bond(Double);
                state.key = "2";
            }
            KeyEvent::Pressed(Key::Num3) => {
                graph.put_bond(Triple);
                state.key = "3";
            }
            KeyEvent::Pressed(Key::Backspace) => {
                graph.clear();
                state.key = "Backspace";
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
                    KeyEvent::Pressed(_) => { "Pressed" }
                    KeyEvent::Released(_) => { "" }
                }
            }
        }
    }
}

pub(crate) fn input_view_mode(app_state: &State, state: &mut AppState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::Esc) => app_state.stop(),
            KeyEvent::Pressed(Key::F8) => state.mode = Mode::Insert,
            _ => ()
        }
    }
}
