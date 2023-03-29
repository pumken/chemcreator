//! # Input
//!
//! The `input` module contains functions that interpret user input.

use crate::groups::debug_branches;
use crate::macros::invoke_macro;
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Element::{Cl, C, H, I, N, O};
use crate::molecule::{Atom, ComponentType, Element};
use crate::naming::name_molecule;
use crate::pointer::Pointer;
use crate::spatial::{FromVec2, GridState};
use crate::{AppState, Mode};
use ruscii::app::State;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::{Direction, Vec2};
use Element::{Br, F};

pub(crate) fn input_insert_mode(app_state: &State, state: &mut AppState, graph: &mut GridState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::B) => update(state, graph, ComponentType::Element(Br)),
            KeyEvent::Pressed(Key::C) => update(state, graph, ComponentType::Element(C)),
            KeyEvent::Pressed(Key::F) => update(state, graph, ComponentType::Element(F)),
            KeyEvent::Pressed(Key::H) => update(state, graph, ComponentType::Element(H)),
            KeyEvent::Pressed(Key::I) => update(state, graph, ComponentType::Element(I)),
            KeyEvent::Pressed(Key::L) => update(state, graph, ComponentType::Element(Cl)),
            KeyEvent::Pressed(Key::N) => update(state, graph, ComponentType::Element(N)),
            KeyEvent::Pressed(Key::O) => update(state, graph, ComponentType::Element(O)),
            KeyEvent::Pressed(Key::F5) => {
                graph.clear_all();
                update(state, graph, ComponentType::None);
            }
            KeyEvent::Pressed(Key::F7) => state.macros_enabled = !state.macros_enabled,
            KeyEvent::Pressed(Key::F12) => {
                state.debug = match debug_branches(graph) {
                    Ok(it) => it.to_string(),
                    Err(it) => it.to_string(),
                }
            }
            KeyEvent::Pressed(Key::Num1) => update(state, graph, ComponentType::Order(Single)),
            KeyEvent::Pressed(Key::Num2) => update(state, graph, ComponentType::Order(Double)),
            KeyEvent::Pressed(Key::Num3) => update(state, graph, ComponentType::Order(Triple)),
            KeyEvent::Pressed(Key::Num0) => {
                state.parent_chain_enabled = !state.parent_chain_enabled
            }
            KeyEvent::Pressed(Key::Backspace) => update(state, graph, ComponentType::None),
            KeyEvent::Pressed(Key::Right) => graph.move_cursor(Direction::Right),
            KeyEvent::Pressed(Key::Left) => graph.move_cursor(Direction::Left),
            KeyEvent::Pressed(Key::Up) => graph.move_cursor(Direction::Up),
            KeyEvent::Pressed(Key::Down) => graph.move_cursor(Direction::Down),
            KeyEvent::Pressed(Key::Esc) => state.mode = Mode::Normal,
            _ => {}
        }
    }
}

pub(crate) fn input_view_mode(app_state: &State, state: &mut AppState) {
    for key_event in app_state.keyboard().last_key_events() {
        match key_event {
            KeyEvent::Pressed(Key::Esc) => app_state.stop(),
            KeyEvent::Pressed(Key::F8) => state.mode = Mode::Insert,
            KeyEvent::Pressed(Key::Num0) => {
                state.parent_chain_enabled = !state.parent_chain_enabled
            }
            _ => (),
        }
    }
}

pub(crate) fn start_mode(app_state: &State, state: &mut AppState) {
    for key_event in app_state.keyboard().last_key_events() {
        if let KeyEvent::Pressed(_) = key_event {
            state.mode = Mode::Insert;
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

    let mut chain = None;

    (state.name, state.err) = match name_molecule(graph, &mut chain) {
        Ok(it) => (it, false),
        Err(it) => (it.to_string(), true),
    };

    if let Some(it) = chain {
        state.parent_chain = Some(get_chain_path(graph, it))
    } else {
        state.parent_chain = None
    }
}

fn get_chain_path(graph: &GridState, chain: Vec<Atom>) -> Vec<Vec2> {
    let mut out = vec![];

    for atoms in chain.windows(2) {
        let direction = Direction::from_points(atoms[0].pos, atoms[1].pos).unwrap();
        let mut ptr = Pointer::new(graph, atoms[0].pos);

        'inner: loop {
            out.push(ptr.borrow().unwrap().pos());
            ptr.move_ptr(direction);

            if ptr.pos == atoms[1].pos {
                break 'inner;
            }
        }
    }

    out.push(chain.last().unwrap().pos);
    out
}
