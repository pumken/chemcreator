use ruscii::keyboard::{Key, KeyEvent};
use crate::grid::GridState;
use crate::molecule::{Bond, Symbol};
use crate::molecule::Atom::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

pub(crate) fn handle_insertion(key_event: &KeyEvent, graph: &mut GridState) -> Option<&'static str> {
    match key_event {
        KeyEvent::Pressed(Key::Tab) => {
            graph.insert(Symbol::Atom(C));
            Some("Tab")
        }
        KeyEvent::Pressed(Key::Enter) => {
            graph.insert(Symbol::Atom(H));
            Some("Enter")
        }
        KeyEvent::Pressed(Key::F4) => {
            graph.insert(Symbol::Atom(O));
            Some("F4")
        }
        KeyEvent::Pressed(Key::F1) => {
            if graph.atom_adjacent() {
                graph.insert(Symbol::Bond(Bond::new(Single, Horiz)));
            } else {
                graph.insert(Symbol::Bond(Bond::new(Single, Vert)));
            }
            Some("F1")
        }
        KeyEvent::Pressed(Key::F2) => {
            if graph.atom_adjacent() {
                graph.insert(Symbol::Bond(Bond::new(Double, Horiz)));
            } else {
                graph.insert(Symbol::Bond(Bond::new(Double, Vert)));
            }
            Some("F2")
        }
        KeyEvent::Pressed(Key::F3) => {
            if graph.atom_adjacent() {
                graph.insert(Symbol::Bond(Bond::new(Triple, Horiz)));
            } else {
                graph.insert(Symbol::Bond(Bond::new(Triple, Vert)));
            }
            Some("F3")
        }
        KeyEvent::Pressed(Key::Backspace) => {
            graph.insert(Symbol::None);
            Some("Backspace")
        }
        _ => None
    }
}

pub(crate) fn handle_movement(key_event: &KeyEvent, graph: &mut GridState) -> Option<&'static str> {
    match key_event {
        KeyEvent::Pressed(Key::Right) => {
            graph.cursor.x += 1;
            Some("->")
        }
        KeyEvent::Pressed(Key::Left) => {
            graph.cursor.x -= 1;
            Some("<-")
        }
        KeyEvent::Pressed(Key::Up) => {
            graph.cursor.y -= 1;
            Some("↑")
        }
        KeyEvent::Pressed(Key::Down) => {
            graph.cursor.y += 1;
            Some("↓")
        }
        _ => None
    }
}