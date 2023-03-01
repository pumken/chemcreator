use ruscii::keyboard::{Key, KeyEvent};
use crate::grid::GridState;
use crate::molecule::{Bond, Symbol};
use crate::molecule::Atom::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};

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
            graph.insert(Symbol::Bond(Bond::adjusted(Single, &graph)));
            Some("F1")
        }
        KeyEvent::Pressed(Key::F2) => {
            graph.insert(Symbol::Bond(Bond::adjusted(Double, &graph)));
            Some("F2")
        }
        KeyEvent::Pressed(Key::F3) => {
            graph.insert(Symbol::Bond(Bond::adjusted(Triple, &graph)));
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