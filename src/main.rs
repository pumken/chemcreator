mod grid;
mod molecule;
mod r#macro;

use ruscii::app::{App, State};
use ruscii::terminal::{Window};
use ruscii::drawing::{Pencil};
use ruscii::keyboard::{KeyEvent, Key};
use ruscii::spatial::{Vec2};
use ruscii::terminal::Color::{LightGrey, Red, White};
use crate::molecule::Atom::{C, H, O};
use crate::grid::GridState;
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::molecule::Symbol;
use crate::molecule::Bond;

fn main() {
    let mut species_name = "add components".to_string();
    let mut graph = GridState::new(20, 10);
    let mut app = App::new();
    let mut key = "start";
    let mut mode = Mode::Insert;

    app.run(|app_state: &mut State, window: &mut Window| {
        for key_event in app_state.keyboard().last_key_events() {
            match mode {
                Mode::Insert => {
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => mode = Mode::Normal,
                        KeyEvent::Pressed(Key::Right) => {
                            graph.cursor.x += 1;
                            key = "->";
                        }
                        KeyEvent::Pressed(Key::Left) => {
                            graph.cursor.x -= 1;
                            key = "<-";
                        }
                        KeyEvent::Pressed(Key::Up) => {
                            graph.cursor.y -= 1;
                            key = "^";
                        }
                        KeyEvent::Pressed(Key::Down) => {
                            graph.cursor.y += 1;
                            key = "down";
                        }
                        KeyEvent::Pressed(Key::Backspace) => {
                            graph.insert(Symbol::None);
                            key = "Backspace";
                        }
                        KeyEvent::Pressed(Key::Tab) => {
                            graph.insert(Symbol::Atom(C));
                            key = "Tab";
                        }
                        KeyEvent::Pressed(Key::Enter) => {
                            graph.insert(Symbol::Atom(H));
                            key = "Enter";
                        }
                        KeyEvent::Pressed(Key::F4) => {
                            graph.insert(Symbol::Atom(O));
                            key = "F4";
                        }
                        KeyEvent::Pressed(Key::F1) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Single, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Single, Vert)));
                            }
                            key = "F1";
                        }
                        KeyEvent::Pressed(Key::F2) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Double, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Double, Vert)));
                            }
                            key = "F2";
                        }
                        KeyEvent::Pressed(Key::F3) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Triple, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Triple, Vert)));
                            }
                            key = "F3";
                        }
                        _ => {
                            key = match key_event {
                                KeyEvent::Pressed(_) => { "Pressed" }
                                KeyEvent::Released(_) => { "Released" }
                            }
                        }
                    }
                }
                Mode::Normal => {
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => app_state.stop(),
                        KeyEvent::Pressed(Key::Q) => app_state.stop(),
                        KeyEvent::Pressed(Key::F8) => {
                            mode = Mode::Insert;
                            key = "F8";
                        }
                        _ => {
                            key = match key_event {
                                KeyEvent::Pressed(_) => { "Pressed" }
                                KeyEvent::Released(_) => { "Released" }
                            }
                        }
                    }
                }
            }
            species_name = graph.simple_counter();
        }


        let mut pencil = Pencil::new(window.canvas_mut());
        for cell in graph.cells.iter() {
            pencil.set_foreground(
                match &cell.sym {
                    Symbol::Atom(it) => match it {
                        C => LightGrey,
                        O => Red,
                        _ => White
                    },
                    _ => White
                }
            );
            pencil.draw_text(&format!("{}", match &cell.sym {
                Symbol::Atom(it) => { it.symbol() }
                Symbol::Bond(it) => { it.symbol() }
                Symbol::None => match mode {
                    Mode::Insert => " • ",
                    Mode::Normal => "   ",
                }
            }), Vec2::xy(cell.pos.x * 3, cell.pos.y));
            pencil.set_foreground(White);
        }
        let current_cell = graph.current_cell();
        if current_cell.pos.x == graph.cursor.x && current_cell.pos.y == graph.cursor.y {
            pencil.draw_text(&format!("cells[{}][{}] = {}", current_cell.pos.x, current_cell.pos.y, match &current_cell.sym {
                Symbol::Atom(it) => { format! {"Atom::{}", it.symbol()} }
                Symbol::Bond(it) => { format! {"Bond::{}", it.symbol()} }
                Symbol::None => { format!("None") }
            }), Vec2::xy(graph.size.x * 3 + 3, 7));
        }
        match mode {
            Mode::Insert => {
                pencil.draw_text("-- INSERT MODE --", Vec2::xy(graph.size.x * 3 + 3, 0));
                pencil.draw_text("<", Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y));
                pencil.draw_text(">", Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y));
            }
            _ => ()
        }
        pencil.draw_text(&format!("Name | {}", species_name), Vec2::xy(graph.size.x * 3 + 3, 1));
        pencil.draw_text(&format!("   x | {}", graph.cursor.x), Vec2::xy(graph.size.x * 3 + 3, 2));
        pencil.draw_text(&format!("   y | {}", graph.cursor.y), Vec2::xy(graph.size.x * 3 + 3, 3));
        pencil.draw_text(&format!(" sym | {}", match graph.current_cell().sym {
            Symbol::Atom(it) => { format!("Atom {}", it.symbol()) }
            Symbol::Bond(it) => {
                match it.order {
                    Single => format!("Bond 1x"),
                    Double => format!("Bond 2x"),
                    Triple => format!("Bond 3x"),
                }
            }
            Symbol::None => { format!("") }
        }), Vec2::xy(graph.size.x * 3 + 3, 4));
        pencil.draw_text(&format!(" key | {}", key), Vec2::xy(graph.size.x * 3 + 3, 5));
    });
}

enum Mode {
    Insert,
    Normal,
}