mod grid;
mod molecule;
mod r#macro;

use ruscii::app::{App, State};
use ruscii::terminal::{Window};
use ruscii::drawing::{Pencil};
use ruscii::keyboard::{KeyEvent, Key};
use ruscii::spatial::{Vec2};
use crate::molecule::Atom::{C, H, O};
use crate::grid::GridState;
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::molecule::Symbol;
use crate::molecule::Bond;

fn main() {
    let mut app = App::new();
    let mut graph = GridState::new(20, 10);
    let mut mode = Mode::Normal;
    let mut menu_name = "                ChemCreator                ".to_string();
    let mut menu_pos  = "      Written in Rust by Gavin Tran.       ".to_string();
    let mut menu_sym  = "To start, enter insert mode by pressing F8.".to_string();
    let mut menu_key   = "          You're using version 1.          ";

    app.run(|app_state: &mut State, window: &mut Window| {
        for key_event in app_state.keyboard().last_key_events() {
            match mode {
                Mode::Insert => {
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => mode = Mode::Normal,
                        KeyEvent::Pressed(Key::Right) => {
                            graph.cursor.x += 1;
                            menu_key = "->";
                        }
                        KeyEvent::Pressed(Key::Left) => {
                            graph.cursor.x -= 1;
                            menu_key = "<-";
                        }
                        KeyEvent::Pressed(Key::Up) => {
                            graph.cursor.y -= 1;
                            menu_key = "^";
                        }
                        KeyEvent::Pressed(Key::Down) => {
                            graph.cursor.y += 1;
                            menu_key = "down";
                        }
                        KeyEvent::Pressed(Key::Backspace) => {
                            graph.insert(Symbol::None);
                            menu_key = "Backspace";
                        }
                        KeyEvent::Pressed(Key::Tab) => {
                            graph.insert(Symbol::Atom(C));
                            menu_key = "Tab";
                        }
                        KeyEvent::Pressed(Key::Enter) => {
                            graph.insert(Symbol::Atom(H));
                            menu_key = "Enter";
                        }
                        KeyEvent::Pressed(Key::F4) => {
                            graph.insert(Symbol::Atom(O));
                            menu_key = "F4";
                        }
                        KeyEvent::Pressed(Key::F1) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Single, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Single, Vert)));
                            }
                            menu_key = "F1";
                        }
                        KeyEvent::Pressed(Key::F2) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Double, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Double, Vert)));
                            }
                            menu_key = "F2";
                        }
                        KeyEvent::Pressed(Key::F3) => {
                            if graph.atom_adjacent() {
                                graph.insert(Symbol::Bond(Bond::new(Triple, Horiz)));
                            } else {
                                graph.insert(Symbol::Bond(Bond::new(Triple, Vert)));
                            }
                            menu_key = "F3";
                        }
                        _ => {
                            menu_key = match key_event {
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
                            menu_key = "F8";
                        }
                        _ => ()
                    }
                }
            }

            menu_name = graph.simple_counter();
            menu_pos = format!("({}, {})", graph.cursor.x, graph.cursor.y);
            menu_sym = match graph.current_cell().sym {
                Symbol::Atom(it) => { format!("Atom {}", it.symbol()) }
                Symbol::Bond(it) => {
                    match it.order {
                        Single => format!("Bond 1x"),
                        Double => format!("Bond 2x"),
                        Triple => format!("Bond 3x"),
                    }
                }
                Symbol::None => { format!("") }
            }
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid
        // Could optimize by not copying the entire graph, figure out later
        for cell in graph.clone().into_iter().flatten() {
            pencil.set_foreground(*&cell.sym.color());
            pencil.draw_text(&format!("{}", match &cell.sym {
                Symbol::Atom(it) => { it.symbol() }
                Symbol::Bond(it) => { it.symbol() }
                Symbol::None => match mode {
                    Mode::Insert => " • ",
                    Mode::Normal => "   ",
                }
            }), Vec2::xy(cell.pos.x * 3, cell.pos.y));
        }

        // Insert mode and cursor
        match mode {
            Mode::Insert => {
                pencil.draw_text("-- INSERT MODE --", Vec2::xy(graph.size.x * 3 + 3, 0));
                pencil.draw_text("<", Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y));
                pencil.draw_text(">", Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y));
            }
            _ => ()
        }

        // Menu
        pencil.draw_text(&format!("name | {}", menu_name), Vec2::xy(graph.size.x * 3 + 3, 1));
        pencil.draw_text(&format!(" pos | {}", menu_pos), Vec2::xy(graph.size.x * 3 + 3, 2));
        pencil.draw_text(&format!(" sym | {}", menu_sym), Vec2::xy(graph.size.x * 3 + 3, 3));
        pencil.draw_text(&format!(" key | {}", menu_key), Vec2::xy(graph.size.x * 3 + 3, 4));
    });
}

enum Mode {
    Insert,
    Normal,
}