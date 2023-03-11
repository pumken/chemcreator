//! # ChemCreator
//!
//! `chemcreator` is a binary crate that allows you to get information about an organic
//! molecule by building it in a text-based user interface.

#![warn(missing_docs)]

use ruscii::app::{App, State};
use ruscii::drawing::Pencil;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::Vec2;
use ruscii::terminal::Color::{Cyan, Red, White};
use ruscii::terminal::Window;
use crate::algorithm::{endpoint_carbons, name_molecule};

use crate::grid::{GridState, Invert};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Cell;
use crate::molecule::Element::{C, H, O};

mod grid;
mod molecule;
mod macros;
mod input;
mod algorithm;

fn main() {
    let mut app = App::new();
    let version = env!("CARGO_PKG_VERSION");
    let mut graph = GridState::new(20, 10);
    let mut mode = Mode::Start;
    let mut menu_name = "                ChemCreator                ".to_string();
    let mut menu_pos = "      Written in Rust by Gavin Tran.       ".to_string();
    let mut menu_sym = "To start, enter insert mode by pressing F8.".to_string();
    let mut menu_key = "";  // Overridden in Start mode to retain &str type
    let mut menu_err = "".to_string();
    let mut menu_debug = "".to_string();

    app.run(|app_state: &mut State, window: &mut Window| {
        match mode {
            Mode::Insert => {
                for key_event in app_state.keyboard().last_key_events() {
                    match key_event {
                        KeyEvent::Pressed(Key::C) => {
                            graph.put_atom(C);
                            menu_key = "C";
                        }
                        KeyEvent::Pressed(Key::H) => {
                            graph.put_atom(H);
                            menu_key = "H";
                        }
                        KeyEvent::Pressed(Key::O) => {
                            graph.put_atom(O);
                            menu_key = "O";
                        }
                        KeyEvent::Pressed(Key::Num1) => {
                            graph.put_bond(Single);
                            menu_key = "1";
                        }
                        KeyEvent::Pressed(Key::Num2) => {
                            graph.put_bond(Double);
                            menu_key = "2";
                        }
                        KeyEvent::Pressed(Key::Num3) => {
                            graph.put_bond(Triple);
                            menu_key = "3";
                        }
                        KeyEvent::Pressed(Key::Backspace) => {
                            graph.clear();
                            menu_key = "Backspace";
                        }
                        KeyEvent::Pressed(Key::Right) => {
                            graph.move_right();
                            menu_key = "→";
                        }
                        KeyEvent::Pressed(Key::Left) => {
                            graph.move_left();
                            menu_key = "←";
                        }
                        KeyEvent::Pressed(Key::Up) => {
                            graph.move_up();
                            menu_key = "↑";
                        }
                        KeyEvent::Pressed(Key::Down) => {
                            graph.move_down();
                            menu_key = "↓";
                        }
                        KeyEvent::Pressed(Key::Esc) => mode = Mode::Normal,
                        _ => {
                            menu_key = match key_event {
                                KeyEvent::Pressed(_) => { "Pressed" }
                                KeyEvent::Released(_) => { "" }
                            }
                        }
                    }
                }

                (menu_name, menu_err) = match name_molecule(&graph) {
                    Ok(it) => (it, "".to_string()),
                    Err(it) => ("unidentified".to_string(), it.to_string()),
                };
                menu_pos = format!("({}, {})", graph.cursor.x, graph.cursor.y);
                menu_sym = match graph.current_cell() {
                    Cell::Atom(it) => { format!("Atom {}", it.symbol()) }
                    Cell::Bond(it) => {
                        match it.order {
                            Single => format!("Bond 1x"),
                            Double => format!("Bond 2x"),
                            Triple => format!("Bond 3x"),
                        }
                    }
                    Cell::None(_) => { format!("") }
                };
                menu_debug = match endpoint_carbons(&graph) {
                    Ok(it) => {
                        it.iter()
                            .fold("".to_string(), |a, b| {
                                format!("{a} ({}, {}),", b.pos().x, b.pos().y)
                            })
                    }
                    Err(it) => it.to_string()
                };
            }
            Mode::Normal => {
                menu_pos = "".to_string();
                menu_sym = "".to_string();
                menu_key = "";
                menu_err = "".to_string();
                for key_event in app_state.keyboard().last_key_events() {
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => app_state.stop(),
                        KeyEvent::Pressed(Key::F8) => mode = Mode::Insert,
                        _ => ()
                    }
                }
            }
            Mode::Start => {
                for key_event in app_state.keyboard().last_key_events() {
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => app_state.stop(),
                        KeyEvent::Pressed(Key::F8) => mode = Mode::Insert,
                        _ => ()
                    }
                }
            }
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid and startup screen
        match mode {
            Mode::Start => {
                pencil.draw_text("┌────┐          ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 2).inv(&graph));
                pencil.draw_text("│  C │hem       ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 1).inv(&graph));
                pencil.draw_text("└────┼────┐     ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 0).inv(&graph));
                pencil.draw_text("     │ Cr │eator", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 1).inv(&graph));
                pencil.draw_text("     └────┘     ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 2).inv(&graph));
                pencil.draw_text("    Release 1   ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 3).inv(&graph));
            }
            _ => {
                for cell in graph.cells.iter().flatten() {
                    pencil.set_foreground(*&cell.color());
                    pencil.draw_text(&format!("{}", match &cell {
                        Cell::Atom(it) => { it.symbol() }
                        Cell::Bond(it) => { it.symbol() }
                        Cell::None(_) => match mode {
                            Mode::Insert => " • ",
                            Mode::Normal => "   ",
                            _ => "   "
                        }
                    }), Vec2::xy(cell.pos().x * 3, cell.pos().y).inv(&graph));
                }
            }
        }

        // Insert mode and cursor
        match mode {
            Mode::Insert => {
                pencil.draw_text("-- INSERT MODE --", Vec2::xy(graph.size.x * 3 + 3, 1));
                pencil.set_foreground(Cyan);
                pencil.draw_text("<", Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y).inv(&graph));
                pencil.draw_text(">", Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y).inv(&graph));
                pencil.set_foreground(White);
            }
            _ => ()
        }

        // Menu
        pencil.draw_text(match mode {
            Mode::Start => "",
            _ => "ChemCreator"
        }, Vec2::xy(graph.size.x * 3 + 3, 0));
        pencil.draw_text(&format!("name | {}", menu_name), Vec2::xy(graph.size.x * 3 + 3, 2));
        pencil.draw_text(&format!(" pos | {}", menu_pos), Vec2::xy(graph.size.x * 3 + 3, 3));
        pencil.draw_text(&format!(" sym | {}", menu_sym), Vec2::xy(graph.size.x * 3 + 3, 4));
        pencil.draw_text(&format!(" key | {}", match mode {
            Mode::Start => format!("        You're using version {version}.          ").to_string(),
            _ => menu_key.to_string()
        }), Vec2::xy(graph.size.x * 3 + 3, 5));
        pencil.set_foreground(Red);
        pencil.draw_text(&format!("{}", menu_err), Vec2::xy(graph.size.x * 3 + 3, 7));
        pencil.draw_text(&format!("{}", menu_debug), Vec2::xy(graph.size.x * 3 + 3, 8));
    });
}

/// Represents the mode the app is in at a given time.
enum Mode {
    Insert,
    Normal,
    Start,
}