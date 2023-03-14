//! # ChemCreator
//!
//! `chemcreator` is a binary crate that allows you to get information about an organic
//! molecule by building it in a text-based user interface.

#![warn(missing_docs)]

use ruscii::app::{App, State};
use ruscii::drawing::Pencil;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::{Direction, Vec2};
use ruscii::terminal::Color::{Cyan, Red, White};
use ruscii::terminal::Window;
use crate::algorithm::{debug_chain, name_molecule};
use crate::input::{input_insert_mode, input_view_mode};

use crate::spatial::{GridState, Invert};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Cell;
use crate::molecule::Element::{C, H, O};

mod spatial;
mod molecule;
mod macros;
mod input;
mod algorithm;
mod pointer;

fn main() {
    let mut app = App::new();
    let version = env!("CARGO_PKG_VERSION");
    let mut graph = GridState::new(20, 10);
    let mut state = AppState::default();

    app.run(|app_state: &mut State, window: &mut Window| {
        match state.mode {
            Mode::Insert => {
                input_insert_mode(app_state, &mut state, &mut graph);

                (state.menu_name, state.menu_err) = match name_molecule(&graph) {
                    Ok(it) => (it, "".to_string()),
                    Err(it) => ("unidentified".to_string(), it.to_string()),
                };
                state.menu_pos = format!("({}, {})", graph.cursor.x, graph.cursor.y);
                state.menu_sym = match graph.current_cell() {
                    Cell::Atom(it) => { format!("Atom {}", it.symbol()) }
                    Cell::Bond(it) => {
                        match it.order {
                            Single => "Bond 1x".to_string(),
                            Double => "Bond 2x".to_string(),
                            Triple => "Bond 3x".to_string(),
                        }
                    }
                    Cell::None(_) => { String::new() }
                };
            }
            Mode::Normal => {
                state.menu_pos = "".to_string();
                state.menu_sym = "".to_string();
                state.menu_key = "";
                state.menu_err = "".to_string();
                input_view_mode(app_state, &mut state)
            }
            Mode::Start => input_view_mode(app_state, &mut state),
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid and startup screen
        match state.mode {
            Mode::Start => {
                pencil.draw_text("┌────┐          ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 2).inv(&graph));
                pencil.draw_text("│  C │hem       ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 1).inv(&graph));
                pencil.draw_text("└────┼────┐     ", Vec2::xy(graph.size.x + 2, graph.size.y / 2).inv(&graph));
                pencil.draw_text("     │ Cr │eator", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 1).inv(&graph));
                pencil.draw_text("     └────┘     ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 2).inv(&graph));
                pencil.draw_text("    Release 1   ", Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 3).inv(&graph));
            }
            _ => {
                for cell in graph.cells.iter().flatten() {
                    pencil.set_foreground(*&cell.color());
                    pencil.draw_text(match &cell {
                        Cell::Atom(it) => { it.symbol() }
                        Cell::Bond(it) => { it.symbol() }
                        Cell::None(_) => match state.mode {
                            Mode::Insert => " • ",
                            Mode::Normal => "   ",
                            _ => "   "
                        }
                    }, Vec2::xy(cell.pos().x * 3, cell.pos().y).inv(&graph));
                }
            }
        }

        // Insert mode and cursor
        if let Mode::Insert = state.mode {
            pencil.draw_text("-- INSERT MODE --", Vec2::xy(graph.size.x * 3 + 3, 1));
            pencil.set_foreground(Cyan);
            pencil.draw_text("<", Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y).inv(&graph));
            pencil.draw_text(">", Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y).inv(&graph));
            pencil.set_foreground(White);
        }

        // Menu
        pencil.draw_text(match state.mode {
            Mode::Start => "",
            _ => "ChemCreator"
        }, Vec2::xy(graph.size.x * 3 + 3, 0));
        pencil.draw_text(&format!("name | {}", state.menu_name), Vec2::xy(graph.size.x * 3 + 3, 2));
        pencil.draw_text(&format!(" pos | {}", state.menu_pos), Vec2::xy(graph.size.x * 3 + 3, 3));
        pencil.draw_text(&format!(" sym | {}", state.menu_sym), Vec2::xy(graph.size.x * 3 + 3, 4));
        pencil.draw_text(&format!(" key | {}", match state.mode {
            Mode::Start => format!("        You're using version {version}.          "),
            _ => state.menu_key.to_string()
        }), Vec2::xy(graph.size.x * 3 + 3, 5));
        pencil.set_foreground(Red);
        pencil.draw_text(&state.menu_err.to_string(), Vec2::xy(graph.size.x * 3 + 3, 7));
        pencil.draw_text(&state.menu_debug.to_string(), Vec2::xy(graph.size.x * 3 + 3, 8));
    });
}

/// Represents the mode the app is in at a given time.
enum Mode {
    Insert,
    Normal,
    Start,
}

struct AppState {
    mode: Mode,
    menu_name: String,
    menu_pos: String,
    menu_sym: String,
    menu_key: &'static str,
    menu_err: String,
    menu_debug: String,
}

impl Default for AppState {
    fn default() -> Self {
        AppState {
            mode: Mode::Start,
            menu_name: "                ChemCreator                ".to_string(),
            menu_pos:  "      Written in Rust by Gavin Tran.       ".to_string(),
            menu_sym:  "To start, enter insert mode by pressing F8.".to_string(),
            menu_key: "",  // Overridden in Start mode to retain &str type
            menu_err: "".to_string(),
            menu_debug: "".to_string(),
        }
    }
}