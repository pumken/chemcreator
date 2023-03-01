use ruscii::app::{App, State};
use ruscii::drawing::Pencil;
use ruscii::keyboard::{Key, KeyEvent};
use ruscii::spatial::Vec2;
use ruscii::terminal::Window;

use crate::grid::GridState;
use crate::input::{handle_insertion, handle_movement};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Symbol;

mod grid;
mod molecule;
mod r#macro;
mod input;

fn main() {
    let mut app = App::new();
    let mut graph = GridState::new(20, 10);
    let mut mode = Mode::Normal;
    let mut menu_name = "                ChemCreator                ".to_string();
    let mut menu_pos  = "      Written in Rust by Gavin Tran.       ".to_string();
    let mut menu_sym  = "To start, enter insert mode by pressing F8.".to_string();
    let mut menu_key   = "          You're using version 1.          ";

    app.run(|app_state: &mut State, window: &mut Window| {
        match mode {
            Mode::Insert => {
                for key_event in app_state.keyboard().last_key_events() {
                    if let Some(it) = handle_insertion(key_event, &mut graph) { menu_key = it; };
                    if let Some(it) = handle_movement(key_event, &mut graph) { menu_key = it; }
                    match key_event {
                        KeyEvent::Pressed(Key::Esc) => mode = Mode::Normal,
                        _ => {
                            menu_key = match key_event {
                                KeyEvent::Pressed(_) => { "Pressed" }
                                KeyEvent::Released(_) => { "Released" }
                            }
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
            Mode::Normal => {
                for key_event in app_state.keyboard().last_key_events() {
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
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid
        // Could optimize by not copying the entire graph, figure out later
        for cell in graph.cells.iter().flatten() {
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