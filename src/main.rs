//! # ChemCreator
//!
//! `chemcreator` is a binary crate that allows you to get information about an organic
//! molecule by building it in a text-based user interface.

#![warn(missing_docs)]

use crate::gui::draw_grid_box;
use crate::input::{input_insert_mode, input_view_mode, start_mode};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::{Cell, ComponentType};
use crate::spatial::GridState;
use crate::Mode::{Display, Insert};
use ruscii::app::{App, State};
use ruscii::drawing::Pencil;
use ruscii::spatial::Vec2;
use ruscii::terminal::Color::{Red, White, Yellow};
use ruscii::terminal::Window;

mod chain;
mod groups;
mod gui;
mod input;
mod macros;
mod molecule;
mod naming;
mod numerics;
mod pointer;
mod spatial;
mod validation;

fn main() {
    let mut app = App::new();
    let version = env!("CARGO_PKG_VERSION");
    let mut graph = GridState::new(27, 11);
    let mut state = AppState::default();

    app.run(|app_state: &mut State, window: &mut Window| {
        let window_size = window.size();

        match state.mode {
            Insert => {
                input_insert_mode(app_state, &mut state, &mut graph);

                state.pos = format!("({}, {})", graph.cursor.x, graph.cursor.y);
                state.sym = match graph
                    .current_cell()
                    .expect("current cell should be within bounds")
                {
                    Cell::Atom(it) => ComponentType::Element(it.element),
                    Cell::Bond(it) => ComponentType::Order(it.order),
                    Cell::None(_) => ComponentType::None,
                };
            }
            Display => {
                state.pos = "".to_string();
                input_view_mode(app_state, &mut state)
            }
            Mode::Start => start_mode(app_state, &mut state),
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid and startup screen
        pencil.set_origin(Vec2::xy(6, 0));
        match state.mode {
            Mode::Start => {
                gui::draw_logo(&mut pencil, graph.size, Vec2::y(2));
                gui::draw_start_message(&mut pencil, Vec2::xy(graph.size.x * 3 + 5, 1));
            }
            _ => gui::draw_grid(
                &mut pencil,
                &mut graph,
                if state.parent_chain_enabled {
                    &state.parent_chain
                } else {
                    &None
                },
                Vec2::y(2),
                state.mode,
            ),
        }

        draw_grid_box(&mut pencil, graph.size, Vec2::xy(-1, 1));

        // Insert mode and cursor
        if let Insert = state.mode {
            gui::draw_insert_mode(&mut pencil, Vec2::xy(graph.size.x * 3 / 2, 1));
            gui::draw_cursor(&mut pencil, &graph, Vec2::y(2));
        } else if let Mode::Start = state.mode {
            pencil.draw_text(
                &format!(" {version} "),
                Vec2::xy(
                    graph.size.x * 3 - version.len() as i32 - 3,
                    graph.size.y + 2,
                ),
            );
        }

        // Menu
        let pos_string = format!(" {} ", state.pos);
        let pos = if state.mode == Insert {
            &pos_string
        } else {
            ""
        };
        pencil
            .draw_center_text(
                match state.mode {
                    Mode::Start => "",
                    _ => " ChemCreator ",
                },
                Vec2::xy(
                    graph.size.x * 3 / 2,
                    if state.mode == Display { 1 } else { 0 },
                ),
            )
            .draw_text(pos, Vec2::xy(0, 1))
            .set_foreground(Red)
            .draw_text(&state.debug.to_string(), Vec2::xy(graph.size.x * 3 + 3, 8))
            .set_foreground(White);

        if matches!(state.mode, Insert) {
            let sym_type = match state.sym {
                ComponentType::Element(_) => " Atom ",
                ComponentType::Order(_) => " Bond ",
                ComponentType::None => "",
            };
            let sym_sym = match state.sym {
                ComponentType::Element(it) => format!("{} ", it.symbol()),
                ComponentType::Order(Single) => "1x ".to_string(),
                ComponentType::Order(Double) => "2x ".to_string(),
                ComponentType::Order(Triple) => "3x ".to_string(),
                ComponentType::None => "".to_string(),
            };
            pencil
                .set_foreground(state.sym.color())
                .draw_text(
                    &sym_sym,
                    Vec2::xy(graph.size.x * 3 - sym_sym.len() as i32 - 1, 1),
                )
                .set_foreground(White)
                .draw_text(
                    sym_type,
                    Vec2::xy(
                        graph.size.x * 3 - sym_type.len() as i32 - sym_sym.len() as i32 - 1,
                        1,
                    ),
                );
        }

        if state.macros_enabled && state.mode == Insert {
            pencil
                .set_foreground(Yellow)
                .draw_text(" macros enabled ", Vec2::y(graph.size.y + 2))
                .set_foreground(White);
        }

        // Statistics
        if state.name.is_empty() {
            return;
        }

        pencil.move_origin(Vec2::xy(graph.size.x * 3 + 5, 1));
        gui::draw_wrapped_name(&mut graph, &mut state, window_size, &mut pencil);

        if matches!(state.mode, Display) {
            gui::draw_statistics(&graph, &mut state, &mut pencil);
        }
    });
}

/// Represents the mode the app is in at a given time.
#[derive(Copy, Clone, Default, Debug, PartialEq)]
enum Mode {
    Insert,
    Display,
    #[default]
    Start,
}

/// Contains the running state of the app not including the grid.
#[derive(Default)]
struct AppState {
    mode: Mode,
    name: String,
    pos: String,
    sym: ComponentType,
    err: bool,
    debug: String,
    macros_enabled: bool,
    parent_chain_enabled: bool,
    parent_chain: Option<Vec<Vec2>>,
}

#[cfg(test)]
mod test_utils {
    use crate::molecule::BondOrder::Single;
    use crate::molecule::Element::{C, H, O};
    use crate::molecule::{Atom, BondOrder, Cell, Element};
    use crate::spatial::GridState;
    use crate::test_utils::GW::{A, B};
    use ruscii::spatial::Vec2;

    /// An enum used to make it easier to construct [`GridState`]s with the `graph_with` macro.
    pub(crate) enum GW {
        A(Element),
        B(BondOrder),
    }

    /// Creates a [`GridState`] with the given `vals` at (`x`, `y`). Used with the [`GW`] enum.
    #[macro_export]
    macro_rules! graph_with {
        ($width:expr, $height:expr) => {
            GridState::new($width, $height)
        };
        ($width:expr, $height:expr, $([$x:expr, $y:expr; $val:expr],)*) => {{
            let mut graph = GridState::new($width, $height);
            $(
            graph.cursor = Vec2::xy($x, $y);
            match $val {
                A(it) => graph.put_atom(it),
                B(it) => graph.put_bond(it),
            }
            )*
            graph
        }};
    }

    #[test]
    fn unwrap_atom_returns_atom() {
        let graph = graph_with!(1, 1,
            [0, 0; A(C)],
        );
        let cell = graph.get(Vec2::zero()).unwrap();
        let atom = cell.unwrap_atom();

        assert_eq!(
            atom,
            Atom {
                element: C,
                pos: Vec2::zero(),
            }
        );
    }

    #[test]
    fn graph_with_empty() {
        let graph = graph_with!(5, 5);

        let any_filled_cells = graph
            .cells
            .iter()
            .flatten()
            .any(|cell| !matches!(cell, Cell::None(_)));

        assert!(!any_filled_cells);
    }

    #[test]
    fn graph_with_generates_gridstate() {
        let a = {
            let mut graph = GridState::new(3, 3);
            graph.cursor = Vec2::xy(0, 1);
            graph.put_atom(C);
            graph.cursor = Vec2::xy(1, 1);
            graph.put_bond(Single);
            graph.cursor = Vec2::xy(2, 1);
            graph.put_atom(O);
            graph.cursor = Vec2::xy(2, 0);
            graph.put_atom(H);
            graph
        };
        let b = graph_with!(3, 3,
            [0, 1; A(C)],
            [1, 1; B(Single)],
            [2, 1; A(O)],
            [2, 0; A(H)],
        );

        assert_eq!(a, b)
    }
}
