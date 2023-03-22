//! # ChemCreator
//!
//! `chemcreator` is a binary crate that allows you to get information about an organic
//! molecule by building it in a text-based user interface.

#![warn(missing_docs)]

use crate::input::{input_insert_mode, input_view_mode};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::Cell;
use crate::spatial::{GridState, Invert};
use crate::Mode::Insert;
use ruscii::app::{App, State};
use ruscii::drawing::Pencil;
use ruscii::spatial::Vec2;
use ruscii::terminal::Color::{Cyan, Red, White};
use ruscii::terminal::{Color, Window};

mod chain;
mod groups;
mod input;
mod macros;
mod molecule;
mod naming;
mod pointer;
mod spatial;
mod validation;

fn main() {
    let mut app = App::new();
    let version = env!("CARGO_PKG_VERSION");
    let mut graph = GridState::new(21, 11);
    let mut state = AppState::default();

    app.run(|app_state: &mut State, window: &mut Window| {
        match state.mode {
            Insert => {
                input_insert_mode(app_state, &mut state, &mut graph);

                state.pos = format!("({}, {})", graph.cursor.x, graph.cursor.y);
                state.sym = match graph
                    .current_cell()
                    .expect("current cell should be within bounds")
                {
                    Cell::Atom(it) => {
                        format!("Atom {}", it.symbol())
                    }
                    Cell::Bond(it) => match it.order {
                        Single => "Bond 1x".to_string(),
                        Double => "Bond 2x".to_string(),
                        Triple => "Bond 3x".to_string(),
                    },
                    Cell::None(_) => String::new(),
                };
            }
            Mode::Normal => {
                state.pos = "".to_string();
                state.sym = "".to_string();
                state.key = "";
                state.err = "".to_string();
                input_view_mode(app_state, &mut state)
            }
            Mode::Start => input_view_mode(app_state, &mut state),
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid and startup screen
        match state.mode {
            Mode::Start => {
                pencil
                    .draw_text(
                        "┌────┐          ",
                        Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 2).inv(&graph),
                    )
                    .draw_text(
                        "│  C │hem       ",
                        Vec2::xy(graph.size.x + 2, graph.size.y / 2 + 1).inv(&graph),
                    )
                    .draw_text(
                        "└────┼────┐     ",
                        Vec2::xy(graph.size.x + 2, graph.size.y / 2).inv(&graph),
                    )
                    .draw_text(
                        "     │ Cr │eator",
                        Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 1).inv(&graph),
                    )
                    .draw_text(
                        "     └────┘     ",
                        Vec2::xy(graph.size.x + 2, graph.size.y / 2 - 2).inv(&graph),
                    );
            }
            _ => {
                for cell in graph.cells.iter().flatten() {
                    pencil.set_foreground(cell.color()).draw_text(
                        match &cell {
                            Cell::Atom(it) => it.symbol(),
                            Cell::Bond(it) => it.symbol(),
                            Cell::None(_) => match state.mode {
                                Mode::Insert => " • ",
                                Mode::Normal => "   ",
                                _ => "   ",
                            },
                        },
                        Vec2::xy(cell.pos().x * 3, cell.pos().y).inv(&graph),
                    );
                }
            }
        }

        // Insert mode and cursor
        if let Mode::Insert = state.mode {
            pencil
                .draw_text("-- INSERT MODE --", Vec2::xy(graph.size.x * 3 + 3, 1))
                .set_foreground(Cyan)
                .draw_text(
                    "<",
                    Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y).inv(&graph),
                )
                .draw_text(
                    ">",
                    Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y).inv(&graph),
                )
                .set_foreground(White);
        }

        // Menu
        pencil
            .draw_text(
                match state.mode {
                    Mode::Start => "",
                    _ => "ChemCreator",
                },
                Vec2::xy(graph.size.x * 3 + 3, 0),
            )
            .draw_text(
                &format!("name | {}", state.name),
                Vec2::xy(graph.size.x * 3 + 3, 2),
            )
            .draw_text(
                &format!(" pos | {}", state.pos),
                Vec2::xy(graph.size.x * 3 + 3, 3),
            )
            .draw_text(
                &format!(" sym | {}", state.sym),
                Vec2::xy(graph.size.x * 3 + 3, 4),
            )
            .draw_text(
                &format!(
                    " key | {}",
                    match state.mode {
                        Mode::Start => format!("        You're using version {version}.          "),
                        _ => state.key.to_string(),
                    }
                ),
                Vec2::xy(graph.size.x * 3 + 3, 5),
            )
            .set_foreground(Red)
            .draw_text(&state.err.to_string(), Vec2::xy(graph.size.x * 3 + 3, 7))
            .draw_text(&state.debug.to_string(), Vec2::xy(graph.size.x * 3 + 3, 8));

        if state.macros_enabled && state.mode == Insert {
            pencil
                .set_foreground(Color::Yellow)
                .draw_text("Macro mode enabled.", Vec2::xy(graph.size.x * 3 + 3, 8));
        }
    });
}

/// Represents the mode the app is in at a given time.
#[derive(Copy, Clone, Debug, PartialEq)]
enum Mode {
    Insert,
    Normal,
    Start,
}

/// Contains the running state of the app not including the grid.
struct AppState {
    mode: Mode,
    name: String,
    pos: String,
    sym: String,
    key: &'static str,
    err: String,
    debug: String,
    macros_enabled: bool,
}

impl Default for AppState {
    fn default() -> Self {
        AppState {
            mode: Mode::Start,
            name: "                ChemCreator                ".to_string(),
            pos: "      Written in Rust by Gavin Tran.       ".to_string(),
            sym: "To start, enter insert mode by pressing F8.".to_string(),
            key: "", // Overridden in Start mode to retain &str type
            err: "".to_string(),
            debug: "".to_string(),
            macros_enabled: false,
        }
    }
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
                pos: Vec2::zero()
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
