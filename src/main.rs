//! # ChemCreator
//!
//! `chemcreator` is a binary crate that allows you to get information about an organic
//! molecule by building it in a text-based user interface.

#![warn(missing_docs)]

use Color::Yellow;
use crate::input::{input_insert_mode, input_view_mode, start_mode};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::{Cell, Element};
use crate::spatial::{EnumAll, GridState, Invert};
use crate::Mode::Insert;
use ruscii::app::{App, State};
use ruscii::drawing::{Pencil, RectCharset};
use ruscii::spatial::Vec2;
use ruscii::terminal::Color::{Cyan, Red, White};
use ruscii::terminal::{Color, Style, Window};
use Mode::Normal;

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
            Normal => {
                state.pos = "".to_string();
                state.sym = "".to_string();
                state.err = "".to_string();
                input_view_mode(app_state, &mut state)
            }
            Mode::Start => start_mode(app_state, &mut state),
        }

        let mut pencil = Pencil::new(window.canvas_mut());

        // Grid and startup screen
        pencil.set_origin(Vec2::xy(6, 0));
        match state.mode {
            Mode::Start => draw_logo(&mut pencil, &graph, Vec2::y(2)),
            _ => draw_grid(&mut pencil, &mut graph, Vec2::y(2), state.mode),
        }

        // Insert mode and cursor
        if let Insert = state.mode {
            draw_insert_mode(&mut pencil, Vec2::xy(graph.size.x * 3 / 2, 1));
            draw_cursor(&mut pencil, &graph, Vec2::y(2));
        } else {
            pencil.draw_rect(&RectCharset::simple_round_lines(), Vec2::xy(-1, 1), Vec2::xy(graph.size.x * 3 + 2, graph.size.y + 2));
            if let Mode::Start = state.mode {
                pencil.draw_text(&format!(" {version} "), Vec2::xy(graph.size.x * 3 - version.len() as i32 - 3, graph.size.y + 2));
            }
        }

        // Menu
        pencil
            .draw_center_text(
                match state.mode {
                    Mode::Start => "",
                    _ => " ChemCreator ",
                },
                Vec2::xy(graph.size.x * 3 / 2, if state.mode == Normal { 1 } else { 0 }),
            )
            .draw_text(
                &state.pos.to_string(),
                Vec2::xy(0, 1),
            )
            .draw_text(
                &state.sym.to_string(),
                Vec2::xy(graph.size.x * 3 - state.sym.len() as i32, 1),
            )
            .set_foreground(Red)
            .draw_text(&state.err.to_string(), Vec2::xy(graph.size.x * 3 + 3, 7))
            .draw_text(&state.debug.to_string(), Vec2::xy(graph.size.x * 3 + 3, 8))
            .set_foreground(White);

        if state.macros_enabled && state.mode == Insert {
            pencil
                .set_foreground(Yellow)
                .draw_text("Macro mode enabled.", Vec2::y(graph.size.y + 2))
                .set_foreground(White);
        }

        // Statistics
        pencil
            .move_origin(Vec2::xy(graph.size.x * 3 + 5, 1))
            .set_style(Style::Bold)
            .draw_text(
                &state.name.to_string(),
                Vec2::zero(),
            )
            .set_style(Style::Plain);

        if !matches!(state.mode, Normal) {
            return
        }

        let mut mass = 0.0;

        let mut missed = 0;
        for (index, element) in Element::all().into_iter().enumerate() {
            let count = graph.count(|it| it.is_atom() && it.unwrap_atom().element == element);

            if count == 0 {
                missed += 1;
                continue
            }

            mass += element.mass() * count as f32;

            pencil
                .set_foreground(element.color())
                .draw_text(element.symbol(), Vec2::y(2 + index - missed))
                .set_foreground(White)
                .draw_text(&format!("| {}", count), Vec2::xy(6, 2 + index - missed));
        }

        pencil.draw_text(&format!("atomic weight | {:.3} amu", mass), Vec2::xy(15, 2));
    });
}

fn draw_logo(pencil: &mut Pencil, graph: &GridState, pos: Vec2) {
    pencil
        .move_origin(pos)
        .draw_text("┌────┐          ", Vec2::xy(graph.size.x * 3 / 2 - 8, graph.size.y / 2 + 2).inv(graph),)
        .draw_text("│  C │hem       ", Vec2::xy(graph.size.x * 3 / 2 - 8, graph.size.y / 2 + 1).inv(graph),)
        .draw_text("└────┼────┐     ", Vec2::xy(graph.size.x * 3 / 2 - 8, graph.size.y / 2).inv(graph),)
        .draw_text("     │ Cr │eator", Vec2::xy(graph.size.x * 3 / 2 - 8, graph.size.y / 2 - 1).inv(graph),)
        .draw_text("     └────┘     ", Vec2::xy(graph.size.x * 3 / 2 - 8, graph.size.y / 2 - 2).inv(graph),)
        .move_origin(-pos);
}

fn draw_grid(pencil: &mut Pencil, graph: &mut GridState, pos: Vec2, mode: Mode) {
    pencil.move_origin(pos);
    for cell in graph.cells.iter().flatten() {
        pencil.set_foreground(cell.color()).draw_text(
            match &cell {
                Cell::Atom(it) => it.symbol(),
                Cell::Bond(it) => it.symbol(),
                Cell::None(_) => match mode {
                    Insert => " • ",
                    Normal => "   ",
                    _ => "   ",
                },
            },
            Vec2::xy(cell.pos().x * 3, cell.pos().y).inv(graph),
        );
    }
    pencil.move_origin(-pos);
}

fn draw_insert_mode(pencil: &mut Pencil, center: Vec2) {
    pencil.draw_center_text("-- INSERT MODE --", center);
}

fn draw_cursor(pencil: &mut Pencil, graph: &GridState, pos: Vec2) {
    let color = *pencil.foreground();
    pencil
        .move_origin(pos)
        .set_foreground(Cyan)
        .draw_text(
            "<",
            Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y).inv(graph),
        )
        .draw_text(
            ">",
            Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y).inv(graph),
        )
        .set_foreground(color)
        .move_origin(-pos);
}

/// Represents the mode the app is in at a given time.
#[derive(Copy, Clone, Default, Debug, PartialEq)]
enum Mode {
    Insert,
    Normal,
    #[default]
    Start,
}

/// Contains the running state of the app not including the grid.
#[derive(Default)]
struct AppState {
    mode: Mode,
    name: String,
    pos: String,
    sym: String,
    err: String,
    debug: String,
    macros_enabled: bool,
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
