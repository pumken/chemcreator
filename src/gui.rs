//! # GUI
//!
//! The `gui` module provides functions for drawing various UI elements.

use crate::molecule::BondOrder::{Double, Triple};
use crate::molecule::{Bond, BondOrientation, Cell, Element};
use crate::spatial::{EnumAll, GridState, InvertVec2};
use crate::Mode::{Insert, Normal};
use crate::{AppState, Mode};
use ruscii::drawing::{Pencil, RectCharset};
use ruscii::spatial::Vec2;
use ruscii::terminal::Color;
use ruscii::terminal::Color::{Cyan, DarkGrey, White};
use Color::Red;

pub(crate) fn draw_grid_box(pencil: &mut Pencil, graph_size: Vec2, pos: Vec2) {
    pencil.draw_rect(
        &RectCharset::simple_round_lines(),
        pos,
        Vec2::xy(graph_size.x * 3 + 2, graph_size.y + 2),
    );
}

pub(crate) fn draw_logo(pencil: &mut Pencil, graph_size: Vec2, pos: Vec2) {
    pencil
        .move_origin(pos)
        .draw_text(
            "┌────┐          ",
            Vec2::xy(graph_size.x * 3 / 2 - 8, graph_size.y / 2 + 2).inv(graph_size.y),
        )
        .draw_text(
            "│  C │hem       ",
            Vec2::xy(graph_size.x * 3 / 2 - 8, graph_size.y / 2 + 1).inv(graph_size.y),
        )
        .draw_text(
            "└────┼────┐     ",
            Vec2::xy(graph_size.x * 3 / 2 - 8, graph_size.y / 2).inv(graph_size.y),
        )
        .draw_text(
            "     │ Cr │eator",
            Vec2::xy(graph_size.x * 3 / 2 - 8, graph_size.y / 2 - 1).inv(graph_size.y),
        )
        .draw_text(
            "     └────┘     ",
            Vec2::xy(graph_size.x * 3 / 2 - 8, graph_size.y / 2 - 2).inv(graph_size.y),
        )
        .move_origin(-pos);
}

pub(crate) fn draw_grid(
    pencil: &mut Pencil,
    graph: &mut GridState,
    chain_highlighting: &Option<Vec<Vec2>>,
    pos: Vec2,
    mode: Mode,
) {
    pencil.move_origin(pos);
    let previous_color = pencil.foreground().to_owned();

    for cell in graph.cells.iter().flatten() {
        let color = if let Some(it) = chain_highlighting {
            if cell.pos() == it[0] {
                Cyan
            } else if it.contains(&cell.pos()) {
                Red
            } else {
                DarkGrey
            }
        } else {
            cell.color()
        };

        pencil.set_foreground(color).draw_text(
            match &cell {
                Cell::Atom(it) => it.symbol(),
                Cell::Bond(it) => it.symbol(),
                Cell::None(_) => match mode {
                    Insert => " • ",
                    Normal => "   ",
                    _ => "   ",
                },
            },
            Vec2::xy(cell.pos().x * 3, cell.pos().y).inv(graph.size.y),
        );
    }

    pencil.set_foreground(previous_color);
    pencil.move_origin(-pos);
}

pub fn draw_insert_mode(pencil: &mut Pencil, center: Vec2) {
    pencil.draw_center_text(" INSERT MODE ", center);
}

pub fn draw_cursor(pencil: &mut Pencil, graph: &GridState, pos: Vec2) {
    let color = *pencil.foreground();
    pencil
        .move_origin(pos)
        .set_foreground(Cyan)
        .draw_text(
            "<",
            Vec2::xy(graph.cursor.x * 3 - 1, graph.cursor.y).inv(graph.size.y),
        )
        .draw_text(
            ">",
            Vec2::xy(graph.cursor.x * 3 + 3, graph.cursor.y).inv(graph.size.y),
        )
        .set_foreground(color)
        .move_origin(-pos);
}

pub fn draw_start_message(pencil: &mut Pencil, pos: Vec2) {
    pencil
        .move_origin(pos)
        .draw_center_text("ChemCreator", Vec2::xy(20, 2))
        .draw_center_text("A text-based tool for", Vec2::xy(20, 4))
        .draw_center_text("identifying organic molecules.", Vec2::xy(20, 5))
        .draw_center_text("Written in Rust by Gavin Tran.", Vec2::xy(20, 7))
        .draw_center_text("Press any key to start.", Vec2::xy(20, 9))
        .move_origin(-pos);
}

pub(crate) fn draw_statistics(graph: &GridState, state: &mut AppState, pencil: &mut Pencil) {
    let mut mass = 0.0;

    let mut missed = 0;
    let mut final_index = 0;
    for (index, element) in Element::all().into_iter().enumerate() {
        let count = graph.count(|it| it.is_atom() && it.unwrap_atom().element == element);

        if count == 0 {
            missed += 1;
            continue;
        }

        mass += element.mass() * count as f32;

        pencil
            .set_foreground(element.color())
            .draw_text(element.symbol(), Vec2::y(2 + index - missed))
            .set_foreground(White)
            .draw_text(&format!("| {}", count), Vec2::xy(6, 2 + index - missed));

        final_index = 2 + index - missed;
    }

    let mut missed = 0;
    for (index, order) in [Double, Triple].into_iter().enumerate() {
        let count = graph.count(|it| it.is_bond() && it.unwrap_bond().order == order);

        if count == 0 {
            missed += 1;
            continue;
        }

        pencil
            .draw_text(
                Bond::new(Vec2::zero(), order, BondOrientation::Horiz).symbol(),
                Vec2::y(2 + final_index + index - missed),
            )
            .draw_text(
                &format!("| {}", count),
                Vec2::xy(6, 2 + final_index + index - missed),
            );
    }

    if state.err {
        return;
    }

    let carbon = graph.count(|it| it.is_atom() && it.unwrap_atom().element == Element::C);
    let nitrogen = graph.count(|it| it.is_atom() && it.unwrap_atom().element == Element::N);
    let hydrogen = graph.count(|it| it.is_atom() && it.unwrap_atom().element == Element::H);
    let halogens = graph.count(|it| {
        it.is_atom()
            && matches!(
                it.unwrap_atom().element,
                Element::Br | Element::Cl | Element::F | Element::I
            )
    });

    let ihd = (2 * carbon + 2 + nitrogen - hydrogen - halogens) / 2;

    pencil
        .draw_text(&format!("atomic weight | {:.3} amu", mass), Vec2::xy(15, 2))
        .draw_text(&format!("IHD           | {}", ihd), Vec2::xy(15, 3))
        .draw_text(
            &format!("name length   | {}", state.name.len()),
            Vec2::xy(15, 5),
        );
}

pub(crate) fn draw_wrapped_name(
    graph: &mut GridState,
    state: &mut AppState,
    window_size: Vec2,
    pencil: &mut Pencil,
) {
    let wrap_length = window_size.x - 14 - graph.size.x * 3;
    let lines = state
        .name
        .chars()
        .collect::<Vec<char>>()
        .chunks(wrap_length as usize)
        .map(|chunk| chunk.iter().collect())
        .collect::<Vec<String>>();
    pencil.set_foreground(if state.err { Red } else { White });

    for line in lines {
        pencil
            .draw_text(&line, Vec2::zero())
            .move_origin(Vec2::y(1));
    }
    pencil.move_origin(Vec2::y(-1)).set_foreground(White);
}
