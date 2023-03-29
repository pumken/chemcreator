//! # GUI
//!
//! The `gui` module provides functions for drawing various UI elements.

use Color::Red;
use crate::molecule::Cell;
use crate::spatial::{GridState, Invert};
use crate::Mode;
use crate::Mode::{Insert, Normal};
use ruscii::drawing::{Pencil, RectCharset};
use ruscii::spatial::Vec2;
use ruscii::terminal::Color;
use ruscii::terminal::Color::Cyan;

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
    mode: Mode
) {
    pencil.move_origin(pos);

    for cell in graph.cells.iter().flatten() {
        let color = if let Some(it) = chain_highlighting {
            if cell.pos() == it[0] {
                Cyan
            } else if it.contains(&cell.pos()) {
                Red
            } else {
                cell.color()
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
