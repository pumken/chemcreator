use std::fmt::Debug;
use ruscii::spatial::Vec2;
use crate::molecule::Atom::{C, H};
use crate::molecule::Symbol;

/// Represents the state of a grid, including all its [Cell]s and the position of the cursor.
#[derive(Clone)]
pub(crate) struct GridState {
    pub(crate) cells: Vec<Vec<Cell>>,
    pub(crate) size: Vec2,
    pub(crate) cursor: Vec2,
}

impl GridState {
    /// Creates a new [GridState] with the given `width` and `height` populated by [Cell]s that
    /// all contain [Symbol::None] and with the cursor set to position (0, 0).
    pub(crate) fn new(width: usize, height: usize) -> GridState {
        let mut cells = vec![];

        for x in 0usize..width {
            let mut row = vec![];
            for y in 0usize..height {
                let new_cell = Cell { sym: Symbol::None, pos: Vec2::xy(x, y) };
                row.push(new_cell);
            }
            cells.push(row);
        }

        GridState {
            cells,
            size: Vec2::xy(width, height),
            cursor: Vec2::xy(0, 0),
        }
    }

    /// Sets the [Symbol] that the current cell contains to the given [Symbol].
    pub(crate) fn insert(&mut self, symbol: Symbol) {
        self.current_cell_mut().sym = symbol;
    }

    /// Returns a reference to the [Cell] to which the cursor is currently pointing.
    pub(crate) fn current_cell(&self) -> &Cell {
        &self.cells[self.cursor.x as usize][self.cursor.y as usize]
    }

    /// Returns a mutable reference to the [Cell] to which the cursor is currently pointing.
    pub(crate) fn current_cell_mut(&mut self) -> &mut Cell {
        &mut self.cells[self.cursor.x as usize][self.cursor.y as usize]
    }

    /// Moves the cursor upwards safely.
    pub(crate) fn move_up(&mut self) {
        if self.cursor.y != self.size.y - 1 {
            self.cursor.y += 1;
        }
    }

    /// Moves the cursor downwards safely.
    pub(crate) fn move_down(&mut self) {
        if self.cursor.y != 0 {
            self.cursor.y -= 1;
        }
    }

    /// Moves the cursor leftwards safely.
    pub(crate) fn move_left(&mut self) {
        if self.cursor.x != 0 {
            self.cursor.x -= 1;
        }
    }

    /// Moves the cursor rightwards safely.
    pub(crate) fn move_right(&mut self) {
        if self.cursor.x != self.size.x - 1 {
            self.cursor.x += 1;
        }
    }

    /// Returns the number of [Cell]s that satisfy the given `predicate`.
    pub(crate) fn count(&self, predicate: fn(&Cell) -> bool) -> i32 {
        let mut count = 0;
        for cell in self.cells.iter().flatten() {
            if predicate(&cell) {
                count += 1;
            }
        }
        count
    }

    /// Determines if there are any [Atom]s that are horizontally adjacent to the current [Cell].
    pub(crate) fn atom_adjacent(&self) -> bool {
        let left = &self.cells[self.cursor.x as usize - 1][self.cursor.y as usize];
        let right = &self.cells[self.cursor.x as usize + 1][self.cursor.y as usize];

        if let Symbol::Atom(_) = left.sym { return true }
        if let Symbol::Atom(_) = right.sym { return true }
        false
    }

    /// A temporary function that identifies the molecule in the [GridState] by counting
    /// occurrences of the `[C]` and `[H]` [Atom].
    ///
    /// This should eventually be replaced because it ignores all bonds and `[O]` [Atom]s and
    /// assumes that the molecule must be an alkane, alkene, or alkyne.
    pub(crate) fn simple_counter(&self) -> Result<String, &str> {
        let carbon = self.count(|element| if let Symbol::Atom(it) = element.sym {
            it == C
        } else {
            false
        });
        let hydrogen = self.count(|element| if let Symbol::Atom(it) = element.sym {
            it == H
        } else {
            false
        });
        let prefix = match carbon {
            0 => return Err("No carbons found."),
            1 => "meth",
            2 => "eth",
            3 => "prop",
            4 => "but",
            5 => "pent",
            6 => "hex",
            7 => "hept",
            8 => "oct",
            9 => "non",
            10 => "dec",
            _ => return Err("Cannot exceed carbon length of 10.")
        };
        let suffix = if hydrogen == 2 * carbon - 2 {
            "yne"
        } else if hydrogen == 2 * carbon {
            "ene"
        } else if hydrogen == 2 * carbon + 2 {
            "ane"
        } else {
            return Err("Cannot determine class other than alkane, alkene, or alkyne.");
        };
        Ok(format!("{}{}", prefix, suffix))
    }
}

/// Represents a cell in the [GridState].
#[derive(Clone, Debug)]
pub(crate) struct Cell {
    pub(crate) sym: Symbol,
    pub(crate) pos: Vec2,
}

/// A trait made specifically for [Vec2] to invert the graph so that the origin is in the
/// bottom-left corner rather than the top left corner.
///
/// The sole function returns a [Vec2] because it is not expected that any other type will use it.
pub(crate) trait Invert {
    fn inv(self, graph: &GridState) -> Vec2;
}

impl Invert for Vec2 {
    /// Inverts the `y` of this [Vec2]. For *why* this exists, see [Invert].
    fn inv(mut self, graph: &GridState) -> Vec2 {
        self.y = graph.size.y - self.y - 1;
        self
    }
}