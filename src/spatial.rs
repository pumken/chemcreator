//! # Spatial
//!
//! The `spatial` module provides functionality for the state and traversal of the grid with which
//! the user interacts, including the [`GridState`] struct.

use std::cmp::Ordering;
use ruscii::spatial::{Direction, Vec2};
use crate::molecule::Element::{C, H};
use crate::molecule::{Atom, Bond, BondOrder, Cell, Element};
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::pointer::Pointer;

/// Represents the state of a grid, including all its [Cell]s and the position of the cursor.
#[derive(Clone, Debug)]
pub struct GridState {
    pub(crate) cells: Vec<Vec<Cell>>,
    pub(crate) size: Vec2,
    pub(crate) cursor: Vec2,
}

impl GridState {
    /// Creates a new [GridState] with the given `width` and `height` populated by [Cell]s that
    /// all contain [Cell::None] and with the cursor set to the center.
    pub(crate) fn new(width: usize, height: usize) -> GridState {
        let mut cells = vec![];

        for x in 0usize..width {
            let mut row = vec![];
            for y in 0usize..height {
                let new_cell = Cell::None(Vec2::xy(x, y));
                row.push(new_cell);
            }
            cells.push(row);
        }

        GridState {
            cells,
            size: Vec2::xy(width, height),
            cursor: Vec2::xy(width / 2, height / 2),
        }
    }

    /// Returns a reference to the cell at the given `pos`.
    ///
    /// ## Errors
    ///
    /// If a `pos` is not a valid point within the graph, this function returns [Err].
    pub fn get(&self, pos: Vec2) -> Result<&Cell, String> {
        if !self.contains(pos) {
            return Err(format!("Invalid access of graph at ({}, {})", pos.x, pos.y));
        }
        Ok(&self.cells[pos.x as usize][pos.y as usize])
    }

    /// Sets the current [Cell] pointed to by the cursor to a [Cell::Atom] with the given [Element].
    pub fn put_atom(&mut self, element: Element) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_atom = Cell::Atom(Atom { element, pos: cursor_pos.to_vec2() });
        self.cells[cursor_pos.0][cursor_pos.1] = new_atom;
    }

    /// Sets the current [Cell] pointed to by the cursor to an [Cell::Bond] with the given
    /// [BondOrder].
    pub fn put_bond(&mut self, order: BondOrder) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_atom = Cell::Bond(Bond {
            pos: cursor_pos.to_vec2(),
            order,
            orient: if self.atom_adjacent() { Horiz } else { Vert },
        });
        self.cells[cursor_pos.0][cursor_pos.1] = new_atom;
    }


    pub fn upgrade_atom(&self, atom: Atom) -> Result<&Cell, String> {
        let cell = self.get(atom.pos)?;

        if !matches!(cell, Cell::Atom(_)) {
            return Err(format!("Attempted to upgrade non-atom cell at ({}, {})", atom.pos.x, atom.pos.y))
        }

        Ok(cell)
    }

    /// Sets the current [Cell] pointed to by the cursor to [Cell::None].
    pub(crate) fn clear(&mut self) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_none = Cell::None(cursor_pos.to_vec2());
        self.cells[cursor_pos.0][cursor_pos.1] = new_none;
    }

    pub fn is_empty(&self) -> bool {
        for cell in self.cells.iter().flatten() {
            match cell {
                Cell::None(_) => {}
                _ => return false
            }
        }
        true
    }

    /// Returns a [`Vec`] of all non-empty (i.e., not-[`Cell::None`]) [`Cell`]s.
    pub fn filled_cells(&self) -> Vec<&Cell> {
        self.find_all(|cell| match cell {
            Cell::Atom(_) | Cell::Bond(_) => true,
            Cell::None(_) => false
        })
    }

    /// Returns a reference to the [Cell] to which the cursor is currently pointing.
    pub(crate) fn current_cell(&self) -> &Cell {
        &self.cells[self.cursor.x as usize][self.cursor.y as usize]
    }

    /// Moves the cursor safely in the given `direction`.
    pub fn move_cursor(&mut self, direction: Direction) {
        let displacement = direction.to_vec2();
        if self.contains(self.cursor + displacement) {
            self.cursor += displacement;
        }
    }

    /// Returns whether the given position is a valid point on the grid or not.
    pub(crate) fn contains(&self, pos: Vec2) -> bool {
        pos.x >= 0 && pos.y >= 0 && pos.x < self.size.x && pos.y < self.size.y
    }

    /// Returns the number of [Cell]s that satisfy the given `predicate`.
    pub(crate) fn count(&self, predicate: fn(&Cell) -> bool) -> i32 {
        let mut count = 0;
        for cell in self.cells.iter().flatten() {
            if predicate(cell) {
                count += 1;
            }
        }
        count
    }

    /// Returns a reference to the first occurrence of a [`Cell`] that satisfies the given
    /// `predicate` or [`None`] if none are found that do.
    pub fn find(&self, predicate: fn(&Cell) -> bool) -> Option<&Cell> {
        self.cells.iter().flatten().find(|&cell| predicate(cell))
    }

    /// Returns a [`Vec`] of references to all [`Cells`] that satisfy the given `predicate`.
    pub fn find_all(&self, predicate: fn(&Cell) -> bool) -> Vec<&Cell> {
        let mut out = vec![];

        for cell in self.cells.iter().flatten() {
            if predicate(&cell) {
                out.push(cell);
            }
        }
        out
    }

    /// Determines if there are any [Atom]s that are horizontally adjacent to the current [Cell].
    fn atom_adjacent(&self) -> bool {
        let mut lptr = Pointer::new(self.current_cell(), self);
        let left = {
            lptr.move_ptr(Direction::Left);
            lptr.borrow()
        };
        let mut rptr = Pointer::new(self.current_cell(), self);
        let right = {
            rptr.move_ptr(Direction::Right);
            rptr.borrow()
        };

        if let Ok(Cell::Atom(_)) = left { return true; }
        if let Ok(Cell::Atom(_)) = right { return true; }
        false
    }

    /// A temporary function that identifies the molecule in the [GridState] by counting
    /// occurrences of the `[C]` and `[H]` [Atom].
    ///
    /// This should eventually be replaced because it ignores all bonds and `[O]` [Atom]s and
    /// assumes that the molecule must be an alkane, alkene, or alkyne.
    pub(crate) fn simple_counter(&self) -> Result<String, &str> {
        let carbon = self.count(|element| if let Cell::Atom(it) = element {
            it.element == C
        } else {
            false
        });
        let hydrogen = self.count(|element| if let Cell::Atom(it) = element {
            it.element == H
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

impl PartialEq for GridState {
    fn eq(&self, other: &Self) -> bool {
        if self.size != other.size {
            return false;
        }
        let a = self.cells.iter().flatten();
        let b = other.cells.iter().flatten().collect::<Vec<&Cell>>();

        for (i, cell) in a.enumerate() {
            if cell != b[i] {
                return false;
            }
        }
        true
    }
}

pub(crate) trait ToVec2 {
    fn to_vec2(self) -> Vec2;
}

impl ToVec2 for (usize, usize) {
    /// Converts the tuple to a [`Vec2`].
    fn to_vec2(self) -> Vec2 {
        Vec2::xy(self.0, self.1)
    }
}

impl ToVec2 for Direction {
    fn to_vec2(self) -> Vec2 {
        match self {
            Direction::Up => Vec2::y(1),
            Direction::Down => Vec2::y(-1),
            Direction::Right => Vec2::x(1),
            Direction::Left => Vec2::x(-1),
            Direction::None => Vec2::zero(),
        }
    }
}

pub(crate) trait Cellular {
    fn pos(&self) -> Vec2;
}

/// A trait made specifically for [`Vec2`] and [`Direction`] to invert the graph so that the origin
/// is in the bottom-left corner rather than the top left corner.
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

/// Creates a two-dimensional [`Vec`] (i.e., a `Vec<Vec<T>>`) with the given `x` and `y` with
/// every element set to `val`.
///
/// ## Examples
///
/// ```rust
/// let a = nested_vec![2; 4; true];
/// let mut b = vec![];
///
/// for x in 0usize..2usize {
///     b.push(vec![4; true]);
/// }
///
/// assert_eq!(a, b);
/// ```
#[macro_export]
macro_rules! nested_vec {
    ($x:expr; $y:expr; $val:expr) => {{
        let mut nested_vec = Vec::with_capacity($x as usize);
        for _ in 0usize..($x as usize) {
            nested_vec.push(vec![$val; $y as usize]);
        }
        nested_vec
    }};
}

pub(crate) trait EnumAll {
    fn all() -> Vec<Self> where Self: Sized;
}

impl EnumAll for Direction {
    /// Returns all [Direction]s except for [Direction::None].
    fn all() -> Vec<Self> {
        vec![Direction::Up, Direction::Down, Direction::Left, Direction::Right]
    }
}

pub(crate) trait FromVec2 {
    fn from_points(first: Vec2, second: Vec2) -> Result<Self, String> where Self: Sized;
}

impl FromVec2 for Direction {
    /// Returns the [`Direction`] from the `first` [`Vec2`] to the `second`.
    ///
    /// ## Errors
    ///
    /// If the two given points do not lie on an orthogonal line, an [`Err`] is returned.
    fn from_points(first: Vec2, second: Vec2) -> Result<Direction, String> {
        let displacement = second - first;

        if displacement == Vec2::zero() {
            return Err(format!("Passed equal points ({}, {}) to from_points.", first.x, first.y))
        }

        match (displacement.x.cmp(&0), displacement.y.cmp(&0)) {
            (Ordering::Less, Ordering::Equal) => Ok(Direction::Left),
            (Ordering::Greater, Ordering::Equal) => Ok(Direction::Right),
            (Ordering::Equal, Ordering::Less) => Ok(Direction::Down),
            (Ordering::Equal, Ordering::Greater) => Ok(Direction::Up),
            _ => Err(format!("Direction from ({}, {}) to ({}, {}) is not orthogonal.",
                             first.x, first.y, second.x, second.y)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::molecule::BondOrder::Single;
    use crate::molecule::Atom;
    use crate::molecule::Element::O;
    use crate::graph_with;
    use crate::test_utils::GW::{A, B};
    use super::*;

    #[test]
    fn gs_get_returns_err_out_of_bounds() {
        let graph = graph_with!(1, 1);
        let err = graph.get(Vec2::xy(-1, -1));

        assert!(matches!(err, Err(_)));
    }

    #[test]
    fn gs_new_creates_sized_grid() {
        let graph = GridState::new(20, 10);

        assert!(graph.contains(Vec2::xy(0, 0)));
        assert!(graph.contains(Vec2::xy(19, 9)));
        assert!(!graph.contains(Vec2::xy(20, 10)));
    }

    #[test]
    fn gs_get_returns_correct_cell() {
        let graph = graph_with!(2, 2,
            [0, 1; A(C)],
            [1, 0; A(O)],
            [1, 1; B(Single)]
        );

        assert_eq!(*graph.get(Vec2::xy(0, 0)).unwrap(), Cell::None(Vec2::xy(0, 0)));
        assert_eq!(*graph.get(Vec2::xy(1, 0)).unwrap(), Cell::Atom(Atom { element: O, pos: Vec2::xy(1, 0) }))
    }

    #[test]
    fn to_vec2_converts_correctly() {
        let x = (8usize, 2usize).to_vec2();
        assert_eq!(x, Vec2::xy(8, 2));
    }

    #[test]
    fn all_returns_every_direction() {
        let x = Direction::all();
        assert_eq!(x, vec![Direction::Up, Direction::Down, Direction::Left, Direction::Right])
    }

    #[test]
    fn from_points_returns_correct_direction() {
        let center = Vec2::zero();
        let up = Vec2::y(1);
        let right = Vec2::x(1);

        let a = Direction::from_points(center, up).unwrap();
        let b = Direction::from_points(center, right).unwrap();

        assert_eq!(a, Direction::Up);
        assert_eq!(b, Direction::Right)
    }
}
