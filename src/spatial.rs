//! # Spatial
//!
//! The `spatial` module provides functionality for the state and traversal of the grid with which
//! the user interacts, including the [`GridState`] struct.

use crate::macros::bond_conversion;
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::molecule::{Atom, Bond, BondOrder, Cell, ComponentType, Element};
use crate::pointer::Pointer;
use ruscii::spatial::{Direction, Vec2};
use std::cmp::Ordering;

/// Represents the state of a grid, including all its [Cell]s and the position of the cursor.
#[derive(Clone, Debug)]
pub struct GridState {
    pub cells: Vec<Vec<Cell>>,
    pub size: Vec2,
    pub cursor: Vec2,
}

impl GridState {
    /// Creates a new [`GridState`] with the given `width` and `height` populated by [`Cell`]s that
    /// all contain [`Cell::None`] and with the cursor set to the center.
    pub fn new(width: usize, height: usize) -> GridState {
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

    /// Returns a reference to the [`Cell`] at the given `pos`.
    ///
    /// ## Errors
    ///
    /// If the `pos` is not a valid point within the graph, this function returns [`Err`].
    pub fn get(&self, pos: Vec2) -> Result<&Cell, String> {
        if !self.contains(pos) {
            return Err(format!("Invalid access of graph at ({}, {})", pos.x, pos.y));
        }
        Ok(&self.cells[pos.x as usize][pos.y as usize])
    }

    /// Returns a reference to the [`Cell`] at the given `pos`.
    ///
    /// ## Errors
    ///
    /// If the `pos` is not a valid point within the graph, this function returns [`Err`].
    pub(crate) fn get_mut(&mut self, pos: Vec2) -> Result<&mut Cell, String> {
        if !self.contains(pos) {
            return Err(format!("Invalid access of graph at ({}, {})", pos.x, pos.y));
        }
        Ok(&mut self.cells[pos.x as usize][pos.y as usize])
    }

    /// Returns a reference to the [`Cell`] to which the cursor is currently pointing.
    ///
    /// ## Errors
    ///
    /// If the `pos` is not a valid point within the graph, this function returns [`Err`].
    pub fn current_cell(&self) -> Result<&Cell, String> {
        self.get(self.cursor)
    }

    /// Returns a mutable reference to the [`Cell`] to which the cursor is currently pointing.
    ///
    /// ## Errors
    ///
    /// If the `pos` is not a valid point within the graph, this function returns [`Err`].
    fn current_cell_mut(&mut self) -> Result<&mut Cell, String> {
        self.get_mut(self.cursor)
    }

    /// Sets the current [`Cell`] pointed to by the cursor to a [`Cell::Atom`] with the given
    /// [`Element`].
    pub fn put_atom(&mut self, element: Element) {
        let atom = Cell::Atom(Atom {
            element,
            pos: self.cursor,
        });
        *self
            .current_cell_mut()
            .expect("cursor should be within bounds") = atom;
    }

    /// Sets the current [`Cell`] pointed to by the cursor to an [`Cell::Bond`] with the given
    /// [`BondOrder`].
    pub fn put_bond(&mut self, order: BondOrder) {
        let orient = if self.atom_adjacent(self.cursor) {
            Horiz
        } else {
            Vert
        };
        let bond = Cell::Bond(Bond {
            pos: self.cursor,
            order,
            orient,
        });
        *self
            .current_cell_mut()
            .expect("cursor should be within bounds") = bond;
        bond_conversion(self, self.cursor, order, orient);
    }

    /// Sets the current [`Cell`] pointed to by the cursor to [`Cell::None`].
    pub fn clear_cell(&mut self) {
        let empty_cell = Cell::None(self.cursor);
        *self
            .current_cell_mut()
            .expect("cursor should be within bounds") = empty_cell;
    }

    /// Sets the [`Cell`] at the given `pos` to the given `comp`.
    pub fn put(&mut self, pos: Vec2, comp: ComponentType) {
        *self.get_mut(pos).unwrap() = match comp {
            ComponentType::Element(it) => Cell::Atom(Atom { element: it, pos }),
            ComponentType::Order(it) => {
                let orient = if self.atom_adjacent(pos) { Horiz } else { Vert };
                bond_conversion(self, self.cursor, it, orient);
                Cell::Bond(Bond {
                    pos,
                    order: it,
                    orient,
                })
            }
            ComponentType::None => Cell::None(pos),
        }
    }

    /// Empties every [`Cell`] on the [`GridState`].
    pub fn clear_all(&mut self) {
        for cell in self.cells.iter_mut().flatten() {
            *cell = Cell::None(cell.pos());
        }
    }

    /// Checks if the [`GridState`] is empty, i.e., all the cells it contains are set to
    /// [`Cell::None`].
    pub fn is_empty(&self) -> bool {
        for cell in self.cells.iter().flatten() {
            match cell {
                Cell::None(_) => {}
                _ => return false,
            }
        }
        true
    }

    /// Returns a [`Vec`] of all non-empty (i.e., not-[`Cell::None`]) [`Cell`]s.
    pub fn filled_cells(&self) -> Vec<&Cell> {
        self.find_all(|cell| match cell {
            Cell::Atom(_) | Cell::Bond(_) => true,
            Cell::None(_) => false,
        })
    }

    /// Moves the cursor safely in the given `direction`.
    pub fn move_cursor(&mut self, direction: Direction) {
        let displacement = direction.to_vec2();
        if self.contains(self.cursor + displacement) {
            self.cursor += displacement;
        }
    }

    /// Returns whether the given `pos` is a valid point on the grid or not.
    pub fn contains(&self, pos: Vec2) -> bool {
        pos.x >= 0 && pos.y >= 0 && pos.x < self.size.x && pos.y < self.size.y
    }

    /// Returns a reference to the first occurrence of a [`Cell`] that satisfies the given
    /// `predicate` or [`None`] if none are found that do.
    pub fn find(&self, predicate: fn(&Cell) -> bool) -> Option<&Cell> {
        self.cells.iter().flatten().find(|&cell| predicate(cell))
    }

    /// Returns a [`Vec`] of references to all [`Cell`]s that satisfy the given `predicate`.
    pub fn find_all(&self, predicate: fn(&Cell) -> bool) -> Vec<&Cell> {
        self.cells
            .iter()
            .flatten()
            .filter(|&cell| predicate(cell))
            .collect::<Vec<&Cell>>()
    }

    pub fn count<F>(&self, predicate: F) -> i32
    where
        F: Fn(&Cell) -> bool,
    {
        self.cells
            .iter()
            .flatten()
            .filter(|&cell| predicate(cell))
            .count() as i32
    }

    /// Determines if there are any [`Atom`]s that are horizontally adjacent to the current
    /// [`Cell`].
    fn atom_adjacent(&self, pos: Vec2) -> bool {
        let mut lptr = Pointer::new(self, pos);
        let left = {
            lptr.move_ptr(Direction::Left);
            lptr.borrow()
        };
        let mut rptr = Pointer::new(self, pos);
        let right = {
            rptr.move_ptr(Direction::Right);
            rptr.borrow()
        };

        GridState::is_atom_or_horizontal(left) || GridState::is_atom_or_horizontal(right)
    }

    fn is_atom_or_horizontal(cell_result: Result<&Cell, String>) -> bool {
        if let Ok(cell) = cell_result {
            cell.is_atom()
                || if let Cell::Bond(it) = cell {
                    it.orient == Horiz
                } else {
                    false
                }
        } else {
            false
        }
    }
}

// Implemented to ignore cursor
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
    fn all() -> Vec<Self>
    where
        Self: Sized;
}

impl EnumAll for Direction {
    /// Returns all [Direction]s except for [Direction::None].
    fn all() -> Vec<Self> {
        vec![
            Direction::Up,
            Direction::Down,
            Direction::Left,
            Direction::Right,
        ]
    }
}

pub(crate) trait FromVec2 {
    fn from_points(first: Vec2, second: Vec2) -> Result<Self, String>
    where
        Self: Sized;
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
            return Err(format!(
                "Passed equal points ({}, {}) to from_points.",
                first.x, first.y
            ));
        }

        match (displacement.x.cmp(&0), displacement.y.cmp(&0)) {
            (Ordering::Less, Ordering::Equal) => Ok(Direction::Left),
            (Ordering::Greater, Ordering::Equal) => Ok(Direction::Right),
            (Ordering::Equal, Ordering::Less) => Ok(Direction::Down),
            (Ordering::Equal, Ordering::Greater) => Ok(Direction::Up),
            _ => Err(format!(
                "Direction from ({}, {}) to ({}, {}) is not orthogonal.",
                first.x, first.y, second.x, second.y
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::Atom;
    use crate::molecule::BondOrder::Single;
    use crate::molecule::Element::{C, O};
    use crate::test_utils::GW::{A, B};

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
            [1, 1; B(Single)],
        );

        assert_eq!(
            *graph.get(Vec2::xy(0, 0)).unwrap(),
            Cell::None(Vec2::xy(0, 0))
        );
        assert_eq!(
            *graph.get(Vec2::xy(1, 0)).unwrap(),
            Cell::Atom(Atom {
                element: O,
                pos: Vec2::xy(1, 0),
            })
        )
    }

    #[test]
    fn all_returns_every_direction() {
        let x = Direction::all();
        assert_eq!(
            x,
            vec![
                Direction::Up,
                Direction::Down,
                Direction::Left,
                Direction::Right,
            ]
        )
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
