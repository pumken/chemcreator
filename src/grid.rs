use ruscii::spatial::{Direction, Vec2};
use crate::algorithm::InvalidGraphError;
use crate::algorithm::InvalidGraphError::IncompleteBond;
use crate::molecule::Element::{C, H};
use crate::molecule::{Atom, Bond, BondOrder, BondOrientation, Cell, Element};
use crate::molecule::BondOrientation::{Horiz, Vert};

/// Represents the state of a grid, including all its [Cell]s and the position of the cursor.
#[derive(Clone)]
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
    pub(crate) fn put_atom(&mut self, element: Element) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_atom = Cell::Atom(Atom { element, pos: cursor_pos.to_vec2() });
        self.cells[cursor_pos.0][cursor_pos.1] = new_atom;
    }

    /// Sets the current [Cell] pointed to by the cursor to an [Cell::Bond] with the given
    /// [BondOrder].
    pub(crate) fn put_bond(&mut self, order: BondOrder) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_atom = Cell::Bond(Bond {
            pos: cursor_pos.to_vec2(),
            order,
            orient: if self.atom_adjacent() { Horiz } else { Vert },
        });
        self.cells[cursor_pos.0][cursor_pos.1] = new_atom;
    }

    /// Sets the current [Cell] pointed to by the cursor to [Cell::None].
    pub(crate) fn clear(&mut self) {
        let cursor_pos = (self.cursor.x as usize, self.cursor.y as usize);
        let new_none = Cell::None(cursor_pos.to_vec2());
        self.cells[cursor_pos.0][cursor_pos.1] = new_none;
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

    /// Returns whether the given position is a valid point on the grid or not.
    pub(crate) fn contains(&self, pos: Vec2) -> bool {
        pos.x >= 0 && pos.y >= 0 && pos.x < self.size.x && pos.y < self.size.y
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

    pub(crate) fn find_all(&self, predicate: fn(&Cell) -> bool) -> Vec<&Cell> {
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
        let left = &self.cells[self.cursor.x as usize - 1][self.cursor.y as usize];
        let right = &self.cells[self.cursor.x as usize + 1][self.cursor.y as usize];

        if let Cell::Atom(_) = left { return true; }
        if let Cell::Atom(_) = right { return true; }
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

trait ToVec2 {
    fn to_vec2(self) -> Vec2;
}

impl ToVec2 for (usize, usize) {
    fn to_vec2(self) -> Vec2 {
        Vec2::xy(self.0, self.1)
    }
}

/// A struct used to move around a [GridState] and borrow __immutable__ references to its cells.
///
/// Not to be confused with the pointer in computer science used to store a memory address.
///
/// > There are only two hard things in Computer Science: cache invalidation and naming things.
/// > — Phil Karlton
pub(crate) struct Pointer<'a> {
    pub(crate) graph: &'a GridState,
    pub(crate) pos: Vec2,
}

impl<'a> Pointer<'a> {
    pub(crate) fn new(cell: &Cell, graph: &'a GridState) -> Pointer<'a> {
        Pointer { graph, pos: cell.pos() }
    }

    pub(crate) fn borrow(&self) -> &Cell {
        &self.graph.cells[self.pos.x as usize][self.pos.y as usize]
    }

    /// Returns a [`Vec`] of references to the [`Cell`]s bonded to the atom currently pointed to.
    /// Assumes that the current cell is a [`Cell::Atom`].
    ///
    /// ## Panics
    ///
    /// When this function is called, the [`Pointer`] is expected to be pointing to a valid
    /// [`Cell::Atom`]. If it is found that this is not the case, this function will panic!
    ///
    /// ## Errors
    ///
    /// If one of the bonds to the current cell is found to be dangling, an
    /// [`IncompleteBond`] will be returned.
    pub(crate) fn bonded(&self) -> Result<Vec<&Cell>, InvalidGraphError> {
        match self.borrow() {
            Cell::Atom(_) => {}
            Cell::Bond(_) => panic!("Pointer cannot access bonded atoms: atom was expected at ({}, {}) but bond
            was found", self.pos.x, self.pos.y),
            Cell::None(_) => panic!("Pointer cannot access bonded atoms: atom was expected at ({}, {}) but none
            was found", self.pos.x, self.pos.y)
        };

        let mut out = vec![];

        for direction in Direction::all() {
            if let Ok(result) = self.graph.get(self.pos + direction.vec2()) {
                if let Cell::Bond(it) = result {
                    if it.orient == BondOrientation::from_direction(direction) {
                        out.push(self.traverse_bond(direction)?);
                    }
                }
            }
        }
        Ok(out)
    }

    pub(crate) fn bonded_carbons(&self) -> Result<usize, InvalidGraphError> {
        let bonded = self.bonded()?;
        Ok(bonded.iter()
            .filter(|atom| if let Cell::Atom(it) = atom {
                it.element == C
            } else {
                false
            })
            .map(|it| it)
            .count())
    }

    /// Gets the [`Cell::Atom`] at the other side of the bond in the given `direction`.
    ///
    /// Assumes that there indeed _is_ a bond in the given `direction`.
    ///
    /// ## Panics
    ///
    /// This function panics if [`Direction::None`] is passed as `direction`.
    ///
    /// ## Errors
    ///
    /// If the [`Pointer`] created to traverse the bond encounters a [`Cell::Bond`] of the incorrect
    /// orientation or a [`Cell::None`], an [`IncompleteBond`] is returned.
    fn traverse_bond(&self, direction: Direction) -> Result<&Cell, InvalidGraphError> {
        if let Direction::None = direction {
            panic!("Cannot pass Direction::None to traverse_bond.");
        }

        let mut traversal_ptr = Pointer { graph: self.graph, pos: self.pos };

        loop {
            traversal_ptr.move_ptr(direction);
            match traversal_ptr.borrow() {
                Cell::Atom(_) => break Ok(  // this doesn't work when using Pointer::borrow
                    &traversal_ptr.graph.cells[traversal_ptr.pos.x as usize][traversal_ptr.pos.y as usize]
                ),
                Cell::Bond(it) => if it.orient == BondOrientation::from_direction(direction) {
                    continue
                } else {
                    break Err(IncompleteBond(traversal_ptr.pos))
                }
                Cell::None(_) => break Err(IncompleteBond(traversal_ptr.pos))
            }
        }
    }

    pub fn move_ptr(&mut self, direction: Direction) -> bool {
        if self.graph.contains(self.pos + direction.vec2()) {
            self.pos += direction.vec2();
            true
        } else {
            false
        }
    }

    pub(crate) fn move_up(&mut self) -> bool {
        if self.pos.y < self.graph.size.y - 1 {
            self.pos.y += 1;
            true
        } else {
            false
        }
    }

    pub(crate) fn move_down(&mut self) -> bool {
        if self.pos.y > 0 {
            self.pos.y -= 1;
            true
        } else {
            false
        }
    }

    pub(crate) fn move_left(&mut self) -> bool {
        if self.pos.x > 0 {
            self.pos.y -= 1;
            true
        } else {
            false
        }
    }

    pub(crate) fn move_right(&mut self) -> bool {
        if self.pos.x < self.graph.size.x - 1 {
            self.pos.x += 1;
            true
        } else {
            false
        }
    }
}

pub(crate) trait Cellular {
    fn pos(&self) -> Vec2;
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

pub(crate) trait EnumAll {
    fn all() -> Vec<Self> where Self: Sized;
}

impl EnumAll for Direction {
    /// Returns all [Direction]s except for [Direction::None].
    fn all() -> Vec<Self> {
        vec![Direction::Up, Direction::Down, Direction::Left, Direction::Right]
    }
}