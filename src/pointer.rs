use ruscii::spatial::{Direction, Vec2};
use crate::algorithm::InvalidGraphError;
use crate::algorithm::InvalidGraphError::IncompleteBond;
use crate::molecule::{Atom, BondOrientation, Cell};
use crate::molecule::Element::C;
use crate::spatial::{EnumAll, GridState};

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

    /// Returns a reference to the cell currently pointed to.
    ///
    /// ## Error
    ///
    /// If this [`Pointer`] does not point to a valid cell, an [`Err`] is returned.
    pub(crate) fn borrow(&self) -> Result<&Cell, String> {
        self.graph.get(self.pos)
    }

    /// Returns a [`Vec`] of references to the non-empty [`Cell`]s adjacent to the [`Cell`]
    /// currently pointed to.
    pub fn connected(&self) -> Vec<&Cell> {
        let mut out = vec![];

        for direction in Direction::all() {
            if let Ok(result) = self.graph.get(self.pos + direction.vec2()) {
                match result {
                    Cell::Atom(_) | Cell::Bond(_) => out.push(result),
                    Cell::None(_) => {}
                }
            }
        }

        out
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
        match self.borrow().unwrap() {
            Cell::Atom(_) => {}
            Cell::Bond(_) => panic!("Pointer cannot access bonded atoms: atom was expected at ({}, {}) but bond
            was found", self.pos.x, self.pos.y),
            Cell::None(_) => panic!("Pointer cannot access bonded atoms: atom was expected at ({}, {}) but none
            was found", self.pos.x, self.pos.y)
        };

        let mut out = vec![];

        for direction in Direction::all() {
            if let Ok(Cell::Bond(it)) = self.graph.get(self.pos + direction.vec2()) {
                if it.orient == BondOrientation::from_direction(direction) {
                    out.push(self.traverse_bond(direction)?);
                }
            }
        }
        Ok(out)
    }

    /// Returns a [`Vec`] of the carbon atoms bonded to the current cell.
    ///
    /// ## Panics
    ///
    /// This function will panic if the current cell is not a [`Cell::Atom`].
    ///
    /// ## Errors
    ///
    /// If one of the bonds to the current cell is found to be dangling, an
    /// [`IncompleteBond`] will be returned.
    pub fn bonded_carbons(&self) -> Result<Vec<Atom>, InvalidGraphError> {
        let bonded = self.bonded()?;
        Ok(bonded.iter()
            .filter(|&&cell| if let Cell::Atom(it) = cell {
                it.element == C
            } else {
                false
            })
            .map(|&it| match it {
                Cell::Atom(it) => it.to_owned(),
                _ => panic!("bonded returned non-atom cell to bonded_carbons")
            })
            .collect::<Vec<Atom>>())
    }

    /// Returns the number of carbon atoms bonded to the current cell.
    ///
    /// ## Panics
    ///
    /// This function will panic if the current cell is not a [`Cell::Atom`].
    ///
    /// ## Errors
    ///
    /// If one of the bonds to the current cell is found to be dangling, an
    /// [`IncompleteBond`] will be returned.
    pub fn bonded_carbon_count(&self) -> Result<usize, InvalidGraphError> {
        let bonded = self.bonded()?;
        Ok(bonded.iter()
            .filter(|&&atom| if let Cell::Atom(it) = atom {
                it.element == C
            } else {
                false
            })
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
                Ok(Cell::Atom(_)) => break Ok(  // this doesn't work when using Pointer::borrow
                                            &traversal_ptr.graph.cells[traversal_ptr.pos.x as usize][traversal_ptr.pos.y as usize]
                ),
                Ok(Cell::Bond(it)) => if it.orient == BondOrientation::from_direction(direction) {
                    continue;
                } else {
                    break Err(IncompleteBond(traversal_ptr.pos));
                }
                Ok(Cell::None(_)) => break Err(IncompleteBond(traversal_ptr.pos)),
                Err(_) => break Err(IncompleteBond(traversal_ptr.pos))
            }
        }
    }

    /// Gets the amount of bonds that the current [`Cell::Atom`] has.
    ///
    /// Assumes that the current cell is a [`Cell::Atom`]. Counts attached atoms as single
    /// bonds.
    ///
    /// ## Panics
    ///
    /// If the current cell is not valid, this function will panic.
    ///
    /// ## Errors
    ///
    /// If the current cell is not a [`Cell::Atom`] a message describing the error and the position
    /// at which it occurred is returned.
    pub(crate) fn bond_count(&self) -> Result<u8, String> {
        let mut out = 0;

        match self.borrow().unwrap() {
            Cell::Atom(_) => {}
            Cell::Bond(_) => return Err(format!("Pointer cannot access bond count: atom was
            expected at ({}, {}) but bond was found", self.pos.x, self.pos.y)),
            Cell::None(_) => return Err(format!("Pointer cannot access bond count: atom was
            expected at ({}, {}) but none was found", self.pos.x, self.pos.y)),
        };

        for direction in Direction::all() {
            if let Ok(cell) = self.graph.get(self.pos + direction.vec2()) {
                match cell {
                    Cell::Atom(_) => out += 1,
                    Cell::Bond(it) => out += it.order.order(),
                    Cell::None(_) => {}
                }
            }
        }
        Ok(out)
    }

    pub fn move_ptr(&mut self, direction: Direction) -> bool {
        if self.graph.contains(self.pos + direction.vec2()) {
            self.pos += direction.vec2();
            true
        } else {
            false
        }
    }
}
