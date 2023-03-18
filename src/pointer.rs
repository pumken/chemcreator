//! # Pointer
//!
//! The `pointer` module contains the [`Pointer`] struct that allows for traversal over a
//! [`GridState`] and molecular structures on it.

use crate::groups::InvalidGraphError;
use crate::groups::InvalidGraphError::IncompleteBond;
use crate::molecule::Element::C;
use crate::molecule::{Atom, BondOrder, BondOrientation, Cell};
use crate::spatial::{EnumAll, GridState, ToVec2};
use ruscii::spatial::{Direction, Vec2};

/// A struct used to move around a [GridState] and borrow __immutable__ references to its cells.
///
/// Not to be confused with the pointer in computer science used to store a memory address.
///
/// > There are only two hard things in Computer Science: cache invalidation and naming things.
/// > — Phil Karlton
#[derive(Debug, Clone)]
pub(crate) struct Pointer<'a> {
    pub(crate) graph: &'a GridState,
    pub(crate) pos: Vec2,
}

impl<'a> Pointer<'a> {
    pub(crate) fn new(graph: &'a GridState, pos: Vec2) -> Pointer<'a> {
        Pointer { graph, pos }
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
            if let Ok(result) = self.graph.get(self.pos + direction.to_vec2()) {
                match result {
                    Cell::Atom(_) | Cell::Bond(_) => out.push(result),
                    Cell::None(_) => {}
                }
            }
        }

        out
    }

    /// Returns a [`Vec`] of [`Directions`] to the non-empty [`Cell`]s adjacent to the [`Cell`]
    /// currently pointed to.
    pub fn connected_directions(&self) -> Vec<Direction> {
        let mut out = vec![];

        for direction in Direction::all() {
            if let Ok(result) = self.graph.get(self.pos + direction.to_vec2()) {
                match result {
                    Cell::Atom(_) | Cell::Bond(_) => out.push(direction),
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
    pub(crate) fn bonded(&self) -> Result<Vec<Atom>, InvalidGraphError> {
        match self.borrow().unwrap() {
            Cell::Atom(_) => {}
            Cell::Bond(_) => panic!(
                "Pointer cannot access bonded atoms: atom was expected at ({}, {}) but bond
            was found",
                self.pos.x, self.pos.y
            ),
            Cell::None(_) => panic!(
                "Pointer cannot access bonded atoms: atom was expected at ({}, {}) but none
            was found",
                self.pos.x, self.pos.y
            ),
        };

        let mut out = vec![];

        for direction in Direction::all() {
            let adjacent_cell = self.graph.get(self.pos + direction.to_vec2());

            match adjacent_cell {
                Ok(Cell::Atom(_)) => out.push(self.traverse_bond(direction)?),
                Ok(Cell::Bond(it)) => {
                    if it.orient == BondOrientation::from_direction(direction) {
                        out.push(self.traverse_bond(direction)?)
                    }
                }
                _ => {}
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
        Ok(bonded
            .iter()
            .filter(|&atom| atom.element == C)
            .map(|it| it.to_owned())
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
        Ok(bonded.iter().filter(|&atom| atom.element == C).count())
    }

    /// Gets the [`Cell::Atom`] at the other side of the bond in the given `direction`.
    ///
    /// Assumes that there indeed _is_ a bond, implicit or not, in the given `direction`.
    ///
    /// ## Panics
    ///
    /// This function panics if [`Direction::None`] is passed as `direction`.
    ///
    /// ## Errors
    ///
    /// If the [`Pointer`] created to traverse the bond encounters a [`Cell::Bond`] of the incorrect
    /// orientation or a [`Cell::None`], an [`IncompleteBond`] is returned.
    pub fn traverse_bond(&self, direction: Direction) -> Result<Atom, InvalidGraphError> {
        if let Direction::None = direction {
            panic!("Cannot pass Direction::None to traverse_bond");
        }

        let mut traversal_ptr = Pointer {
            graph: self.graph,
            pos: self.pos,
        };

        loop {
            traversal_ptr.move_ptr(direction);
            match traversal_ptr.borrow() {
                Ok(Cell::Atom(it)) => break Ok(it.to_owned()),
                Ok(Cell::Bond(it)) => {
                    if it.orient == BondOrientation::from_direction(direction) {
                        continue;
                    } else {
                        break Err(IncompleteBond(traversal_ptr.pos));
                    }
                }
                Ok(Cell::None(_)) => break Err(IncompleteBond(traversal_ptr.pos)),
                Err(_) => break Err(IncompleteBond(traversal_ptr.pos)),
            }
        }
    }

    /// Returns the [`BondOrder`] of bond connected to the [`Cell::Atom`] currently pointed to in
    /// the given `direction`. If there is none or the cell is invalid, [`None`] is returned.
    ///
    /// ## Panics
    ///
    /// If the current cell is not a [`Cell::Atom`] or [`Direction::None`] is given, this function
    /// will panic.
    pub fn bond_order(&self, direction: Direction) -> Option<BondOrder> {
        let mut traversal_ptr = Pointer {
            graph: self.graph,
            pos: self.pos,
        };

        if let Direction::None = direction {
            panic!("Cannot pass Direction::None to bond_order")
        }

        if !traversal_ptr.move_ptr(direction) {
            return None;
        };
        match traversal_ptr.borrow().unwrap() {
            Cell::Atom(_) => Some(BondOrder::Single),
            Cell::Bond(it) => Some(it.order),
            _ => None,
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
            Cell::Bond(_) => {
                return Err(format!(
                    "Pointer cannot access bond count: atom was
            expected at ({}, {}) but bond was found",
                    self.pos.x, self.pos.y
                ))
            }
            Cell::None(_) => {
                return Err(format!(
                    "Pointer cannot access bond count: atom was
            expected at ({}, {}) but none was found",
                    self.pos.x, self.pos.y
                ))
            }
        };

        for direction in Direction::all() {
            if let Ok(cell) = self.graph.get(self.pos + direction.to_vec2()) {
                match cell {
                    Cell::Atom(_) => out += 1,
                    Cell::Bond(it) => out += it.order.order(),
                    Cell::None(_) => {}
                }
            }
        }
        Ok(out)
    }

    /// Moves the pointer safely in the given `direction`.
    ///
    /// If moving in the given `direction` would result in moving into an invalid point outside
    /// the [`GridState`], the move does not occur and `false` is returned.
    pub fn move_ptr(&mut self, direction: Direction) -> bool {
        if self.graph.contains(self.pos + direction.to_vec2()) {
            self.pos += direction.to_vec2();
            true
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::molecule::BondOrder::{Double, Single, Triple};
    use crate::molecule::Element::{C, H, O};
    use crate::test_utils::unwrap_atom;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn borrow_returns_correct_cell() {
        let graph = graph_with!(2, 2);
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let cell = ptr.borrow().unwrap();

        assert_eq!(cell, graph.get(Vec2::xy(1, 1)).unwrap())
    }

    #[test]
    fn borrow_out_of_bounds_returns_err() {
        let graph = graph_with!(1, 1);
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(-1, -1),
        };
        let cell = ptr.borrow();

        assert!(matches!(cell, Err(_)))
    }

    #[test]
    fn connected_returns_non_empty_neighbors() {
        let graph = graph_with!(3, 3,
            [1, 0; A(H)],
            [0, 1; A(H)],
            [1, 1; A(C)],
            [2, 1; B(Single)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let a = ptr.connected();
        let b = vec![
            graph.get(Vec2::xy(1, 0)).unwrap(),
            graph.get(Vec2::xy(0, 1)).unwrap(),
            graph.get(Vec2::xy(2, 1)).unwrap(),
        ];

        assert_eq!(a, b);
    }

    #[test]
    fn connected_behaves_at_boundaries() {
        let graph = graph_with!(2, 2,
            [0, 0; A(C)],
            [1, 0; A(H)],
            [0, 1; B(Triple)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let a = ptr.connected();
        let b = vec![
            graph.get(Vec2::xy(0, 1)).unwrap(),
            graph.get(Vec2::xy(1, 0)).unwrap(),
        ];

        assert_eq!(a, b);
    }

    #[test]
    fn connected_returns_directions() {
        let graph = graph_with!(4, 4,
            [0, 1; A(C)],
            [1, 0; A(H)], [1, 1; A(C)], [1, 2; B(Double)], [1, 3; A(O)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let directions = ptr.connected_directions();

        assert_eq!(
            directions,
            vec![Direction::Up, Direction::Down, Direction::Left]
        );
    }

    #[test]
    fn bonded() {
        let graph = graph_with!(5, 5,
            [2, 0; A(H)],
            [2, 1; B(Single)],
            [2, 2; A(C)],
            [3, 2; B(Single)],
            [4, 2; A(H)],
            [2, 3; B(Double)],
            [2, 4; A(O)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(2, 2),
        };
        let a = ptr.bonded().unwrap();
        let b = vec![
            graph.get(Vec2::xy(2, 4)).unwrap(), // ordered by direction u/d/l/r
            graph.get(Vec2::xy(2, 0)).unwrap(),
            graph.get(Vec2::xy(4, 2)).unwrap(),
        ]
        .iter()
        .map(|&cell| unwrap_atom(cell))
        .collect::<Vec<Atom>>();

        assert_eq!(a, b);
    }

    #[test]
    fn bonded_implicit() {
        let graph = graph_with!(3, 3,
            [1, 0; A(H)],
            [0, 1; A(H)],
            [1, 1; A(C)],
            [2, 1; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let a = ptr.bonded().unwrap();
        let b = vec![
            graph.get(Vec2::xy(1, 0)).unwrap(),
            graph.get(Vec2::xy(0, 1)).unwrap(),
            graph.get(Vec2::xy(2, 1)).unwrap(),
        ]
        .iter()
        .map(|&cell| unwrap_atom(cell))
        .collect::<Vec<Atom>>();

        assert_eq!(a, b);
    }

    #[test]
    fn bonded_behaves_at_boundaries() {
        let graph = graph_with!(3, 3,
            [0, 0; A(C)],
            [0, 1; A(H)],
            [1, 0; B(Single)],
            [2, 0; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let bonded = ptr.bonded().unwrap();
        let expected = vec![
            graph.get(Vec2::xy(0, 1)).unwrap(),
            graph.get(Vec2::xy(2, 0)).unwrap(),
        ]
        .iter()
        .map(|&cell| unwrap_atom(cell))
        .collect::<Vec<Atom>>();

        assert_eq!(bonded, expected);
    }

    #[test]
    fn bonded_carbons() {
        let graph = graph_with!(3, 3,
            [1, 0; A(H)],
            [0, 1; A(C)],
            [1, 1; A(C)],
            [2, 1; A(C)],
            [1, 2; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let a = ptr.bonded_carbons().unwrap();
        let b = vec![
            graph.get(Vec2::xy(0, 1)).unwrap(),
            graph.get(Vec2::xy(2, 1)).unwrap(),
        ]
        .iter()
        .map(|&cell| unwrap_atom(cell))
        .collect::<Vec<Atom>>();

        assert_eq!(a, b);
    }

    #[test]
    fn bonded_carbon_count() {
        let graph = graph_with!(3, 3,
            [1, 0; A(C)],
            [0, 1; A(C)],
            [1, 1; A(C)],
            [2, 1; A(C)],
            [1, 2; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::xy(1, 1),
        };
        let carbons = ptr.bonded_carbon_count().unwrap();

        assert_eq!(carbons, 3usize);
    }

    #[test]
    fn traverse_bond() {
        let graph = graph_with!(1, 3,
            [0, 0; A(C)],
            [0, 1; B(Single)],
            [0, 2; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let atom = ptr.traverse_bond(Direction::Up).unwrap();

        assert_eq!(
            atom,
            Atom {
                element: H,
                pos: Vec2::y(2)
            }
        );
    }

    #[test]
    fn traverse_bond_variable_length() {
        let graph = graph_with!(1, 10,
            [0, 0; A(C)],
            [0, 1; B(Single)],
            [0, 2; B(Single)],
            [0, 3; B(Single)],
            [0, 4; B(Single)],
            [0, 5; B(Single)],
            [0, 6; B(Single)],
            [0, 7; B(Single)],
            [0, 8; B(Single)],
            [0, 9; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let atom = ptr.traverse_bond(Direction::Up).unwrap();

        assert_eq!(
            atom,
            Atom {
                element: H,
                pos: Vec2::xy(0, 9)
            }
        );
    }

    #[test]
    fn traverse_bond_implicit() {
        let graph = graph_with!(1, 2,
            [0, 0; A(C)],
            [0, 1; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let atom = ptr.traverse_bond(Direction::Up).unwrap();

        assert_eq!(
            atom,
            Atom {
                element: H,
                pos: Vec2::y(1)
            }
        );
    }

    #[test]
    fn traverse_bond_out_of_bounds_returns_err() {
        let graph = graph_with!(3, 1,
            [0, 0; A(C)],
            [1, 0; B(Single)],
            [2, 0; B(Single)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let err = ptr.traverse_bond(Direction::Right);

        assert_eq!(err, Err(IncompleteBond(Vec2::x(2))));
    }

    #[test]
    fn bond_order() {
        let graph = graph_with!(2, 2,
            [0, 0; A(C)],
            [0, 1; B(Triple)],
            [1, 0; B(Double)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let a = ptr.bond_order(Direction::Up).unwrap();
        let b = ptr.bond_order(Direction::Right).unwrap();

        assert_eq!(a, Triple);
        assert_eq!(b, Double);
    }

    #[test]
    fn bond_order_implicit() {
        let graph = graph_with!(2, 2,
            [0, 0; A(C)],
            [0, 1; A(O)],
            [1, 0; A(H)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let a = ptr.bond_order(Direction::Up).unwrap();
        let b = ptr.bond_order(Direction::Right).unwrap();
        let c = ptr.bond_order(Direction::Left);

        assert_eq!(a, Single);
        assert_eq!(b, Single);
        assert_eq!(c, None);
    }

    #[test]
    fn bond_order_out_of_bounds() {
        let graph = graph_with!(1, 1,
            [0, 0; A(C)]
        );
        let ptr = Pointer {
            graph: &graph,
            pos: Vec2::zero(),
        };
        let option = ptr.bond_order(Direction::Right);

        assert_eq!(option, None)
    }
}
