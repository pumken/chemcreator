use std::fmt::{Display, Formatter};
use ruscii::spatial::{Direction, Vec2};
use ruscii::terminal::Color;
use ruscii::terminal::Color::{LightGrey, Red, White};
use crate::grid::Cellular;
use crate::molecule::Element::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

/// Represents a functional group, containing the [Cell]s that comprise it.
pub(crate) struct Group {
    pub(crate) cells: Vec<Cell>,
    pub(crate) class: GroupType,
}

/// Represents a type of functional group on a molecule.
pub(crate) enum GroupType {
    /* Alkyl groups */
    Methyl,
    Ethyl,
    Propyl,
    Isopropyl,
    Butyl,
    Pentyl,
    Hexyl,
    Heptyl,
    Octyl,
    Nonyl,
    Decyl,
    /* Alkenyl groups in future */
    /* Alkynyl groups in future */
    /* Halide groups */
    Bromo,
    Chloro,
    Fluoro,
    Iodo,
    /* General groups */
    Hydroxyl,
    Carbonyl,
    Carboxyl,
    /* Phenyl later */ Ester,
    Ether,
}

/// Represents a molecular component.
#[derive(Clone, Debug)]
pub enum Cell {
    Atom(Atom),
    Bond(Bond),
    None(Vec2),
}

impl Cell {
    pub(crate) fn color(&self) -> Color {
        match &self {
            Cell::Atom(it) => match it.element {
                C => LightGrey,
                O => Red,
                _ => White
            },
            _ => White
        }
    }

    pub(crate) fn pos(&self) -> Vec2 {
        match self {
            Cell::Atom(it) => it.pos,
            Cell::Bond(it) => it.pos,
            Cell::None(it) => *it
        }
    }

    pub(crate) fn is_not_atom(&self) -> bool {
        if let Cell::Atom(_) = self {
            false
        } else {
            true
        }
    }
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub(crate) element: Element,
    pub(crate) pos: Vec2
}

impl Atom {
    pub(crate) const fn symbol(&self) -> &str {
        match self.element {
            C => "[C]",
            H => "[H]",
            O => "[O]"
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Element {
    C,
    H,
    O,
}

impl Element {
    pub(crate) const fn bond_number(&self) -> i32 {
        match *self {
            C => 4,
            H => 1,
            O => 2
        }
    }
}

impl Display for Element {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match &self {
            C => "Carbon",
            H => "Hydrogen",
            O => "Oxygen"
        })
    }
}

#[derive(Copy, Debug)]
pub struct Bond {
    pub(crate) pos: Vec2,
    pub(crate) order: BondOrder,
    pub(crate) orient: BondOrientation,
}

impl Bond {
    pub(crate) fn symbol(&self) -> &str {
        match (&self.order, &self.orient) {
            (Single, Horiz) => "———",
            (Single, Vert) => " | ",
            (Double, Horiz) => "===",
            (Double, Vert) => " ‖ ",
            (Triple, Horiz) => "≡≡≡",
            (Triple, Vert) => " T ",
        }
    }
}

impl Cellular for Bond {
    fn pos(&self) -> Vec2 {
        self.pos
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BondOrientation {
    Vert,
    Horiz,
}

impl BondOrientation {
    /// Returns the related [`BondOrientation`] to the given `direction`.
    ///
    /// ## Panics
    ///
    /// If [`Direction::None`] is passed to this function.
    pub const fn from_direction(direction: Direction) -> BondOrientation {
        match direction {
            Direction::Up | Direction::Down => Vert,
            Direction::Left | Direction::Right => Horiz,
            Direction::None => panic!("Attempted to pass Direction::None to from_direction.")
        }
    }
}

impl Clone for Bond {
    fn clone(&self) -> Bond {
        Bond { pos: self.pos.clone(), order: self.order.clone(), orient: self.orient.clone() }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}