use std::fmt::{Display, Formatter};
use ruscii::terminal::Color;
use ruscii::terminal::Color::{LightGrey, Red, White};
use crate::grid::{Cell, GridState};
use crate::molecule::Atom::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

/// Represents the current molecule.
pub(crate) struct Molecule {
    pub(crate) name: String,
    pub(crate) groups: Vec<Group>,
}

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
pub(crate) enum Symbol {
    Atom(Atom),
    Bond(Bond),
    None,
}

impl Symbol {
    pub(crate) fn color(&self) -> Color {
        match &self {
            Symbol::Atom(it) => match it {
                C => LightGrey,
                O => Red,
                _ => White
            },
            _ => White
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Atom {
    C,
    H,
    O,
}

impl Atom {
    pub(crate) const fn symbol(&self) -> &str {
        match *self {
            C => "[C]",
            H => "[H]",
            O => "[O]"
        }
    }

    pub(crate) const fn bond_number(&self) -> i32 {
        match *self {
            C => 4,
            H => 1,
            O => 2
        }
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match &self {
            C => "Carbon",
            H => "Hydrogen",
            O => "Oxygen"
        })
    }
}

#[derive(Copy, Debug)]
pub(crate) struct Bond {
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

    pub(crate) fn adjusted(order: BondOrder, graph: &GridState) -> Bond {
        Bond {
            order,
            orient: if graph.atom_adjacent() {
                Horiz
            } else {
                Vert
            },
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum BondOrientation {
    Vert,
    Horiz,
}

impl Clone for Bond {
    fn clone(&self) -> Self {
        Self { order: self.order.clone(), orient: self.orient.clone() }
    }
}