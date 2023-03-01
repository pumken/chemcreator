use ruscii::terminal::Color;
use ruscii::terminal::Color::{LightGrey, Red, White};
use crate::grid::GridState;
use crate::molecule::Atom::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

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

#[derive(Clone)]
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

#[derive(Clone, Copy)]
pub(crate) enum Atom {
    C,
    H,
    O,
}

#[derive(Copy)]
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

#[derive(Clone, Copy)]
pub(crate) enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Clone, Copy)]
pub(crate) enum BondOrientation {
    Vert,
    Horiz,
}

impl Atom {
    pub(crate) fn symbol(&self) -> &str {
        match *self {
            C => "[C]",
            H => "[H]",
            O => "[O]"
        }
    }
}

impl Clone for Bond {
    fn clone(&self) -> Self {
        Self { order: self.order.clone(), orient: self.orient.clone() }
    }
}