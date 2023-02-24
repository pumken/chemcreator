use crate::molecule::Atom::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

#[derive(Clone)]
pub enum Symbol {
    Atom(Atom),
    Bond(Bond),
    None,
}

#[derive(Clone, Copy)]
pub enum Atom {
    C,
    H,
    O,
}

impl Atom {
    pub fn symbol(&self) -> &str {
        match *self {
            C => "[C]",
            H => "[H]",
            O => "[O]"
        }
    }
}

#[derive(Copy)]
pub struct Bond {
    pub order: BondOrder,
    pub orient: BondOrientation
}

#[derive(Clone, Copy)]
pub enum BondOrder {
    Single, Double, Triple
}

#[derive(Clone, Copy)]
pub enum BondOrientation {
    Vert, Horiz
}

impl Bond {
    pub(crate) fn new(order: BondOrder, orient: BondOrientation) -> Bond {
        Bond { order, orient }
    }

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

    pub(crate) fn is_horizontal(&self) -> bool {
        match self.orient {
            Horiz => true,
            Vert => false
        }
    }

    pub(crate) fn is_vertical(&self) -> bool {
        match self.orient {
            Horiz => false,
            Vert => true
        }
    }
}

impl Clone for Bond {
    fn clone(&self) -> Self {
        Self { order: self.order.clone(), orient: self.orient.clone() }
    }
}