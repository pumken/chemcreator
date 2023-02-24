use crate::molecule::Atom::{C, H, O};
use crate::molecule::Bond::{DoubleH, DoubleV, SingleH, SingleV, TripleH, TripleV};

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

#[derive(Clone, Copy)]
pub enum Bond {
    SingleV,
    SingleH,
    DoubleV,
    DoubleH,
    TripleV,
    TripleH,
}

impl Bond {
    pub(crate) fn symbol(&self) -> &str {
        match *self {
            SingleH => "———",
            SingleV => " | ",
            DoubleH => "===",
            DoubleV => " ‖ ",
            TripleH => "≡≡≡",
            TripleV => " T ",
        }
    }
}