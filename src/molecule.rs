//! # Molecule
//!
//! The `molecule` module provides functionality for representing molecular components.

use std::fmt::{Display, Formatter};
use ruscii::spatial::{Direction, Vec2};
use ruscii::terminal::Color;
use ruscii::terminal::Color::{LightGrey, Red, White};
use crate::spatial::Cellular;
use crate::molecule::Element::{C, H, O};
use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};

/// Represents a type of functional group on a molecule.
#[derive(Clone, PartialEq)]
pub enum Group {
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

impl ToString for Group {
    fn to_string(&self) -> String {
        match *self {
            Group::Methyl => "methyl",
            Group::Ethyl => "ethyl",
            Group::Propyl => "propyl",
            Group::Isopropyl => "isopropyl",
            Group::Butyl => "butyl",
            Group::Pentyl => "pentyl",
            Group::Hexyl => "hexyl",
            Group::Heptyl => "heptyl",
            Group::Octyl => "octyl",
            Group::Nonyl => "nonyl",
            Group::Decyl => "decyl",
            Group::Bromo => "bromo",
            Group::Chloro => "chloro",
            Group::Fluoro => "fluoro",
            Group::Iodo => "iodo",
            Group::Hydroxyl => "hydroxyl",
            Group::Carbonyl => "carbonyl",
            Group::Carboxyl => "carboxyl",
            Group::Ester => "ester",
            Group::Ether => "ether",
        }.to_string()
    }
}

/// Represents a molecular component.
#[derive(Clone, Debug, PartialEq)]
pub enum Cell {
    Atom(Atom),
    Bond(Bond),
    None(Vec2),
}

impl Cell {
    pub fn color(&self) -> Color {
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

    pub(crate) fn is_atom(&self) -> bool {
        matches!(self, Cell::Atom(_))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub(crate) element: Element,
    pub(crate) pos: Vec2,
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
pub enum Element {
    C,
    H,
    O,
}

impl Element {
    /// Returns the number of bonds the current [`Element`] should have.
    pub(crate) const fn bond_number(&self) -> u8 {
        match *self {
            C => 4,
            H => 1,
            O => 2
        }
    }

    pub(crate) const fn id(&self) -> &str {
        match *self {
            C => "C",
            H => "H",
            O => "O",
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

#[derive(Copy, Debug, PartialEq)]
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

impl BondOrder {
    pub const fn order(&self) -> u8 {
        match self {
            Single => 1,
            Double => 2,
            Triple => 3,
        }
    }

    pub fn id(&self) -> String {
        format!("{}", self.order())
    }
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
        Bond { pos: self.pos, order: self.order, orient: self.orient }
    }
}

pub enum Substituent {
    Branch(Branch),
    Group(Group),
}

impl ToString for Substituent {
    fn to_string(&self) -> String {
        match self {
            Substituent::Branch(_) => "".to_string(),
            Substituent::Group(it) => it.to_string(),
        }
    }
}

pub struct Branch {
    pub chain: Vec<Atom>,
    pub groups: Vec<Vec<Substituent>>,
}

impl ToString for Branch {
    fn to_string(&self) -> String {
        let mut index = 0;
        self.groups.iter()
            .fold("".to_string(), |a, b| {
                let str = format!("{a}{index}: {} ", b.iter().fold("".to_string(), |c, d| {
                    format!("{c}{} ", d.to_string())
                }));
                index += 1;
                str
            })
    }
}

impl Branch {
    pub fn new(chain: Vec<Atom>) -> Branch {
        let len = chain.len();
        Branch {
            chain,
            groups: Vec::with_capacity(len),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GroupNode {
    pub bond: BondOrder,
    pub atom: Element,
    pub next: Vec<GroupNode>,
}

impl ToString for GroupNode {
    fn to_string(&self) -> String {
        let start = format!("{}{}", self.bond.id(), self.atom.id());

        let mut mid = self.next.clone();
        mid.sort_by_key(|a| a.to_string());
        mid.iter()
            .fold(start, |a, b| format!("{a}({})", b.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use crate::groups::group_node_tree;
    use crate::graph_with;
    use crate::molecule::Group::{Bromo, Carbonyl, Hydroxyl};
    use crate::spatial::GridState;
    use crate::test_utils::GW::{A, B};
    use super::*;

    #[test]
    fn group_to_string() {
        let groups = vec![Group::Carboxyl, Group::Carbonyl, Group::Hydroxyl];

        assert_eq!(groups[0].to_string(), "carboxyl");
        assert_eq!(groups[1].to_string(), "carbonyl");
        assert_eq!(groups[2].to_string(), "hydroxyl");
    }

    #[test]
    fn cell_color() {
        let graph = graph_with!(2, 2,
            [0, 0; A(C)],
            [1, 0; A(H)],
            [0, 1; A(O)],
            [1, 1; B(Single)]
        );

        assert_eq!(graph.get(Vec2::xy(0, 0)).unwrap().color(), LightGrey);
        assert_eq!(graph.get(Vec2::xy(1, 0)).unwrap().color(), White);
        assert_eq!(graph.get(Vec2::xy(0, 1)).unwrap().color(), Red);
        assert_eq!(graph.get(Vec2::xy(1, 1)).unwrap().color(), White);
    }

    #[test]
    fn cell_pos() {
        let graph = graph_with!(2, 2,
            [1, 0; A(H)],
            [0, 1; A(O)],
            [1, 1; B(Single)]
        );

        assert_eq!(graph.cells[0][1].pos(), Vec2::xy(0, 1));
        assert_eq!(graph.cells[1][1].pos(), Vec2::xy(1, 1));
    }

    #[test]
    fn substituent_to_string_for_group() {
        let a = Substituent::Group(Hydroxyl);
        let b = Substituent::Group(Bromo);

        assert_eq!(a.to_string(), "hydroxyl".to_string());
        assert_eq!(b.to_string(), "bromo".to_string())
    }

    #[test]
    fn branch_to_string_groups_only() {
        let branch = Branch {
            chain: vec![/* Not used in to_string() */],
            groups: vec![vec![Substituent::Group(Hydroxyl), Substituent::Group(Carbonyl)], vec![Substituent::Group(Bromo)]],
        };

        assert_eq!(branch.to_string(), "0: hydroxyl carbonyl  1: bromo  ")
    }

    #[test]
    fn group_node_to_string() {
        let graph = graph_with!(1, 3,
            [0, 0; A(C)],
            [0, 1; A(O)],
            [0, 2; A(H)]
        );
        let group_node = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up);

        assert_eq!(group_node.to_string(), "1O(1H)");
    }
}
