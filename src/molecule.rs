//! # Molecule
//!
//! The `molecule` module provides functionality for representing molecular components.

use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::molecule::Element::{C, H, N, O};
use ruscii::spatial::{Direction, Vec2};
use ruscii::terminal::Color;
use ruscii::terminal::Color::{Blue, Green, LightGrey, Magenta, Red, White, Xterm, Yellow};
use std::cmp::Ordering;
use std::fmt;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

/// Represents a type of functional group on a molecule.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Group {
    Alkane,
    Alkene,
    Alkyne,
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

impl Group {
    /// Returns the priority level of each functional group for identifying the main functional
    /// group of a molecule.
    ///
    /// A higher priority is indicated by a higher number. The lowest
    /// priority group is assigned zero. A value of [`None`] indicates that the group is _never
    /// the main group_ (i.e., always a prefix).
    pub const fn priority(self) -> Option<i32> {
        let priority = match self {
            Group::Alkane | Group::Alkene | Group::Alkyne => 0,
            Group::Methyl
            | Group::Ethyl
            | Group::Propyl
            | Group::Isopropyl
            | Group::Butyl
            | Group::Pentyl
            | Group::Hexyl
            | Group::Heptyl
            | Group::Octyl
            | Group::Nonyl
            | Group::Decyl => return None,
            Group::Bromo | Group::Chloro | Group::Fluoro | Group::Iodo => return None,
            Group::Hydroxyl => 3,
            Group::Carbonyl => 4,
            Group::Carboxyl => 6,
            Group::Ester => 5,
            Group::Ether => return None,
        };
        Some(priority)
    }

    pub const fn is_chain_group(self) -> bool {
        matches!(self, Group::Alkane | Group::Alkene | Group::Alkyne)
    }
}

impl Display for Group {
    /// Displays the ionic prefix for the current [`Group`].
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        let str = match *self {
            Group::Alkane => "an",
            Group::Alkene => "en",
            Group::Alkyne => "yn",
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
            Group::Carbonyl => "oxo",
            Group::Carboxyl => "carboxyl",
            Group::Ester => "ester",
            Group::Ether => "ether",
        };
        write!(f, "{str}")
    }
}

impl PartialOrd for Group {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        todo!()
    }
}

impl FromStr for Group {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let out = match s {
            "alkane" => Group::Alkane,
            "alkene" => Group::Alkene,
            "alkyne" => Group::Alkyne,
            "bromo" => Group::Bromo,
            "chloro" => Group::Chloro,
            "fluoro" => Group::Fluoro,
            "iodo" => Group::Iodo,
            "hydroxyl" => Group::Hydroxyl,
            "oxo" => Group::Carbonyl,
            "carboxyl" => Group::Carboxyl,
            "ester" => Group::Ester,
            "ether" => Group::Ether,
            _ => return Err(()),
        };
        Ok(out)
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
                Element::Br => Xterm(1),
                C => LightGrey,
                Element::Cl => Green,
                Element::F => Yellow,
                Element::I => Magenta,
                N => Blue,
                O => Red,
                _ => White,
            },
            _ => White,
        }
    }

    pub fn pos(&self) -> Vec2 {
        match self {
            Cell::Atom(it) => it.pos,
            Cell::Bond(it) => it.pos,
            Cell::None(it) => *it,
        }
    }

    pub fn is_atom(&self) -> bool {
        matches!(self, Cell::Atom(_))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub pos: Vec2,
}

impl Atom {
    pub const fn symbol(&self) -> &str {
        match self.element {
            Element::Br => "[Br]",
            C => "[C]",
            Element::Cl => "[Cl]",
            Element::F => "[F]",
            H => "[H]",
            Element::I => "[I]",
            N => "[N]",
            O => "[O]",
        }
    }
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            element: C,
            pos: Vec2::zero(),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Element {
    Br,
    C,
    Cl,
    F,
    H,
    I,
    N,
    O,
}

impl Element {
    /// Returns the number of bonds the current [`Element`] should have.
    pub const fn bond_number(&self) -> u8 {
        match *self {
            C => 4,
            N => 3,
            O => 2,
            H | Element::Br | Element::Cl | Element::F | Element::I => 1,
        }
    }

    /// Returns the single-character identifier of the [`Element`].
    pub const fn id(&self) -> &str {
        match *self {
            Element::Br => "B",
            C => "C",
            Element::Cl => "L",
            Element::F => "F",
            H => "H",
            Element::I => "I",
            N => "N",
            O => "O",
        }
    }
}

impl Display for Element {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Element::Br => "Bromine",
                C => "Carbon",
                Element::Cl => "Chlorine",
                Element::F => "Fluorine",
                H => "Hydrogen",
                Element::I => "Iodine",
                N => "Nitrogen",
                O => "Oxygen",
            }
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub pos: Vec2,
    pub order: BondOrder,
    pub orient: BondOrientation,
}

impl Bond {
    pub fn symbol(&self) -> &str {
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
pub enum BondOrientation {
    Vert,
    Horiz,
}

impl From<Direction> for BondOrientation {
    /// Returns the related [`BondOrientation`] to the given `direction`.
    ///
    /// ## Panics
    ///
    /// If [`Direction::None`] is passed to this function.
    fn from(value: Direction) -> BondOrientation {
        match value {
            Direction::Up | Direction::Down => Vert,
            Direction::Left | Direction::Right => Horiz,
            Direction::None => panic!("Attempted to pass Direction::None to from_direction."),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Substituent {
    Branch(Branch),
    Group(Group),
}

impl Display for Substituent {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Substituent::Branch(_) => "".to_string(),
                Substituent::Group(it) => it.to_string(),
            }
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Branch {
    pub chain: Vec<Atom>,
    pub groups: Vec<Vec<Substituent>>,
}

impl Branch {
    /// Creates a new [`Branch`] with the given `chain` and with a `groups` field with the same
    /// capacity as the `chain`.
    pub fn new(chain: Vec<Atom>) -> Branch {
        let len = chain.len();
        Branch {
            chain,
            groups: Vec::with_capacity(len),
        }
    }
}

impl Display for Branch {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        let out = self
            .groups
            .iter()
            .enumerate()
            .map(|(index, it)| {
                format!(
                    "{index}: {}",
                    it.iter()
                        .map(|it| it.to_string())
                        .collect::<Vec<String>>()
                        .join(", ")
                )
            })
            .collect::<Vec<String>>()
            .join("; ");
        write!(f, "{out}")
    }
}

impl FromStr for Branch {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let groups = s
            .split("; ")
            .map(|it| {
                it.trim_start_matches(|c: char| c.is_ascii_digit() || c == ':' || c == ' ')
                    .split(", ")
                    .map(|str| Substituent::Group(Group::from_str(str).unwrap()))
                    .collect::<Vec<Substituent>>()
            })
            .collect::<Vec<Vec<Substituent>>>();

        let len = 0usize;
        let out = Branch {
            chain: vec![Atom::default(); len],
            groups,
        };
        Ok(out)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GroupNode {
    pub bond: BondOrder,
    pub atom: Element,
    pub next: Vec<GroupNode>,
}

impl Display for GroupNode {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        let primary = format!("{}{}", self.bond.id(), self.atom.id());

        let mut tree = self.next.clone();
        tree.sort_by_key(|node| node.to_string());
        let out = tree.iter().fold(primary, |a, b| format!("{a}({b})"));
        write!(f, "{}", out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_with;
    use crate::groups::group_node_tree;
    use crate::molecule::Group::{Bromo, Carbonyl, Hydroxyl};
    use crate::spatial::GridState;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn group_to_string() {
        let groups = vec![Group::Carboxyl, Carbonyl, Hydroxyl];

        assert_eq!(groups[0].to_string(), "carboxyl");
        assert_eq!(groups[1].to_string(), "oxo");
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
            groups: vec![
                vec![Substituent::Group(Hydroxyl), Substituent::Group(Carbonyl)],
                vec![Substituent::Group(Bromo)],
            ],
        };

        assert_eq!(branch.to_string(), "0: hydroxyl, oxo; 1: bromo")
    }

    #[test]
    fn branch_from_string_parses_str() {
        let branch = Branch::from_str("0: hydroxyl").unwrap();
        let expected = vec![vec![Substituent::Group(Hydroxyl)]];

        assert_eq!(branch.groups, expected);
    }

    #[test]
    fn branch_from_string_multiple_groups() {
        let branch = Branch::from_str("0: hydroxyl; 1: hydroxyl").unwrap();
        let expected = vec![
            vec![Substituent::Group(Hydroxyl)],
            vec![Substituent::Group(Hydroxyl)],
        ];

        assert_eq!(branch.groups, expected);
    }

    #[test]
    fn branch_from_string_multiple_compound_groups() {
        let branch = Branch::from_str("0: hydroxyl, oxo; 1: hydroxyl, oxo").unwrap();
        let expected = vec![
            vec![Substituent::Group(Hydroxyl), Substituent::Group(Carbonyl)],
            vec![Substituent::Group(Hydroxyl), Substituent::Group(Carbonyl)],
        ];

        assert_eq!(branch.groups, expected);
    }

    #[test]
    fn branch_new_contains_sized_vec() {
        let chain = vec![
            Atom {
                element: C,
                pos: Vec2::zero(),
            },
            Atom {
                element: C,
                pos: Vec2::zero(),
            },
        ];
        let branch = Branch::new(chain);

        assert_eq!(branch.groups.capacity(), 2usize);
    }

    #[test]
    fn group_node_to_string() {
        let graph = graph_with!(1, 3,
            [0, 0; A(C)],
            [0, 1; A(O)],
            [0, 2; A(H)]
        );
        let group_node = group_node_tree(&graph, Vec2::xy(0, 0), Direction::Up).unwrap();

        assert_eq!(group_node.to_string(), "1O(1H)");
    }
}
