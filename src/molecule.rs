//! # Molecule
//!
//! The `molecule` module provides functionality for representing molecular components.

use crate::molecule::BondOrder::{Double, Single, Triple};
use crate::molecule::BondOrientation::{Horiz, Vert};
use crate::molecule::Element::{C, H, N, O};
use crate::naming::process_name;
use crate::spatial::EnumAll;
use ruscii::spatial::{Direction, Vec2};
use ruscii::terminal::Color;
use ruscii::terminal::Color::{Blue, Green, LightGrey, Magenta, Red, White, Xterm, Yellow};
use std::cmp::Ordering;
use std::fmt;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

/// Represents a type of functional group on a molecule.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Group {
    /* Not a group, but useful */
    Hydrogen,
    /* Multiple bond groups */
    Alkane,
    Alkene,
    Alkyne,
    /* Halide groups */
    Bromo,
    Chloro,
    Fluoro,
    Iodo,
    /* General groups */
    Hydroxyl,
    Carbonyl,
    /* Compound groups */
    Aldehyde,
    AcidHalide(Halogen),
    Carboxyl,
    /* Chain groups */
    Ester,
    Ether,
    /* Nitrogen groups */
    Nitrile,
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
            Group::Carboxyl => 12,
            Group::Ester => 11,
            Group::AcidHalide(halogen) => match halogen {
                Halogen::Fluorine => 10,
                Halogen::Chlorine => 9,
                Halogen::Bromine => 8,
                Halogen::Iodine => 7,
            },
            Group::Nitrile => 6,
            Group::Aldehyde => 5,
            Group::Carbonyl => 4,
            Group::Hydroxyl => 3,
            Group::Alkane | Group::Alkene | Group::Alkyne => 0,
            Group::Hydrogen
            | Group::Bromo
            | Group::Chloro
            | Group::Fluoro
            | Group::Iodo
            | Group::Ether => return None,
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
            Group::Hydrogen => "",
            Group::Alkane => "an",
            Group::Alkene => "en",
            Group::Alkyne => "yn",
            Group::Bromo => "bromo",
            Group::Chloro => "chloro",
            Group::Fluoro => "fluoro",
            Group::Iodo => "iodo",
            Group::Hydroxyl => "hydroxy",
            Group::Carbonyl => "oxo",
            Group::Aldehyde => "formyl",
            Group::Carboxyl => "carboxy",
            Group::AcidHalide(it) => return write!(f, "{}carbonyl", it.associated_group()),
            Group::Ester => "ester",
            Group::Ether => "ether",
            Group::Nitrile => "cyano",
        };
        write!(f, "{str}")
    }
}

impl PartialOrd<Self> for Group {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let out = match (self.priority(), other.priority()) {
            (Some(a), Some(b)) => a.cmp(&b),
            (Some(_), None) => Ordering::Greater,
            (None, Some(_)) => Ordering::Less,
            (None, None) => Ordering::Equal,
        };
        Some(out)
    }
}

impl FromStr for Group {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let out = match s {
            "hydrogen" => Group::Hydrogen,
            "alkane" => Group::Alkane,
            "alkene" => Group::Alkene,
            "alkyne" => Group::Alkyne,
            "bromo" => Group::Bromo,
            "chloro" => Group::Chloro,
            "fluoro" => Group::Fluoro,
            "iodo" => Group::Iodo,
            "hydroxy" => Group::Hydroxyl,
            "oxo" => Group::Carbonyl,
            "carboxy" => Group::Carboxyl,
            "ester" => Group::Ester,
            "ether" => Group::Ether,
            "cyano" => Group::Nitrile,
            _ => return Err(()),
        };
        Ok(out)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Halogen {
    Fluorine,
    Chlorine,
    Bromine,
    Iodine,
}

impl Halogen {
    pub fn associated_group(&self) -> Group {
        match self {
            Halogen::Fluorine => Group::Fluoro,
            Halogen::Chlorine => Group::Chloro,
            Halogen::Bromine => Group::Bromo,
            Halogen::Iodine => Group::Iodo,
        }
    }
}

impl Display for Halogen {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let str = match self {
            Halogen::Fluorine => "fluoride",
            Halogen::Chlorine => "chloride",
            Halogen::Bromine => "bromide",
            Halogen::Iodine => "iodide",
        };
        write!(f, "{str}")
    }
}

impl EnumAll for Halogen {
    fn all() -> Vec<Self>
    where
        Self: Sized,
    {
        vec![
            Halogen::Fluorine,
            Halogen::Chlorine,
            Halogen::Bromine,
            Halogen::Iodine,
        ]
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

    /// Returns the component type of the atom.
    pub fn comp(&self) -> ComponentType {
        match self {
            Cell::Atom(it) => ComponentType::Element(it.element),
            Cell::Bond(it) => ComponentType::Order(it.order),
            Cell::None(_) => ComponentType::None,
        }
    }

    pub fn is_atom(&self) -> bool {
        matches!(self, Cell::Atom(_))
    }

    pub fn is_empty(&self) -> bool {
        !matches!(self, Cell::None(_))
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum ComponentType {
    Element(Element),
    Order(BondOrder),
    None,
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

impl Substituent {
    pub fn priority(&self) -> Option<i32> {
        match self {
            Substituent::Group(it) => it.priority(),
            Substituent::Branch(_) => None,
        }
    }

    pub fn is_chain_group(&self) -> bool {
        match self {
            Substituent::Group(it) => it.is_chain_group(),
            Substituent::Branch(_) => false,
        }
    }
}

impl Display for Substituent {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Substituent::Branch(it) => it.to_string(),
                Substituent::Group(it) => it.to_string(),
            }
        )
    }
}

#[derive(Clone, Debug)]
pub struct Branch {
    pub chain: Vec<Atom>,
    pub groups: Vec<Vec<Substituent>>,
    pub parent_alpha: Option<Atom>,
}

impl Branch {
    /// Creates a new [`Branch`] with the given `chain` and with a `groups` field with the same
    /// capacity as the `chain`.
    pub fn new(chain: Vec<Atom>) -> Branch {
        let len = chain.len();
        Branch {
            chain,
            groups: Vec::with_capacity(len),
            parent_alpha: None,
        }
    }

    pub fn reversed(&self) -> Branch {
        let (mut chain, mut groups, parent_alpha) = {
            let clone = self.clone();
            (clone.chain, clone.groups, clone.parent_alpha)
        };
        chain.reverse();

        groups = Branch::shift_chain_groups(groups);
        groups.reverse();

        Branch {
            chain,
            groups,
            parent_alpha,
        }
    }

    fn shift_chain_groups(mut groups: Vec<Vec<Substituent>>) -> Vec<Vec<Substituent>> {
        for outer in (1..groups.len()).rev() {
            let link = &mut groups[outer - 1];
            let mut cache = vec![];

            for index in (0..link.len()).rev() {
                if link[index].is_chain_group() {
                    let element = link.remove(index);
                    cache.push(element);
                }
            }

            groups[outer].append(&mut cache);
        }
        groups
    }
}

impl PartialEq for Branch {
    fn eq(&self, other: &Self) -> bool {
        self.chain == other.chain && self.groups == other.groups
    }
}

impl Display for Branch {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        let mut name = match process_name(self.clone()) {
            Ok(it) => it,
            Err(_) => panic!("unable to process {:?}", self),
        };
        name.truncate(name.len() - 3);
        write!(f, "{name}yl")
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
                    .map(|str| {
                        Substituent::Group(match Group::from_str(str) {
                            Ok(it) => it,
                            Err(_) => panic!("\"{str}\" is not a valid group"),
                        })
                    })
                    .collect::<Vec<Substituent>>()
            })
            .collect::<Vec<Vec<Substituent>>>();

        let len = 0usize;
        let out = Branch {
            chain: vec![Atom::default(); len],
            groups,
            parent_alpha: None,
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
    use crate::molecule::Group::{Alkene, Bromo, Carbonyl, Hydroxyl};
    use crate::spatial::GridState;
    use crate::test_utils::GW::{A, B};

    #[test]
    fn group_to_string() {
        let groups = vec![Group::Carboxyl, Carbonyl, Hydroxyl];

        assert_eq!(groups[0].to_string(), "carboxy");
        assert_eq!(groups[1].to_string(), "oxo");
        assert_eq!(groups[2].to_string(), "hydroxy");
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

        assert_eq!(a.to_string(), "hydroxy".to_string());
        assert_eq!(b.to_string(), "bromo".to_string())
    }

    #[test]
    fn branch_reversed_correctly_reverses_chain() {
        let branch = Branch::from_str("0: hydroxy, bromo, alkene; 1: chloro").unwrap();
        let reversed = branch.reversed();
        let expected = Branch::from_str("0: chloro, alkene; 1: hydroxy, bromo").unwrap();

        assert_eq!(reversed, expected)
    }

    #[test]
    fn branch_shift_chain_groups() {
        let mut groups = vec![vec![Substituent::Group(Alkene)], vec![]];
        groups = Branch::shift_chain_groups(groups);
        let expected = vec![vec![], vec![Substituent::Group(Alkene)]];

        assert_eq!(groups, expected);
    }

    #[test]
    fn branch_to_string_groups_only() {
        let branch = Branch {
            chain: vec![/* Not used in to_string() */],
            groups: vec![
                vec![Substituent::Group(Hydroxyl), Substituent::Group(Carbonyl)],
                vec![Substituent::Group(Bromo)],
            ],
            parent_alpha: None,
        };
        let expected = Branch::from_str("0: hydroxy, oxo; 1: bromo").unwrap();

        assert_eq!(branch, expected);
    }

    #[test]
    fn branch_from_string_parses_str() {
        let branch = Branch::from_str("0: hydroxy").unwrap();
        let expected = vec![vec![Substituent::Group(Hydroxyl)]];

        assert_eq!(branch.groups, expected);
    }

    #[test]
    fn branch_from_string_multiple_groups() {
        let branch = Branch::from_str("0: hydroxy; 1: hydroxy").unwrap();
        let expected = vec![
            vec![Substituent::Group(Hydroxyl)],
            vec![Substituent::Group(Hydroxyl)],
        ];

        assert_eq!(branch.groups, expected);
    }

    #[test]
    fn branch_from_string_multiple_compound_groups() {
        let branch = Branch::from_str("0: hydroxy, oxo; 1: hydroxy, oxo").unwrap();
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
