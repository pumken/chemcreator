use std::iter::Cloned;
use itertools::Itertools;
use std::ops::{Index, IndexMut};
use std::slice::Iter;
use std::vec::IntoIter;
use ruscii::spatial::Vec2;
use crate::molecule::Atom::{C, H};
use crate::molecule::{GroupType, Symbol};

#[derive(Clone)]
pub(crate) struct GridState {
    cells: Vec<Cell>,
    pub(crate) size: Vec2,
    pub(crate) cursor: Vec2,
}

impl GridState {
    pub(crate) fn new(width: usize, height: usize) -> GridState {
        let mut cells = Vec::with_capacity(width * height);

        for x in (0 as usize)..width {
            for y in (0 as usize)..height {
                let new_cell = Cell { sym: Symbol::None, pos: Vec2::xy(x, y) };
                cells.push(new_cell);
            }
        }

        GridState {
            cells,
            size: Vec2::xy(width, height),
            cursor: Vec2::xy(0, 0),
        }
    }

    pub(crate) fn insert(&mut self, symbol: Symbol) {
        let x = self.cursor.x as usize;
        let y = self.cursor.y as usize;
        self[x][y].sym = symbol;
    }

    pub(crate) fn current_cell(&self) -> &Cell {
        &self[self.cursor.x as usize][self.cursor.y as usize]
    }

    pub(crate) fn find<F>(&self, predicate: F) -> Option<&Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        self.cells.iter()
            .find(|&element| predicate(element))
    }

    pub(crate) fn find_all<F>(&self, predicate: F) -> Vec<Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        let mut cells: Vec<Cell> = vec![];
        for cell in &self.cells[..] {
            if predicate(&cell) {
                cells.push(cell.clone());
            }
        }
        cells
    }

    pub(crate) fn count<F>(&self, predicate: F) -> i32
        where
            F: Fn(&Cell) -> bool,
    {
        let mut count = 0;
        for cell in &self.cells[..] {
            if predicate(&cell) {
                count += 1;
            }
        }
        count
    }

    fn grouped(&self) -> Vec<Vec<Cell>> {
        self.cells
            .iter()
            .group_by(|a| a.pos.x)
            .into_iter()
            .map(|(_, group)| group.map(|c| c.to_owned()).collect())
            .collect::<Vec<Vec<Cell>>>()
    }

    pub(crate) fn into_iter(self) -> IntoIter<Vec<Cell>> {
        self.grouped().into_iter()
    }

    /// Returns the length of a column of the grid.
    pub(crate) fn len(&self) -> usize {
        self.size.x as usize
    }

    pub(crate) fn atom_adjacent(&self) -> bool {
        let left = &self[self.cursor.x as usize - 1][self.cursor.y as usize];
        let right = &self[self.cursor.x as usize + 1][self.cursor.y as usize];
        match left.sym {
            Symbol::Atom(_) => true,
            _ => match right.sym {
                Symbol::Atom(_) => true,
                _ => false
            }
        }
    }

    pub(crate) fn chain_endpoints(&self) -> Vec<Cell> {
        let mut endpoints: Vec<Cell> = vec![];
        for cell in self.find_all(|element| match element.sym {
            Symbol::Atom(it) => match it {
                C => true,
                _ => false
            },
            _ => false
        }) {
            if Neighbors::get(&self, cell.pos).count(|element| match element.sym {
                Symbol::Atom(it) => match it {
                    C => true,
                    _ => false
                },
                _ => false
            }) <= 1 {
                endpoints.push(cell);
            }
        };
        endpoints
    }

    pub(crate) fn simply_connected(&self) -> bool {
        for cell in self.cells.iter() {
            let temp;
            match cell.sym {
                Symbol::Atom(C) => temp = cell,
                _ => continue
            }
            let neighbors = Neighbors::get(&self, temp.pos);
        };
        false
    }

    pub(crate) fn simple_counter(&self) -> String {
        let mut success = true;
        {
            let carbon = self.count(|element| match element.sym {
                Symbol::Atom(it) => return match it {
                    C => true,
                    _ => false
                },
                _ => false
            });
            let hydrogen = self.count(|element| match element.sym {
                Symbol::Atom(it) => return match it {
                    H => true,
                    _ => false
                },
                _ => false
            });
            let prefix = match carbon {
                1 => "meth",
                2 => "eth",
                3 => "prop",
                4 => "but",
                5 => "pent",
                6 => "hex",
                7 => "hept",
                8 => "oct",
                9 => "non",
                10 => "dec",
                _ => {
                    success = false;
                    ""
                }
            };
            let suffix = if hydrogen == 2 * carbon - 2 {
                "yne"
            } else if hydrogen == 2 * carbon {
                "ene"
            } else if hydrogen == 2 * carbon + 2 {
                "ane"
            } else {
                success = false;
                ""
            };
            if success {
                format!("{}{}", prefix, suffix)
            } else {
                "unidentified".to_string()
            }
        }
    }
}

impl Index<usize> for GridState {
    type Output = [Cell];

    fn index(&self, index: usize) -> &Self::Output {
        let begin = *&self.size.y as usize * index;
        let end = *&self.size.y as usize * (index + 1);

        if end > self.cells.len() {
            panic!("index out of bounds: the len is {} but the index is {index}", &self.size.x)
        }
        &self.cells[begin..end]
    }
}

impl IndexMut<usize> for GridState {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let begin = *&self.size.y as usize * index;
        let end = *&self.size.y as usize * (index + 1);

        if end > self.cells.len() {
            panic!("index out of bounds: the len is {} but the index is {index}", &self.size.x)
        }
        &mut self.cells[begin..end]
    }
}

pub(crate) struct Molecule {
    pub(crate) groups: Vec<Group>,
}

pub(crate) struct Group {
    pub(crate) cells: Vec<Cell>,
    pub(crate) class: GroupType,
}

#[derive(Clone)]
pub(crate) struct Cell {
    pub(crate) sym: Symbol,
    pub(crate) pos: Vec2,
}

pub(crate) struct Neighbors {
    pub(crate) pos: Vec2,
    pub(crate) neighbors: Vec<Cell>,
}

impl Neighbors {
    pub(crate) fn get(graph: &GridState, pos: Vec2) -> Neighbors {
        let adjacents = vec![
            (1, 0, true),
            (0, 1, false),
            (-1, 0, true),
            (0, -1, false),
        ];
        let mut neighbors: Vec<Cell> = vec!();

        for (x, y, is_horizontal) in adjacents {
            match &graph[(pos.x + x) as usize][(pos.y + y) as usize].sym {
                Symbol::Bond(it) => {
                    if is_horizontal == it.is_horizontal() {
                        neighbors.push(graph[(pos.x + x * 2) as usize][(pos.y + y * 2) as usize].clone())
                    }
                }
                _ => ()
            }
        }
        Neighbors { pos, neighbors }
    }

    pub(crate) fn count<F>(&self, predicate: F) -> i32
        where
            F: Fn(&Cell) -> bool,
    {
        let mut count = 0;
        for cell in &self.neighbors {
            if predicate(&cell) {
                count += 1;
            }
        }
        count
    }
}