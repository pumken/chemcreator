use std::ops::{Index, IndexMut};
use ruscii::spatial::Vec2;
use crate::molecule::Atom::{C, H};
use crate::molecule::{GroupType, Symbol};

pub struct GridState {
    pub cells: Vec<Cell>,
    // privatize when i figure out how to implement an iterator
    pub size: Vec2,
    pub cursor: Vec2,
}

impl GridState {
    pub fn new(width: usize, height: usize) -> GridState {
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

    pub fn insert(&mut self, symbol: Symbol) {
        let x = self.cursor.x as usize;
        let y = self.cursor.y as usize;
        self[x][y].sym = symbol;
    }

    pub fn current_cell(&self) -> &Cell {
        &self[self.cursor.x as usize][self.cursor.y as usize]
    }

    pub fn find<F>(&self, predicate: F) -> Option<&Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        self.cells.iter()
            .find(|&element| predicate(element))
    }

    pub fn find_all<F>(&self, predicate: F) -> Vec<Cell>
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

    pub fn count<F>(&self, predicate: F) -> i32
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

    pub fn atom_adjacent(&self) -> bool {
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

    pub fn chain_endpoints(&self) -> Vec<Cell> {
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

    pub fn simply_connected(&self) -> bool {
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

    pub fn simple_counter(&self) -> String {
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
                    return "".to_string();
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

        if begin < 0 || end >= self.cells.len() {
            panic!("index out of bounds: the len is {} but the index is {index}", &self.size.x)
        }
        &self.cells[begin..end]
    }
}

impl IndexMut<usize> for GridState {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let begin = *&self.size.y as usize * index;
        let end = *&self.size.y as usize * (index + 1);

        if begin < 0 || end >= self.cells.len() {
            panic!("index out of bounds: the len is {} but the index is {index}", &self.size.x)
        }
        &mut self.cells[begin..end]
    }
}

// Fix later
// impl Iterator for GridState {
//     type Item = Vec<&mut Cell>;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         self.cells
//             .group_by(|a, b| a.pos.x == b.pos.x)
//             .collect()
//             .next()
//     }
// }

pub struct Molecule {
    pub groups: Vec<Group>,
}

pub struct Group {
    pub cells: Vec<Cell>,
    pub class: GroupType,
}

#[derive(Clone)]
pub struct Cell {
    pub sym: Symbol,
    pub pos: Vec2,
}

pub struct Neighbors {
    pub pos: Vec2,
    pub neighbors: Vec<Cell>,
}

impl Neighbors {
    pub fn get(graph: &GridState, pos: Vec2) -> Neighbors {
        let adjacents = vec![
            (1, 0, true),
            (0, 1, false),
            (-1, 0, true),
            (0, -1, false),
        ];
        let mut neighbors: Vec<Cell> = vec!();

        for (x, y, isHorizontal) in adjacents {
            match &graph[(pos.x + x) as usize][(pos.y + y) as usize].sym {
                Symbol::Bond(it) => {
                    if isHorizontal == it.is_horizontal() {
                        neighbors.push(graph[(pos.x + x * 2) as usize][(pos.y + y * 2) as usize].clone())
                    }
                }
                _ => ()
            }
        }
        Neighbors { pos, neighbors }
    }

    pub fn count<F>(&self, predicate: F) -> i32
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