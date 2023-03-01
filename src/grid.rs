use ruscii::spatial::Vec2;
use crate::molecule::Atom::{C, H};
use crate::molecule::{GroupType, Symbol};

#[derive(Clone)]
pub(crate) struct GridState {
    pub(crate) cells: Vec<Vec<Cell>>,
    pub(crate) size: Vec2,
    pub(crate) cursor: Vec2,
}

impl GridState {
    pub(crate) fn new(width: usize, height: usize) -> GridState {
        let mut cells = vec![];

        for x in (0 as usize)..width {
            let mut row = vec![];
            for y in (0 as usize)..height {
                let new_cell = Cell { sym: Symbol::None, pos: Vec2::xy(x, y) };
                row.push(new_cell);
            }
            cells.push(row);
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
        self.cells[x][y].sym = symbol;
    }

    pub(crate) fn current_cell(&self) -> &Cell {
        &self.cells[self.cursor.x as usize][self.cursor.y as usize]
    }

    pub(crate) fn count(&self, predicate: fn(&Cell) -> bool) -> i32 {
        let mut count = 0;
        for cell in self.cells.iter().flatten() {
            if predicate(&cell) {
                count += 1;
            }
        }
        count
    }

    pub(crate) fn atom_adjacent(&self) -> bool {
        let left = &self.cells[self.cursor.x as usize - 1][self.cursor.y as usize];
        let right = &self.cells[self.cursor.x as usize + 1][self.cursor.y as usize];
        match left.sym {
            Symbol::Atom(_) => true,
            _ => match right.sym {
                Symbol::Atom(_) => true,
                _ => false
            }
        }
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