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

    pub(crate) fn find<F>(&self, predicate: F) -> Option<&Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        self.cells.iter()
            .flatten()
            .find(|&element| predicate(element))
    }

    pub(crate) fn find_all<F>(&self, predicate: F) -> Vec<Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        let mut cells: Vec<Cell> = vec![];
        for cell in self.cells.iter().flatten() {
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
        for cell in self.cells.iter().flatten() {
            if predicate(&cell) {
                count += 1;
            }
        }
        count
    }

    /// Returns the length of a column of the grid.
    pub(crate) fn len(&self) -> usize {
        self.size.x as usize
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
            match &graph.cells[(pos.x + x) as usize][(pos.y + y) as usize].sym {
                Symbol::Bond(it) => {
                    if is_horizontal == it.is_horizontal() {
                        neighbors.push(graph.cells[(pos.x + x * 2) as usize][(pos.y + y * 2) as usize].clone())
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