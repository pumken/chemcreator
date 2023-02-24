use ruscii::spatial::Vec2;
use crate::molecule::Atom::{C, H};
use crate::molecule::Symbol;

pub struct GridState {
    pub cells: Vec<Vec<Cell>>,
    pub size: Vec2,
    pub cursor: Vec2,
}

#[derive(Clone)]
pub struct Cell {
    pub sym: Symbol,
    pub pos: Vec2,
}

impl GridState {
    pub fn new(width: i32, height: i32) -> GridState {
        let mut temp: Vec<Vec<Cell>> = Vec::with_capacity(width as usize);

        for x in (0 as usize)..(width as usize) {
            let mut temp1 = Vec::with_capacity(height as usize);
            for y in (0 as usize)..(height as usize) {
                let new_cell = Cell { sym: Symbol::None, pos: Vec2::xy(x, y) };
                temp1.push(new_cell);
            }
            temp.push(temp1);
        }

        GridState {
            cells: temp,
            size: Vec2::xy(width, height),
            cursor: Vec2::xy(0, 0),
        }
    }

    pub fn insert(&mut self, symbol: Symbol) {
        let x = self.cursor.x as usize;
        let y = self.cursor.y as usize;
        self.cells[x][y].sym = symbol;
    }

    pub fn current_cell(&self) -> &Cell {
        &self.cells[self.cursor.x as usize][self.cursor.y as usize]
    }

    pub fn atom_adjacent(&self) -> bool {
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

    pub fn find<F>(&self, predicate: F) -> Option<&Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        self.cells.iter()
            .flatten()
            .find(|&element| predicate(element))
    }

    pub fn count<F>(&self, predicate: F) -> i32
        where
            F: Fn(&Cell) -> bool,
    {
        let mut count = 0;
        for row in &self.cells[..] {
            for cell in row {
                if predicate(&cell) {
                    count += 1;
                }
            }
        }
        count
    }

    pub fn find_all<F>(&self, predicate: F) -> Vec<Cell>
        where
            F: Fn(&Cell) -> bool,
    {
        let mut cells: Vec<Cell> = vec![];
        for row in &self.cells[..] {
            for cell in row {
                if predicate(&cell) {
                    cells.push(cell.clone());
                }
            }
        }
        cells
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
        for row in self.cells.iter() {
            for cell in row {
                let temp;
                match cell.sym {
                    Symbol::Atom(C) => temp = cell,
                    _ => continue
                }
                let neighbors = Neighbors::get(&self, temp.pos);
            }
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
            match &graph.cells[(pos.x + x) as usize][(pos.y + y) as usize].sym {
                Symbol::Bond(it) => {
                    if isHorizontal == it.is_horizontal() {
                        neighbors.push(graph.cells[(pos.x + x * 2) as usize][(pos.y + y * 2) as usize].clone())
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