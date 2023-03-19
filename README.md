<div align="center">
<h1>ChemCreator</h1>

<img src="https://img.shields.io/github/v/release/pumken/chemcreator?include_prereleases"></img>
[![CI](https://github.com/pumken/chemcreator/actions/workflows/CI.yml/badge.svg)](https://github.com/pumken/chemcreator/actions/workflows/CI.yml)
<img src="https://img.shields.io/github/last-commit/pumken/chemcreator"></img>
<img src="https://img.shields.io/github/languages/code-size/pumken/chemcreator"></img>

A text-based tool for identifying organic molecules.

Written in Rust with the [`ruscii`](https://github.com/lemunozm/ruscii) library.
</div>

## Preface

ChemCreator is a program you can open in your terminal to generate the name of
an organic molecule. You recreate the structure in a grid, manually inputting
atoms and bonds by using your keyboard. When a valid molecule is found,
ChemCreator will try its best to give you the molecule's IUPAC systematic name.

ChemCreator is a personal project that I started to explore the Rust programming
language and also try out a moderately difficult problem in traversing a grid to
find certain structures, determine patterns, etc.

It's not perfect; it may get a few names wrong here and there at times, so if
you try it out and find something named incorrectly, please open an issue!

## Features

- [x] Easy-to-use text-based UI
- [ ] Variable-size grid

### Supported Molecule Types

- [x] Alkanes
- [x] Alkenes
- [x] Alkynes
- [x] Alcohols
- [ ] Aldehydes
- [x] Carboxylic acids
- [ ] Ethers
- [x] Ketones
- [ ] Halogenoalkanes
- [ ] Esters
