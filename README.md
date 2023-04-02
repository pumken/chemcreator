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
- [x] Aldehydes
- [x] Amides
- [x] Amines
- [x] Carboxylic acids
- [ ] Esters
- [ ] Ethers
- [x] Halogenoalkanes
- [x] Ketones
- [x] Nitriles

## Usage

### Insert Mode

In insert mode, you can edit the grid to create a molecule. If the structure is valid and ChemCreator recognizes it,
it'll give you a name for it.

![](res/screenshot-insert.png)

## Display Mode

In display mode, a clean picture is shown of your molecule as well as some basic statistics including atom count, atomic
weight, index of hydrogen deficiency (IHD), and the length of the name in characters.

![](res/screenshot-normal.png)

> ⓘ  These screenshots were taken from CLion's inbuilt terminal, but you could just as easily use a terminal of your 
> choice.
