use thiserror::Error;
use crate::molecule::Cell;

#[derive(Error, Debug)]
pub(crate) enum InvalidGraphError<'a> {
    #[error("Molecule is not continuous.")]
    Discontinuity,
    #[error("Cell at ({}, {}) is missing bonds.", .0.pos().x, .0.pos().y)]
    UnfilledValence(&'a Cell),
    #[error("Cell at ({}, {}) has too many bonds.", .0.pos().x, .0.pos().y)]
    OverfilledValence(&'a Cell),
    #[error("Bond at ({}, {}) is incomplete.", .0.pos().x, .0.pos().y)]
    IncompleteBond(&'a Cell),
    #[error("This combination of groups is not supported.")]
    UnsupportedGroups,
}