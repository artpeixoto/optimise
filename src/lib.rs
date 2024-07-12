mod base;
pub use base::*;

pub mod algoritmos;
pub mod ferramentas;
use std::iter;
use nalgebra::{DMatrix, SMatrix};
use replace_with::replace_with_or_abort_and_return;