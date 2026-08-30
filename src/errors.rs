use crate::{AnchorError, ChainCigarError, ReadError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MapError {
    InvalidRead(ReadError),
    Anchor(AnchorError),
    Cigar(ChainCigarError),
}

impl std::fmt::Display for MapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRead(error) => write!(f, "invalid read: {error}"),
            Self::Anchor(error) => write!(f, "anchor discovery failed: {error}"),
            Self::Cigar(error) => write!(f, "chain CIGAR assembly failed: {error}"),
        }
    }
}

impl std::error::Error for MapError {}
