use crate::ReadError;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MapError {
    InvalidRead(ReadError),
    AlgorithmNotReady,
}

impl std::fmt::Display for MapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRead(error) => write!(f, "invalid read: {error}"),
            Self::AlgorithmNotReady => f.write_str("RS-LRA mapping pipeline is not ported yet"),
        }
    }
}

impl std::error::Error for MapError {}
