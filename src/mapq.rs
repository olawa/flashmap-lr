//! Dataset-specific conservative recalibration of the mapper's raw MAPQ.
//! Tables contain exactly one cap for each raw score 0..=60. They can only
//! lower confidence, preserving all search-completeness caps and exact ties.
use std::{
    io::{self, BufRead},
    path::Path,
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MapqCalibration {
    caps: [u8; 61],
}
impl MapqCalibration {
    pub fn from_caps(caps: [u8; 61]) -> Result<Self, String> {
        for (raw, &cap) in caps.iter().enumerate() {
            if cap as usize > raw {
                return Err("calibration cannot increase raw MAPQ".into());
            }
            if raw > 0 && cap < caps[raw - 1] {
                return Err("calibration caps must be monotone".into());
            }
        }
        Ok(Self { caps })
    }
    pub fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let input = io::BufReader::new(std::fs::File::open(path)?);
        let invalid = |s| io::Error::new(io::ErrorKind::InvalidData, s);
        let mut caps = [0; 61];
        let mut next = 0;
        for line in input.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() != 2 {
                return Err(invalid("expected raw_mapq calibrated_cap"));
            }
            let raw: usize = fields[0].parse().map_err(|_| invalid("invalid raw MAPQ"))?;
            let cap: u8 = fields[1].parse().map_err(|_| invalid("invalid MAPQ cap"))?;
            if raw != next || raw > 60 {
                return Err(invalid("expected ordered raw MAPQ rows 0..60"));
            }
            caps[raw] = cap;
            next += 1;
        }
        if next != 61 {
            return Err(invalid("calibration must contain all 61 MAPQ rows"));
        }
        Self::from_caps(caps).map_err(|message| io::Error::new(io::ErrorKind::InvalidData, message))
    }
    pub fn apply(&self, raw: u8) -> u8 {
        self.caps[usize::from(raw.min(60))].min(raw)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn calibration_preserves_caps_and_rejects_confidence_increases() {
        let caps = std::array::from_fn(|q| (q as u8).min(20));
        let table = MapqCalibration::from_caps(caps).unwrap();
        assert_eq!(table.apply(60), 20);
        assert_eq!(table.apply(5), 5);
        assert_eq!(table.apply(0), 0);
        let mut bad = caps;
        bad[0] = 1;
        assert!(MapqCalibration::from_caps(bad).is_err());
        bad = caps;
        bad[30] = 10;
        assert!(MapqCalibration::from_caps(bad).is_err());
    }
}
