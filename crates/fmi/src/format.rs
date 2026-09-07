//! Fixed v14 format constants and collision-table validation shared by writer
//! and readers. Each collision record is two little-endian u64s, never a
//! native Rust tuple cast to bytes.
pub const VERSION: u32 = 14;
pub const PRIMARY_COLLISIONS: u32 = 19;
pub const SECONDARY_COLLISIONS: u32 = 20;
pub const SJDB_COLLISIONS: u32 = 21;
pub fn decode_collisions(bytes: &[u8], range_count: usize) -> Option<Vec<(u64, u64)>> {
    if !bytes.len().is_multiple_of(16) {
        return None;
    }
    let mut pairs = Vec::with_capacity(bytes.len() / 16);
    for chunk in bytes.chunks_exact(16) {
        let index = u64::from_le_bytes(chunk[..8].try_into().ok()?);
        let code = u64::from_le_bytes(chunk[8..].try_into().ok()?);
        if index >= range_count as u64 || pairs.last().is_some_and(|p: &(u64, u64)| p.0 >= index) {
            return None;
        }
        pairs.push((index, code));
    }
    Some(pairs)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_truncated_duplicate_unordered_and_out_of_bounds_records() {
        let encode = |pairs: &[(u64, u64)]| {
            pairs
                .iter()
                .flat_map(|&(i, c)| i.to_le_bytes().into_iter().chain(c.to_le_bytes()))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            decode_collisions(&encode(&[(0, 9), (2, 7)]), 3),
            Some(vec![(0, 9), (2, 7)])
        );
        assert!(decode_collisions(&[0; 15], 3).is_none());
        assert!(decode_collisions(&encode(&[(0, 9), (0, 7)]), 3).is_none());
        assert!(decode_collisions(&encode(&[(2, 9), (1, 7)]), 3).is_none());
        assert!(decode_collisions(&encode(&[(3, 9)]), 3).is_none());
    }
}
