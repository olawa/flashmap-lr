//! Small DNA primitives shared by seeding and alignment refinement.
//!
//! These helpers deliberately implement only unambiguous DNA. Ambiguous
//! bases do not encode as k-mers and compare as mismatches in the HiFi path.

/// Encode one unambiguous DNA base using the two-bit A/C/G/T alphabet.
pub(crate) fn base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Encode an unambiguous DNA k-mer of at most 32 bases.
pub(crate) fn encode_kmer(sequence: &[u8]) -> Option<u64> {
    if sequence.is_empty() || sequence.len() > 32 {
        return None;
    }
    let mut code = 0u64;
    for &base in sequence {
        code = (code << 2) | u64::from(base_code(base)?);
    }
    Some(code)
}

/// Count case-insensitive substitutions over the shared span.
pub(crate) fn mismatch_count(query: &[u8], reference: &[u8]) -> usize {
    query
        .iter()
        .zip(reference)
        .filter(|(query, reference)| !query.eq_ignore_ascii_case(reference))
        .count()
}

/// Case-insensitive substitution rate over the shared span.
pub(crate) fn mismatch_rate(query: &[u8], reference: &[u8]) -> f64 {
    let length = query.len().min(reference.len());
    if length == 0 {
        return 1.0;
    }
    mismatch_count(&query[..length], &reference[..length]) as f64 / length as f64
}

pub(crate) fn bases_equal(query: u8, reference: u8) -> bool {
    query.eq_ignore_ascii_case(&reference)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_encoding_is_case_insensitive_and_rejects_ambiguous_bases() {
        assert_eq!(encode_kmer(b"ACGT"), encode_kmer(b"acgt"));
        assert_eq!(encode_kmer(b"ACGT"), Some(0b00_01_10_11));
        assert_eq!(encode_kmer(b""), None);
        assert_eq!(encode_kmer(b"ACNT"), None);
    }

    #[test]
    fn mismatch_count_treats_case_as_equal_and_ambiguity_as_a_mismatch() {
        assert_eq!(mismatch_count(b"AcgN", b"ACGT"), 1);
        assert_eq!(mismatch_rate(b"AcgN", b"ACGT"), 0.25);
    }
}
