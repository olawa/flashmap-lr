//! The fixed local-DP backend used by the RS-LRA default path.
//!
//! This is intentionally a small KSW2 adapter.  It keeps the C library's
//! representation and thread-local mutable aligner behind this module; the
//! rest of RS-LRA only sees a validated neutral CIGAR.

use crate::{types::cigar_edit_distance, Cigar, CigarOp};
use std::cell::RefCell;

const MAX_WINDOW: usize = 8_192;
/// Ceiling for a banded whole-read alignment.
///
/// [`MAX_WINDOW`] bounds a gap fill, where cost grows with the product of the
/// two spans. A banded pass costs span times band, so a read many times longer
/// than that window is still cheap provided the band stays narrow -- which is
/// the whole reason to run one.
const MAX_BANDED_SPAN: usize = 262_144;
const MAX_BANDED_CELLS: usize = 64_000_000;
const MAX_CELLS: usize = 16_000_000;
pub(crate) const MATCH_SCORE: i8 = 2;
pub(crate) const MISMATCH_PENALTY: i8 = 4;
pub(crate) const GAP_OPEN: i8 = 6;
pub(crate) const GAP_EXTEND: i8 = 1;
pub(crate) const GAP_OPEN_DUAL: i8 = 6;
pub(crate) const GAP_EXTEND_DUAL: i8 = 2;
pub(crate) const GAP_OPEN2_DUAL: i8 = 24;
pub(crate) const GAP_EXTEND2_DUAL: i8 = 1;

thread_local! {
    static KSW2_ALIGNER: RefCell<ksw2rs::Aligner> = RefCell::new(ksw2rs::Aligner::new());
    static QUERY_DNA5: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static REFERENCE_DNA5: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalAlignment {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
    pub edit_distance: u32,
}

/// Run the default KSW2 local extension over one bounded query/reference gap.
///
/// The returned coordinates are relative to the supplied slices.  KSW2's
/// extension CIGAR is trimmed of trailing indels, matching FlashMap's LR
/// adapter; callers must still validate that the resulting consumption fits
/// their anchor geometry.
pub fn align_local(query: &[u8], reference: &[u8], band_width: usize) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    if query.len().max(reference.len()) > MAX_WINDOW
        || query.len().saturating_mul(reference.len()) > MAX_CELLS
    {
        return None;
    }

    let band_width = band_width.clamp(1, MAX_WINDOW);
    let (score, raw_cigar) = KSW2_ALIGNER.with(|aligner_cell| {
        QUERY_DNA5.with(|query_cell| {
            REFERENCE_DNA5.with(|reference_cell| {
                let mut query_dna5 = query_cell.borrow_mut();
                let mut reference_dna5 = reference_cell.borrow_mut();
                encode_dna5(query, &mut query_dna5);
                encode_dna5(reference, &mut reference_dna5);

                let matrix = dna5_matrix();
                let input = ksw2rs::Extz2Input {
                    query: &query_dna5,
                    target: &reference_dna5,
                    m: 5,
                    mat: &matrix,
                    q: GAP_OPEN,
                    e: GAP_EXTEND,
                    w: band_width as i32,
                    zdrop: 100,
                    end_bonus: 0,
                    flag: ksw2rs::KSW_EZ_EXTZ_ONLY,
                };
                let mut aligner = aligner_cell.borrow_mut();
                let extension = aligner.align(&input);
                (extension.max, extension.cigar.clone())
            })
        })
    });

    let ops = raw_cigar_to_ops(&raw_cigar);
    if ops.is_empty() {
        return None;
    }
    let cigar = Cigar::new(ops).ok()?;
    let query_consumed = cigar.query_len() as usize;
    let ref_consumed = cigar.reference_len() as usize;
    if query_consumed == 0
        || ref_consumed == 0
        || query_consumed > query.len()
        || ref_consumed > reference.len()
    {
        return None;
    }
    let query_start = query.len() - query_consumed;
    let ref_start = reference.len() - ref_consumed;
    let query_slice = &query[query_start..];
    let ref_slice = &reference[ref_start..];
    let edit_distance = cigar_edit_distance(&cigar, query_slice, ref_slice)?;

    Some(LocalAlignment {
        score: score as i32,
        query_start,
        query_end: query.len(),
        ref_start,
        ref_end: reference.len(),
        cigar,
        edit_distance,
    })
}

/// Run KSW2 dual-affine local extension over one bounded query/reference gap.
pub fn align_local_dual_affine(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
) -> Option<LocalAlignment> {
    align_local_dual_affine_with_scoring(
        query,
        reference,
        band_width,
        GAP_OPEN_DUAL,
        GAP_EXTEND_DUAL,
        GAP_OPEN2_DUAL,
        GAP_EXTEND2_DUAL,
    )
}

/// Run KSW2 dual-affine local extension with explicit scoring penalties.
pub fn align_local_dual_affine_with_scoring(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
    gap_open2: i8,
    gap_extend2: i8,
) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    if query.len().max(reference.len()) > MAX_WINDOW
        || query.len().saturating_mul(reference.len()) > MAX_CELLS
    {
        return None;
    }

    let band_width = band_width.clamp(1, MAX_WINDOW);
    let (score, raw_cigar) = KSW2_ALIGNER.with(|aligner_cell| {
        QUERY_DNA5.with(|query_cell| {
            REFERENCE_DNA5.with(|reference_cell| {
                let mut query_dna5 = query_cell.borrow_mut();
                let mut reference_dna5 = reference_cell.borrow_mut();
                encode_dna5(query, &mut query_dna5);
                encode_dna5(reference, &mut reference_dna5);

                let matrix = dna5_matrix();
                let input = ksw2rs::Extd2Input {
                    query: &query_dna5,
                    target: &reference_dna5,
                    m: 5,
                    mat: &matrix,
                    q: gap_open,
                    e: gap_extend,
                    q2: gap_open2,
                    e2: gap_extend2,
                    w: band_width as i32,
                    zdrop: 100,
                    end_bonus: 0,
                    flag: ksw2rs::KSW_EZ_EXTZ_ONLY,
                };
                let mut aligner = aligner_cell.borrow_mut();
                let extension = aligner.align_extd2(&input);
                (extension.max, extension.cigar.clone())
            })
        })
    });

    let ops = raw_cigar_to_ops(&raw_cigar);
    if ops.is_empty() {
        return None;
    }
    let cigar = Cigar::new(ops).ok()?;
    let query_consumed = cigar.query_len() as usize;
    let ref_consumed = cigar.reference_len() as usize;
    if query_consumed == 0
        || ref_consumed == 0
        || query_consumed > query.len()
        || ref_consumed > reference.len()
    {
        return None;
    }
    let query_start = query.len() - query_consumed;
    let ref_start = reference.len() - ref_consumed;
    let query_slice = &query[query_start..];
    let ref_slice = &reference[ref_start..];
    let edit_distance = cigar_edit_distance(&cigar, query_slice, ref_slice)?;

    Some(LocalAlignment {
        score: score as i32,
        query_start,
        query_end: query.len(),
        ref_start,
        ref_end: reference.len(),
        cigar,
        edit_distance,
    })
}

/// Run the fixed KSW2 affine-gap alignment while requiring both supplied
/// slices to be consumed end-to-end.
///
/// This is the gap-assembly companion to [`align_local`].  It uses the same
/// scoring matrix and gap penalties, but leaves KSW2's extension-only flag
/// off so traceback is anchored at both ends.  Returning `None` on a
/// truncated/soft-clipped traceback keeps callers from silently emitting a
/// CIGAR that does not cover the complete internal gap.
pub fn align_full(query: &[u8], reference: &[u8], band_width: usize) -> Option<LocalAlignment> {
    align_full_with_scoring(query, reference, band_width, GAP_OPEN, GAP_EXTEND)
}

/// Run KSW2 affine-gap alignment with explicit gap opening and extension penalties.
pub fn align_full_with_scoring(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    if query.len().max(reference.len()) > MAX_WINDOW
        || query.len().saturating_mul(reference.len()) > MAX_CELLS
    {
        return None;
    }
    run_extz2(
        query,
        reference,
        band_width.clamp(1, MAX_WINDOW),
        gap_open,
        gap_extend,
    )
}

/// Run KSW2 dual-affine gap alignment consuming both supplied slices end-to-end.
pub fn align_full_dual_affine(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
) -> Option<LocalAlignment> {
    align_full_dual_affine_with_scoring(
        query,
        reference,
        band_width,
        GAP_OPEN_DUAL,
        GAP_EXTEND_DUAL,
        GAP_OPEN2_DUAL,
        GAP_EXTEND2_DUAL,
    )
}

/// Run KSW2 dual-affine gap alignment with explicit gap opening and extension penalties.
pub fn align_full_dual_affine_with_scoring(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
    gap_open2: i8,
    gap_extend2: i8,
) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    if query.len().max(reference.len()) > MAX_WINDOW
        || query.len().saturating_mul(reference.len()) > MAX_CELLS
    {
        return None;
    }
    run_extd2(
        query,
        reference,
        band_width.clamp(1, MAX_WINDOW),
        gap_open,
        gap_extend,
        gap_open2,
        gap_extend2,
    )
}

/// Align a whole read against a known reference window inside a narrow band.
///
/// This is [`align_full`] without the quadratic guard: the caller must have
/// established the diagonal already, so the band -- not the span product --
/// bounds the work. `None` means the span or the band is too large to be
/// worth it, and the caller should fall back to anchor discovery.
pub fn align_banded(query: &[u8], reference: &[u8], band_width: usize) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    let span = query.len().max(reference.len());
    if span > MAX_BANDED_SPAN {
        return None;
    }
    let band = band_width.clamp(1, MAX_WINDOW);
    if span.saturating_mul(band.saturating_mul(2).saturating_add(1)) > MAX_BANDED_CELLS {
        return None;
    }
    run_extz2(query, reference, band, GAP_OPEN, GAP_EXTEND)
}

/// Align a whole read against a known reference window inside a narrow band using dual-affine gap penalties.
pub fn align_banded_dual_affine(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
) -> Option<LocalAlignment> {
    align_banded_dual_affine_with_scoring(
        query,
        reference,
        band_width,
        GAP_OPEN_DUAL,
        GAP_EXTEND_DUAL,
        GAP_OPEN2_DUAL,
        GAP_EXTEND2_DUAL,
    )
}

/// Align a whole read against a known reference window inside a narrow band using dual-affine gap penalties with explicit scoring.
pub fn align_banded_dual_affine_with_scoring(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
    gap_open2: i8,
    gap_extend2: i8,
) -> Option<LocalAlignment> {
    if query.is_empty() || reference.is_empty() {
        return None;
    }
    let span = query.len().max(reference.len());
    if span > MAX_BANDED_SPAN {
        return None;
    }
    let band = band_width.clamp(1, MAX_WINDOW);
    if span.saturating_mul(band.saturating_mul(2).saturating_add(1)) > MAX_BANDED_CELLS {
        return None;
    }
    run_extd2(
        query,
        reference,
        band,
        gap_open,
        gap_extend,
        gap_open2,
        gap_extend2,
    )
}

/// Shared KSW2 extz2 call and CIGAR validation.
///
/// Callers own the size policy: a gap fill and a banded whole-read pass have
/// different reasons to decline, but the alignment itself is the same.
fn run_extz2(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
) -> Option<LocalAlignment> {
    let (score, raw_cigar) = KSW2_ALIGNER.with(|aligner_cell| {
        QUERY_DNA5.with(|query_cell| {
            REFERENCE_DNA5.with(|reference_cell| {
                let mut query_dna5 = query_cell.borrow_mut();
                let mut reference_dna5 = reference_cell.borrow_mut();
                encode_dna5(query, &mut query_dna5);
                encode_dna5(reference, &mut reference_dna5);

                let matrix = dna5_matrix();
                let input = ksw2rs::Extz2Input {
                    query: &query_dna5,
                    target: &reference_dna5,
                    m: 5,
                    mat: &matrix,
                    q: gap_open,
                    e: gap_extend,
                    w: band_width as i32,
                    zdrop: 100,
                    end_bonus: 0,
                    flag: 0,
                };
                let mut aligner = aligner_cell.borrow_mut();
                let extension = aligner.align(&input);
                (extension.score, extension.cigar.clone())
            })
        })
    });

    let ops = raw_cigar_to_ops_full(&raw_cigar);
    if ops.is_empty() || ops.iter().any(|op| matches!(op, CigarOp::SoftClip(_))) {
        return None;
    }
    let cigar = Cigar::new(ops).ok()?;
    if cigar.query_len() as usize != query.len()
        || cigar.reference_len() as usize != reference.len()
    {
        return None;
    }
    let edit_distance = cigar_edit_distance(&cigar, query, reference)?;

    Some(LocalAlignment {
        score,
        query_start: 0,
        query_end: query.len(),
        ref_start: 0,
        ref_end: reference.len(),
        cigar,
        edit_distance,
    })
}

/// Shared KSW2 extd2 dual-affine call and CIGAR validation.
fn run_extd2(
    query: &[u8],
    reference: &[u8],
    band_width: usize,
    gap_open: i8,
    gap_extend: i8,
    gap_open2: i8,
    gap_extend2: i8,
) -> Option<LocalAlignment> {
    let (score, raw_cigar) = KSW2_ALIGNER.with(|aligner_cell| {
        QUERY_DNA5.with(|query_cell| {
            REFERENCE_DNA5.with(|reference_cell| {
                let mut query_dna5 = query_cell.borrow_mut();
                let mut reference_dna5 = reference_cell.borrow_mut();
                encode_dna5(query, &mut query_dna5);
                encode_dna5(reference, &mut reference_dna5);

                let matrix = dna5_matrix();
                let input = ksw2rs::Extd2Input {
                    query: &query_dna5,
                    target: &reference_dna5,
                    m: 5,
                    mat: &matrix,
                    q: gap_open,
                    e: gap_extend,
                    q2: gap_open2,
                    e2: gap_extend2,
                    w: band_width as i32,
                    zdrop: 100,
                    end_bonus: 0,
                    flag: 0,
                };
                let mut aligner = aligner_cell.borrow_mut();
                let extension = aligner.align_extd2(&input);
                (extension.score, extension.cigar.clone())
            })
        })
    });

    let ops = raw_cigar_to_ops_full(&raw_cigar);
    if ops.is_empty() || ops.iter().any(|op| matches!(op, CigarOp::SoftClip(_))) {
        return None;
    }
    let cigar = Cigar::new(ops).ok()?;
    if cigar.query_len() as usize != query.len()
        || cigar.reference_len() as usize != reference.len()
    {
        return None;
    }
    let edit_distance = cigar_edit_distance(&cigar, query, reference)?;

    Some(LocalAlignment {
        score,
        query_start: 0,
        query_end: query.len(),
        ref_start: 0,
        ref_end: reference.len(),
        cigar,
        edit_distance,
    })
}

const DNA5_TABLE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

fn encode_dna5(sequence: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.reserve(sequence.len());
    for &base in sequence {
        output.push(DNA5_TABLE[base as usize]);
    }
}

fn dna5_matrix() -> [i8; 25] {
    let mut matrix = [-MISMATCH_PENALTY; 25];
    for base in 0..5 {
        matrix[base * 5 + base] = MATCH_SCORE;
    }
    // Ambiguous bases neither reward nor penalize an aligned comparison.
    matrix[24] = 0;
    matrix
}

fn raw_cigar_to_ops(raw: &[u32]) -> Vec<CigarOp> {
    raw_cigar_to_ops_with_terminal_policy(raw, true)
}

fn raw_cigar_to_ops_full(raw: &[u32]) -> Vec<CigarOp> {
    raw_cigar_to_ops_with_terminal_policy(raw, false)
}

fn raw_cigar_to_ops_with_terminal_policy(raw: &[u32], trim_terminal_indels: bool) -> Vec<CigarOp> {
    let mut end = raw.len();
    if trim_terminal_indels {
        while end > 0 {
            let op = raw[end - 1] & 0xf;
            if op == 1 || op == 2 {
                end -= 1;
            } else {
                break;
            }
        }
    }

    let mut ops = Vec::new();
    for &packed in &raw[..end] {
        let len = packed >> 4;
        if len == 0 {
            continue;
        }
        match packed & 0xf {
            0 => ops.push(CigarOp::Match(len)),
            1 => ops.push(CigarOp::Ins(len)),
            2 => ops.push(CigarOp::Del(len)),
            3 => ops.push(CigarOp::SoftClip(len)),
            _ => return Vec::new(),
        }
    }
    ops
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna5_encoding_accepts_lowercase_and_n() {
        let mut encoded = Vec::new();
        encode_dna5(b"aCgTNx", &mut encoded);
        assert_eq!(encoded, vec![0, 1, 2, 3, 4, 4]);
    }

    #[test]
    fn raw_cigar_decoder_trims_terminal_indels() {
        assert_eq!(
            raw_cigar_to_ops(&[(3 << 4), (2 << 4) | 1, (1 << 4) | 2]),
            vec![CigarOp::Match(3)]
        );
    }

    #[test]
    fn exact_window_maps_with_fixed_ksw2_scoring() {
        let alignment = align_local(b"ACGTACGT", b"ACGTACGT", 32).unwrap();
        assert_eq!(alignment.cigar.ops(), &[CigarOp::Match(8)]);
        assert_eq!(alignment.edit_distance, 0);
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.ref_start, 0);
    }

    #[test]
    fn mismatch_stays_inside_match_operation_and_counts_in_nm() {
        let alignment = align_local(b"ACGTACGTTCGTACGT", b"ACGTACGTACGTACGT", 32).unwrap();
        assert_eq!(alignment.cigar.ops(), &[CigarOp::Match(16)]);
        assert_eq!(alignment.edit_distance, 1);
    }

    #[test]
    fn short_insertion_is_reported_as_an_indel() {
        let alignment = align_full(b"ACGTACGTTTACGTACGT", b"ACGTACGTACGTACGT", 32).unwrap();
        assert!(alignment
            .cigar
            .ops()
            .iter()
            .any(|op| matches!(op, CigarOp::Ins(2))));
        assert_eq!(alignment.edit_distance, 2);
    }

    #[test]
    fn full_alignment_consumes_both_slices() {
        let alignment = align_full(b"AAA", b"AATAA", 32).unwrap();
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.ref_end, 5);
        assert_eq!(alignment.cigar.query_len(), 3);
        assert_eq!(alignment.cigar.reference_len(), 5);
        assert_eq!(alignment.edit_distance, 2);
    }

    #[test]
    fn full_alignment_keeps_terminal_indels() {
        let alignment = align_full(b"AAAA", b"AAA", 32).unwrap();
        assert_eq!(alignment.cigar.query_len(), 4);
        assert_eq!(alignment.cigar.reference_len(), 3);
        assert!(alignment
            .cigar
            .ops()
            .iter()
            .any(|op| matches!(op, CigarOp::Ins(1))));
    }

    #[test]
    fn dual_affine_local_and_full_alignment() {
        let local = align_local_dual_affine(b"ACGTACGT", b"ACGTACGT", 32).unwrap();
        assert_eq!(local.cigar.ops(), &[CigarOp::Match(8)]);
        assert_eq!(local.edit_distance, 0);

        let full = align_full_dual_affine(b"ACGTACGTTTACGTACGT", b"ACGTACGTACGTACGT", 32).unwrap();
        assert!(full
            .cigar
            .ops()
            .iter()
            .any(|op| matches!(op, CigarOp::Ins(2))));
        assert_eq!(full.edit_distance, 2);

        let banded = align_banded_dual_affine(b"ACGTACGT", b"ACGTACGT", 32).unwrap();
        assert_eq!(banded.cigar.ops(), &[CigarOp::Match(8)]);
        assert_eq!(banded.edit_distance, 0);
    }
}

#[cfg(test)]
mod band_tests {
    use super::*;

    /// KSW2 restricts the alignment to `|i - j| <= band`, so a band narrower
    /// than the length difference has no path that spells the indel. The
    /// medium gap stage used to clamp its band at 256 while admitting a delta
    /// of up to 512, which is this test's gap.
    #[test]
    fn a_band_narrower_than_the_indel_cannot_spell_it() {
        let unit = b"GATTACACGCTAGCTTACGGTCAAGCTTGAC";
        let reference: Vec<u8> = unit.iter().cycle().take(1_000).copied().collect();
        let mut query = reference[..300].to_vec();
        query.extend_from_slice(&reference[700..1_000]);

        let deletion = |band: usize| -> Option<u32> {
            Some(
                align_full(&query, &reference, band)?
                    .cigar
                    .ops()
                    .iter()
                    .filter_map(|op| match op {
                        crate::CigarOp::Del(len) => Some(*len),
                        _ => None,
                    })
                    .sum(),
            )
        };

        // Not a truncated alignment: no alignment at all, so the stage that
        // asked for it is skipped and the gap falls through to a coarser one.
        assert_eq!(deletion(256), None, "a 256 band has no path to 400 bases");
        assert_eq!(deletion(432), Some(400), "a band above the delta has one");
    }
}
