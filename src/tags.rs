//! SAM optional tag and base-modification (methylation `MM`/`ML`) processing.
//!
//! Standard PacBio HiFi and Oxford Nanopore FASTQ headers carry SAM optional
//! tags in the header comment, e.g.:
//! `@read_name MM:Z:C+m?,1,7,...; ML:B:C,250,180,... RG:Z:sample1`
//!
//! When a read maps to the forward strand, tags are emitted directly.
//! When a read maps to the reverse strand (FLAG 0x10):
//! - `MM:Z` / `Mm:Z` modification strand is inverted (`+` <-> `-`), and deltas are
//!   recalculated for the 5'->3' direction of the reverse-complemented sequence.
//! - `ML:B:C` / `Ml:B:C` probability arrays are reversed within each modification track.
//! - Other SAM tags (`RG:Z`, `CB:Z`, etc.) are preserved verbatim.

/// Transform SAM optional tags from a FASTQ comment for a reverse-strand alignment.
pub fn transform_tags_for_reverse_strand(tags: &str, forward_sequence: &[u8]) -> String {
    let tokens: Vec<&str> = tags.split_whitespace().collect();
    if tokens.is_empty() {
        return String::new();
    }

    // Find MM and ML tags if present
    let mm_idx = tokens
        .iter()
        .position(|t| t.starts_with("MM:Z:") || t.starts_with("Mm:Z:"));
    let ml_idx = tokens
        .iter()
        .position(|t| t.starts_with("ML:B:C,") || t.starts_with("Ml:B:C,"));

    let (transformed_mm, transformed_ml) = match mm_idx {
        Some(m_i) => {
            let mm_tag = tokens[m_i];
            let ml_tag = ml_idx.map(|l_i| tokens[l_i]);
            reverse_base_modifications(forward_sequence, mm_tag, ml_tag)
        }
        None => (None, None),
    };

    let mut result_parts = Vec::with_capacity(tokens.len());
    for (i, token) in tokens.into_iter().enumerate() {
        if Some(i) == mm_idx {
            if let Some(ref new_mm) = transformed_mm {
                result_parts.push(new_mm.as_str());
            } else {
                result_parts.push(token);
            }
        } else if Some(i) == ml_idx {
            if let Some(ref new_ml) = transformed_ml {
                result_parts.push(new_ml.as_str());
            } else {
                result_parts.push(token);
            }
        } else {
            result_parts.push(token);
        }
    }

    result_parts.join("\t")
}

/// Reverse `MM:Z` and `ML:B:C` tags for reverse-complemented sequence.
fn reverse_base_modifications(
    forward_sequence: &[u8],
    mm_tag: &str,
    ml_tag: Option<&str>,
) -> (Option<String>, Option<String>) {
    let tag_prefix = if mm_tag.starts_with("Mm:Z:") {
        "Mm:Z:"
    } else {
        "MM:Z:"
    };
    let mm_payload = &mm_tag[tag_prefix.len()..];

    // Parse ML probabilities if present
    let ml_values: Vec<u8> = ml_tag
        .and_then(|ml| {
            let payload = ml
                .strip_prefix("ML:B:C,")
                .or_else(|| ml.strip_prefix("Ml:B:C,"))?;
            Some(
                payload
                    .split(',')
                    .filter_map(|s| s.trim().parse::<u8>().ok())
                    .collect(),
            )
        })
        .unwrap_or_default();

    let mut ml_cursor = 0usize;
    let mut reversed_tracks = Vec::new();
    let mut reversed_ml_all = Vec::new();

    // MM tracks are separated by semicolons
    for track in mm_payload.split(';').filter(|t| !t.trim().is_empty()) {
        let mut parts = track.split(',');
        let Some(header) = parts.next() else {
            continue;
        };
        if header.len() < 2 {
            reversed_tracks.push(track.to_string());
            continue;
        }

        let canonical_base = header.as_bytes()[0];
        let strand_char = header.as_bytes()[1] as char;
        let mod_and_flag = &header[2..];

        let deltas: Vec<usize> = parts
            .filter_map(|p| p.trim().parse::<usize>().ok())
            .collect();
        let num_mods = deltas.len();

        let track_ml = if ml_cursor + num_mods <= ml_values.len() {
            let slice = &ml_values[ml_cursor..ml_cursor + num_mods];
            ml_cursor += num_mods;
            Some(slice)
        } else {
            None
        };

        if num_mods == 0 {
            let new_strand = if strand_char == '+' { '-' } else { '+' };
            reversed_tracks.push(format!(
                "{}{}{}",
                canonical_base as char, new_strand, mod_and_flag
            ));
            continue;
        }

        // Count total occurrences of canonical base in original forward sequence
        let total_canonical = forward_sequence
            .iter()
            .filter(|&&b| b.eq_ignore_ascii_case(&canonical_base))
            .count();

        // Compute 0-indexed positions of modified bases among the canonical base occurrences
        let mut k_positions = Vec::with_capacity(num_mods);
        let mut cur = 0usize;
        for (idx, &d) in deltas.iter().enumerate() {
            if idx == 0 {
                cur = d;
            } else {
                cur = cur.saturating_add(1).saturating_add(d);
            }
            k_positions.push(cur);
        }

        // Compute new deltas on reverse strand
        let last_k = *k_positions.last().unwrap_or(&0);
        let mut new_deltas = Vec::with_capacity(num_mods);
        // First delta on reverse strand: count of unmodified canonical bases after the last mod on forward
        new_deltas.push(total_canonical.saturating_sub(1).saturating_sub(last_k));

        for i in (1..k_positions.len()).rev() {
            let d_rev = k_positions[i]
                .saturating_sub(k_positions[i - 1])
                .saturating_sub(1);
            new_deltas.push(d_rev);
        }

        let new_strand = if strand_char == '+' { '-' } else { '+' };
        let new_header = format!("{}{}{}", canonical_base as char, new_strand, mod_and_flag);
        let new_deltas_str = new_deltas
            .iter()
            .map(|d| d.to_string())
            .collect::<Vec<_>>()
            .join(",");
        reversed_tracks.push(format!("{new_header},{new_deltas_str}"));

        if let Some(ml_slice) = track_ml {
            for &p in ml_slice.iter().rev() {
                reversed_ml_all.push(p);
            }
        }
    }

    let new_mm = format!("{}{};", tag_prefix, reversed_tracks.join(";"));
    let new_ml = if !reversed_ml_all.is_empty() {
        let prefix = if ml_tag.map_or(false, |t| t.starts_with("Ml:B:C,")) {
            "Ml:B:C,"
        } else {
            "ML:B:C,"
        };
        Some(format!(
            "{}{}",
            prefix,
            reversed_ml_all
                .iter()
                .map(|p| p.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ))
    } else {
        ml_tag.map(|s| s.to_string())
    };

    (Some(new_mm), new_ml)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_single_mm_ml_track() {
        // Sequence with 6 'C's: indices 1, 3, 4, 5, 6, 7 in ACACCCCCAA
        let seq = b"ACACCCCCAA";
        // Mods at C-index 1 (2nd C) and C-index 4 (5th C)
        // deltas: 1, 2. (k0=1, k1=4). total C=6.
        // On reverse: total C=6. k'_0 = 6 - 1 - 4 = 1. k'_1 = 4 - 1 - 1 = 2.
        // deltas_rev: 1, 2.
        let tags = "MM:Z:C+m?,1,2; ML:B:C,100,200 RG:Z:sample1";
        let transformed = transform_tags_for_reverse_strand(tags, seq);
        assert!(transformed.contains("MM:Z:C-m?,1,2;"));
        assert!(transformed.contains("ML:B:C,200,100"));
        assert!(transformed.contains("RG:Z:sample1"));
    }

    #[test]
    fn test_round_trip_reversal() {
        let seq = b"ACACCCCCAA";
        let tags = "MM:Z:C+m?,1,2; ML:B:C,100,200";
        let rev1 = transform_tags_for_reverse_strand(tags, seq);
        let rev2 = transform_tags_for_reverse_strand(&rev1, seq);
        assert_eq!(rev2, "MM:Z:C+m?,1,2;\tML:B:C,100,200");
    }

    #[test]
    fn test_multi_track_pacbio_tags() {
        let seq = b"AACCGGTTAACCGGTT";
        let tags = "MM:Z:A+a.,0;C+m?,1; ML:B:C,255,200 RG:Z:cell1";
        let transformed = transform_tags_for_reverse_strand(tags, seq);
        assert!(transformed.contains("MM:Z:A-a.,3;C-m?,2;"));
        assert!(transformed.contains("ML:B:C,255,200"));
        assert!(transformed.contains("RG:Z:cell1"));
    }
}
