//! SAM optional-field validation and normalization.
//!
//! FASTA/FASTQ header comments are whitespace-delimited, while SAM optional
//! fields are tab-delimited `TAG:TYPE:VALUE` tokens. The parser therefore keeps
//! only syntactically valid SAM fields and normalizes their separators.
//! In particular, MM/ML tags are preserved verbatim: SAM defines their
//! coordinates in the as-sequenced orientation, so they must not be reversed
//! merely because the SAM record carries FLAG 0x10.

/// Return valid SAM optional fields from a FASTA/FASTQ comment, separated by
/// tabs. Free-form comment words and malformed fields are ignored so that a
/// normal FASTQ description can never become an invalid SAM optional field.
///
/// The first occurrence of a tag wins. SAM records may contain each tag at
/// most once; keeping the first occurrence also makes this helper safe for
/// callers that pass already-normalized input more than once.
pub fn normalize_optional_fields(tags: &str) -> String {
    normalize_optional_fields_excluding(tags, &[])
}

/// As [`normalize_optional_fields`], but omit fields whose two-character tag
/// name is listed in `excluded`. SAM writers use this for fields such as NM
/// and AS that are generated from the computed alignment.
pub(crate) fn normalize_optional_fields_excluding(tags: &str, excluded: &[&str]) -> String {
    let mut seen = std::collections::HashSet::new();
    tags.split_whitespace()
        .filter_map(|token| {
            let (name, _, _) = parse_optional_field(token)?;
            if excluded.contains(&name) || !seen.insert(name) {
                return None;
            }
            Some(token)
        })
        .collect::<Vec<_>>()
        .join("\t")
}

fn parse_optional_field(token: &str) -> Option<(&str, char, &str)> {
    let mut parts = token.splitn(3, ':');
    let name = parts.next()?;
    let type_field = parts.next()?;
    let value = parts.next()?;
    if name.len() != 2
        || !name.as_bytes()[0].is_ascii_alphabetic()
        || !name.as_bytes()[1].is_ascii_alphanumeric()
        || type_field.len() != 1
    {
        return None;
    }
    let field_type = type_field.as_bytes()[0] as char;
    is_valid_value(field_type, value).then_some((name, field_type, value))
}

fn is_valid_value(field_type: char, value: &str) -> bool {
    match field_type {
        // SAM's A type is one printable ASCII character.
        'A' => value.len() == 1 && value.as_bytes()[0].is_ascii_graphic(),
        // Rust's signed integer parser accepts the decimal representation used
        // by SAM and rejects whitespace and trailing characters.
        'i' => !value.is_empty() && value.parse::<i64>().is_ok(),
        // Reject NaN and infinities even though Rust's parser accepts them;
        // SAM f values are finite decimal floating-point values.
        'f' => value
            .parse::<f64>()
            .map(|number| number.is_finite())
            .unwrap_or(false),
        // Z values cannot contain SAM field separators or control bytes. A
        // space cannot arrive through a FASTX comment token, but is accepted
        // here for callers that provide a pre-tokenized value.
        'Z' => value.bytes().all(|byte| (0x20..=0x7e).contains(&byte)),
        // H is an even-length hexadecimal string. Empty H values are allowed
        // by the general SAM grammar and are harmless to pass through.
        'H' => value.len().is_multiple_of(2) && value.bytes().all(|byte| byte.is_ascii_hexdigit()),
        'B' => is_valid_array(value),
        _ => false,
    }
}

fn is_valid_array(value: &str) -> bool {
    let Some((subtype, values)) = value.split_once(',') else {
        return false;
    };
    if subtype.len() != 1 || values.is_empty() {
        return false;
    }
    values
        .split(',')
        .all(|value| match subtype.as_bytes()[0] as char {
            'c' => value.parse::<i8>().is_ok(),
            'C' => value.parse::<u8>().is_ok(),
            's' => value.parse::<i16>().is_ok(),
            'S' => value.parse::<u16>().is_ok(),
            'i' => value.parse::<i32>().is_ok(),
            'I' => value.parse::<u32>().is_ok(),
            'f' => value
                .parse::<f64>()
                .map(|number| number.is_finite())
                .unwrap_or(false),
            _ => false,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_valid_fields_and_drops_free_form_comment() {
        let tags = "description MM:Z:C+m?,1,2; ML:B:C,100,200 RG:Z:sample  bad:x:y";
        assert_eq!(
            normalize_optional_fields(tags),
            "MM:Z:C+m?,1,2;\tML:B:C,100,200\tRG:Z:sample"
        );
    }

    #[test]
    fn accepts_common_sam_types_and_rejects_malformed_values() {
        let valid = "NM:i:3 AS:i:-12 ts:A:+ XX:f:1.5 HX:H:deadbeef BQ:B:C,0,255";
        assert_eq!(normalize_optional_fields(valid), valid.replace(' ', "\t"));

        let invalid = "N:i:3 N_:i:3 NM:x:3 NM:i:nope XX:A:ab HX:H:abc BQ:B:C,999";
        assert!(normalize_optional_fields(invalid).is_empty());
    }

    #[test]
    fn keeps_first_duplicate_and_can_exclude_generated_fields() {
        let tags = "RG:Z:first RG:Z:second NM:i:99 AS:i:1 CB:Z:cell";
        assert_eq!(
            normalize_optional_fields(tags),
            "RG:Z:first\tNM:i:99\tAS:i:1\tCB:Z:cell"
        );
        assert_eq!(
            normalize_optional_fields_excluding(tags, &["NM", "AS"]),
            "RG:Z:first\tCB:Z:cell"
        );
    }

    #[test]
    fn preserves_multitrack_mm_ml_without_orientation_changes() {
        let tags = "MM:Z:A+a.,0;C+m?,1; ML:B:C,255,200 RG:Z:cell1";
        assert_eq!(normalize_optional_fields(tags), tags.replace(' ', "\t"));
    }
}
