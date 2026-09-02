//! What to store for a seed that occurs more often than the frequency cap.
//!
//! The cap decides which reference positions a query can ever reach, so this
//! is the one build decision a mapper cannot work around: a position stored
//! nowhere has to be rediscovered from the reference, or missed. Ported from
//! FlashMap's `index::cap_policy` so the writer and the reader agree about
//! what a capped range means.

/// Occurrences to retain for one seed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KeepPlan {
    /// Store every occurrence; the seed is under the cap.
    All,
    /// Store none. The mapper still learns the true multiplicity, so it knows
    /// the seed is repetitive rather than absent, but gets no positions to
    /// chain.
    None,
    /// Store the first `n` occurrences in sorted order.
    First(usize),
    /// Store `n` occurrences spread evenly across the group.
    Spaced(usize),
}

/// How a build treats a seed above the frequency cap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeedCapPolicy {
    /// Keep the first `max_freq` occurrences. The sort is by `(ref_id, pos)`,
    /// so these cluster on the lowest contig rather than representing the
    /// repeat family.
    First,
    /// Store no positions for an over-cap seed.
    Drop,
    /// Keep `max_freq` occurrences spread evenly across the group, so the
    /// stored subset spans the whole repeat family.
    SampleSpaced,
    /// Graded by how repetitive the seed is: keep everything up to `medium`,
    /// spaced-sample up to `high`, store nothing above.
    ///
    /// This exists because the two ends of the distribution want opposite
    /// treatment. A seed with 20 copies keeps 16 of them, so the copy a read
    /// needs usually survives; a seed with 20000 keeps 16, and the sample is
    /// noise. One threshold cannot serve both.
    Adaptive { medium: usize, high: usize },
}

/// A cap policy could not be parsed from its name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnknownCapPolicy(pub String);

impl std::fmt::Display for UnknownCapPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "unknown cap policy {:?}; expected first, drop, sample, sample-spaced or adaptive",
            self.0
        )
    }
}

impl std::error::Error for UnknownCapPolicy {}

impl SeedCapPolicy {
    /// Parse a policy name together with its optional thresholds.
    ///
    /// `medium` and `high` are only meaningful for `adaptive`; without them it
    /// degenerates to `sample-spaced`, which is the least surprising reading
    /// of "adapt, but I gave you nothing to adapt to".
    pub fn parse(
        name: &str,
        max_freq: usize,
        medium: Option<usize>,
        high: Option<usize>,
    ) -> Result<Self, UnknownCapPolicy> {
        Ok(match name {
            "first" => Self::First,
            "drop" => Self::Drop,
            // `sample` was accepted historically and never distinguished from
            // the spaced variant; keep taking it rather than breaking scripts.
            "sample" | "sample-spaced" => Self::SampleSpaced,
            "adaptive" => Self::Adaptive {
                medium: medium.unwrap_or(max_freq),
                high: high.unwrap_or(usize::MAX),
            },
            other => return Err(UnknownCapPolicy(other.to_owned())),
        })
    }

    /// The name this policy is recorded under in index metadata.
    pub fn name(&self) -> &'static str {
        match self {
            Self::First => "first",
            Self::Drop => "drop",
            Self::SampleSpaced => "sample-spaced",
            Self::Adaptive { .. } => "adaptive",
        }
    }

    /// Which occurrences to store for a seed seen `count` times.
    pub fn plan(&self, count: usize, max_freq: usize) -> KeepPlan {
        if count <= max_freq {
            return KeepPlan::All;
        }
        match *self {
            Self::First => KeepPlan::First(max_freq),
            Self::Drop => KeepPlan::None,
            Self::SampleSpaced => KeepPlan::Spaced(max_freq),
            Self::Adaptive { medium, high } => {
                if count <= medium {
                    KeepPlan::All
                } else if count <= high {
                    KeepPlan::Spaced(max_freq)
                } else {
                    KeepPlan::None
                }
            }
        }
    }
}

/// Offsets within a group of `count` occurrences, per the plan.
///
/// Spaced sampling steps by `count / n` so the kept copies span the whole
/// group; `First` returns a prefix. Both are deterministic, which index
/// reproducibility depends on.
pub fn selected_offsets(plan: KeepPlan, count: usize) -> Vec<usize> {
    match plan {
        KeepPlan::All => (0..count).collect(),
        KeepPlan::None => Vec::new(),
        KeepPlan::First(n) => (0..n.min(count)).collect(),
        KeepPlan::Spaced(n) => {
            let n = n.min(count);
            if n == 0 {
                return Vec::new();
            }
            let step = count as f64 / n as f64;
            let mut out = Vec::with_capacity(n);
            let mut last = usize::MAX;
            for index in 0..n {
                let offset = ((index as f64 * step) as usize).min(count - 1);
                // Guard against a repeated index when step rounds below 1.
                let offset = if offset == last { offset + 1 } else { offset };
                if offset >= count {
                    break;
                }
                out.push(offset);
                last = offset;
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_seed_under_the_cap_is_stored_whole_by_every_policy() {
        for policy in [
            SeedCapPolicy::First,
            SeedCapPolicy::Drop,
            SeedCapPolicy::SampleSpaced,
            SeedCapPolicy::Adaptive {
                medium: 16,
                high: 128,
            },
        ] {
            assert_eq!(policy.plan(16, 16), KeepPlan::All, "{policy:?}");
        }
    }

    #[test]
    fn adaptive_grades_by_how_repetitive_the_seed_is() {
        let policy = SeedCapPolicy::Adaptive {
            medium: 32,
            high: 128,
        };
        assert_eq!(policy.plan(20, 16), KeepPlan::All);
        assert_eq!(policy.plan(64, 16), KeepPlan::Spaced(16));
        assert_eq!(policy.plan(5000, 16), KeepPlan::None);
    }

    #[test]
    fn a_spaced_sample_spans_the_group_rather_than_its_start() {
        let offsets = selected_offsets(KeepPlan::Spaced(4), 100);
        assert_eq!(offsets, vec![0, 25, 50, 75]);
        // The prefix policy is what this exists to avoid: sorted by
        // (ref_id, pos), it takes every copy from the lowest contig.
        assert_eq!(selected_offsets(KeepPlan::First(4), 100), vec![0, 1, 2, 3]);
    }

    #[test]
    fn a_sample_never_repeats_or_overruns_its_group() {
        for count in 1..200usize {
            for keep in 1..40usize {
                let offsets = selected_offsets(KeepPlan::Spaced(keep), count);
                assert!(offsets.len() <= keep.min(count));
                assert!(offsets.iter().all(|offset| *offset < count));
                assert!(
                    offsets.windows(2).all(|pair| pair[0] < pair[1]),
                    "count {count} keep {keep} gave {offsets:?}"
                );
            }
        }
    }

    #[test]
    fn adaptive_without_thresholds_is_a_spaced_sample() {
        let policy = SeedCapPolicy::parse("adaptive", 16, None, None).unwrap();
        assert_eq!(policy.plan(1000, 16), KeepPlan::Spaced(16));
    }

    #[test]
    fn a_policy_name_round_trips() {
        for name in ["first", "drop", "sample-spaced", "adaptive"] {
            let policy = SeedCapPolicy::parse(name, 16, None, None).unwrap();
            assert_eq!(policy.name(), name);
        }
        assert!(SeedCapPolicy::parse("nonsense", 16, None, None).is_err());
    }
}
