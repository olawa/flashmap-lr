//! Assembly-scoped DP reuse. Entries only refer to subslices of the two
//! immutable sequences borrowed by Scope; addresses never escape that scope.
use crate::{config::ScoringPolicy, LocalAlignment};
use std::{cell::RefCell, marker::PhantomData};

#[derive(Clone, PartialEq, Eq)]
struct Key {
    q: usize,
    qlen: usize,
    r: usize,
    rlen: usize,
    band: usize,
    scoring: ScoringPolicy,
}
struct Cache {
    query: (usize, usize),
    reference: (usize, usize),
    entries: Vec<(Key, Result<LocalAlignment, crate::dp::DpFailure>)>,
    next: usize,
    recording: bool,
}
thread_local! { static CACHE: RefCell<Option<Cache>> = const { RefCell::new(None) }; }
pub(super) struct Scope<'a> {
    previous: Option<Cache>,
    _borrow: PhantomData<&'a [u8]>,
}
impl<'a> Scope<'a> {
    pub(super) fn new(query: &'a [u8], reference: &'a [u8]) -> Self {
        let cache = Cache {
            query: bounds(query),
            reference: bounds(reference),
            entries: Vec::new(),
            next: 0,
            recording: true,
        };
        Self {
            previous: CACHE.with(|slot| slot.replace(Some(cache))),
            _borrow: PhantomData,
        }
    }
}
impl Scope<'_> {
    /// Only repair probes populate the cache; the ordinary assembly consumes
    /// their results without cloning every unrelated one-off gap.
    pub(super) fn finish_recording(&self) {
        CACHE.with(|slot| {
            if let Some(cache) = slot.borrow_mut().as_mut() {
                cache.recording = false;
            }
        });
    }
}
impl Drop for Scope<'_> {
    fn drop(&mut self) {
        CACHE.with(|slot| {
            slot.replace(self.previous.take());
        });
    }
}
fn bounds(s: &[u8]) -> (usize, usize) {
    (s.as_ptr() as usize, s.len())
}
fn contained(s: &[u8], parent: (usize, usize)) -> bool {
    (s.as_ptr() as usize)
        .checked_sub(parent.0)
        .is_some_and(|offset| offset <= parent.1 && s.len() <= parent.1 - offset)
}
/// The flags, z-drop, and full-alignment algorithm are fixed by align_full;
/// scoring and band are explicit parts of the key. Failures are reusable only
/// for precisely the same request, never for another band's search space.
pub(super) fn align(
    query: &[u8],
    reference: &[u8],
    band: usize,
    scoring: &ScoringPolicy,
) -> (bool, Result<LocalAlignment, crate::dp::DpFailure>) {
    let key = Key {
        q: query.as_ptr() as usize,
        qlen: query.len(),
        r: reference.as_ptr() as usize,
        rlen: reference.len(),
        band,
        scoring: *scoring,
    };
    let cached = CACHE.with(|slot| {
        let borrow = slot.borrow();
        let cache = borrow.as_ref()?;
        if !contained(query, cache.query) || !contained(reference, cache.reference) {
            return None;
        }
        cache
            .entries
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, result)| result.clone())
    });
    if let Some(result) = cached {
        return (true, result);
    }
    let result = crate::dp::align_full_outcome(query, reference, band, scoring);
    // Bound retained traceback output independently of sequence length.
    if result
        .as_ref()
        .map_or(true, |a| a.cigar.ops().len() <= 4096)
    {
        CACHE.with(|slot| {
            let mut borrow = slot.borrow_mut();
            let Some(cache) = borrow.as_mut() else {
                return;
            };
            if !cache.recording
                || !contained(query, cache.query)
                || !contained(reference, cache.reference)
            {
                return;
            }
            let entry = (key, result.clone());
            if cache.entries.len() < 32 {
                cache.entries.push(entry);
            } else {
                cache.entries[cache.next] = entry;
                cache.next = (cache.next + 1) % 32;
            }
        });
    }
    (false, result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scope_reuses_only_identical_immutable_intervals_and_policy() {
        let policy = crate::config::ResolvedMapperPolicy::from_mapper_config(
            &crate::MapperConfig::default(),
        )
        .unwrap()
        .scoring;
        let q = b"AACCGGTTAACCGGTT";
        let r = b"AACCGGTAACCGGTT";
        let result;
        {
            let _scope = Scope::new(q, r);
            let first = align(q, r, 8, &policy);
            assert!(!first.0);
            assert!(first.1.is_ok());
            let second = align(q, r, 8, &policy);
            assert!(second.0);
            assert_eq!(first.1, second.1);
            assert!(!align(q, r, 9, &policy).0);
            let mut changed = policy;
            changed.dual_affine = true;
            changed.gap_open2 = 24;
            changed.gap_extend2 = 1;
            assert!(!align(q, r, 8, &changed).0);
            assert!(!align(&q[1..], r, 8, &policy).0);
            result = first.1;
        }
        let outside = align(q, r, 8, &policy);
        assert!(!outside.0);
        assert_eq!(outside.1, result);
    }
}
