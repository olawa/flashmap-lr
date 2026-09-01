//! Read segmentation shared by the LR seeding phases.

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Segment {
    pub index: usize,
    pub read_start: usize,
    pub read_end: usize,
}

impl Segment {
    pub fn len(&self) -> usize {
        self.read_end.saturating_sub(self.read_start)
    }

    pub fn is_empty(&self) -> bool {
        self.read_start >= self.read_end
    }
}

/// Partition a read into overlapping windows without leaving uncovered bases.
pub fn segment_read(sequence: &[u8], segment_size: usize, overlap: usize) -> Vec<Segment> {
    if sequence.is_empty() || segment_size == 0 {
        return Vec::new();
    }
    if sequence.len() <= segment_size {
        return vec![Segment {
            index: 0,
            read_start: 0,
            read_end: sequence.len(),
        }];
    }

    let overlap = overlap.min(segment_size.saturating_sub(1));
    let step = segment_size - overlap;
    let mut segments = Vec::new();
    let mut start = 0;

    while start < sequence.len() {
        let end = start.saturating_add(segment_size);
        if end >= sequence.len() {
            let adjusted_start = sequence.len() - segment_size;
            if segments
                .last()
                .is_some_and(|segment: &Segment| adjusted_start <= segment.read_start)
            {
                break;
            }
            segments.push(Segment {
                index: segments.len(),
                read_start: adjusted_start,
                read_end: sequence.len(),
            });
            break;
        }

        segments.push(Segment {
            index: segments.len(),
            read_start: start,
            read_end: end,
        });
        start = start.saturating_add(step);
    }

    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn covers_the_entire_read() {
        let sequence = vec![b'A'; 5000];
        let segments = segment_read(&sequence, 2048, 256);
        assert!(!segments.is_empty());
        let mut covered = vec![false; sequence.len()];
        for segment in &segments {
            assert!(!segment.is_empty());
            assert_eq!(segment.len(), segment.read_end - segment.read_start);
            for base in &mut covered[segment.read_start..segment.read_end] {
                *base = true;
            }
        }
        assert!(covered.into_iter().all(|covered| covered));
    }

    #[test]
    fn overlap_is_preserved() {
        let sequence = vec![b'A'; 4000];
        let segments = segment_read(&sequence, 2048, 256);
        assert_eq!(segments.len(), 3);
        assert_eq!(segments[0].read_end - segments[1].read_start, 256);
    }

    #[test]
    fn short_and_empty_reads_are_handled() {
        assert_eq!(segment_read(b"ACGT", 2048, 256).len(), 1);
        assert!(segment_read(b"", 2048, 256).is_empty());
        assert!(segment_read(b"ACGT", 0, 0).is_empty());
    }
}
