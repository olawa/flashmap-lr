//! The fixed RS-LRA worker-pool runner.
//!
//! This module mirrors the LR execution shape used by FlashMap, but keeps the
//! execution boundary independent of the aligner itself:
//!
//! ```text
//! source iterator -> bounded owned batches -> mapper workers -> ordered sink
//! ```
//!
//! A caller injects only the per-read mapper and the output sink.  In
//! particular, this module does not know about SAM/BAM, indexes, or the LR
//! algorithm.  The implementation uses `std::sync::mpsc::sync_channel`, so a
//! standalone RS-LRA build does not need a channel/runtime dependency.

use std::collections::{BTreeMap, VecDeque};
use std::convert::Infallible;
use std::fmt;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::thread;

use crate::config::WorkerPoolConfig;

type RawQueue<R, E> = Arc<BoundedQueue<Result<ReadBatch<R>, E>>>;
type MappedQueue<T, SE, ME> = Arc<BoundedQueue<Result<MappedBatch<T>, InternalFailure<SE, ME>>>>;

struct BoundedQueue<T> {
    state: Mutex<QueueState<T>>,
    not_empty: Condvar,
    not_full: Condvar,
    capacity: usize,
}

struct QueueState<T> {
    items: VecDeque<T>,
    closed: bool,
}

impl<T> BoundedQueue<T> {
    fn new(capacity: usize) -> Self {
        Self {
            state: Mutex::new(QueueState {
                items: VecDeque::with_capacity(capacity),
                closed: false,
            }),
            not_empty: Condvar::new(),
            not_full: Condvar::new(),
            capacity: capacity.max(1),
        }
    }

    fn push(&self, item: T) -> Result<(), ()> {
        let mut state = match self.state.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        while state.items.len() >= self.capacity && !state.closed {
            state = match self.not_full.wait(state) {
                Ok(guard) => guard,
                Err(poisoned) => poisoned.into_inner(),
            };
        }
        if state.closed {
            return Err(());
        }
        state.items.push_back(item);
        self.not_empty.notify_one();
        Ok(())
    }

    fn pop(&self) -> Option<T> {
        let mut state = match self.state.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        while state.items.is_empty() && !state.closed {
            state = match self.not_empty.wait(state) {
                Ok(guard) => guard,
                Err(poisoned) => poisoned.into_inner(),
            };
        }
        if let Some(item) = state.items.pop_front() {
            self.not_full.notify_one();
            Some(item)
        } else {
            None
        }
    }

    fn close(&self) {
        let mut state = match self.state.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        state.closed = true;
        self.not_empty.notify_all();
        self.not_full.notify_all();
    }
}

/// An owned, numbered batch handed from the reader stage to mapper workers.
///
/// The batch number is assigned by the reader and is used only for
/// resequencing.  The reads themselves remain in source order within the
/// batch.
#[derive(Debug)]
pub struct ReadBatch<R> {
    pub batch_id: u64,
    pub reads: Vec<R>,
}

impl<R> ReadBatch<R> {
    pub fn new(batch_id: u64, reads: Vec<R>) -> Self {
        Self { batch_id, reads }
    }

    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    pub fn into_reads(self) -> Vec<R> {
        self.reads
    }
}

/// A numbered batch produced by one mapper worker.
///
/// `results` is in the same order as the reads in the corresponding
/// [`ReadBatch`].  Batches may arrive out of order at the sink, but
/// [`WorkerPool::run`] always invokes the sink in increasing `batch_id` order.
#[derive(Debug)]
pub struct MappedBatch<T> {
    pub batch_id: u64,
    pub results: Vec<T>,
    /// Bytes the worker already encoded for this batch, if an encoder was
    /// given. The sink then only has to write them, in batch order.
    pub encoded: Vec<u8>,
}

impl<T> MappedBatch<T> {
    pub fn new(batch_id: u64, results: Vec<T>) -> Self {
        Self {
            batch_id,
            results,
            encoded: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.results.len()
    }

    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }

    pub fn into_results(self) -> Vec<T> {
        self.results
    }
}

/// Counters returned after a successful pool run.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct WorkerPoolStats {
    /// Number of mapper workers started for this run.
    pub workers: usize,
    /// Maximum number of reads in a mapper batch.
    pub chunk_size: usize,
    /// Number of reads collected by the reader before splitting into chunks.
    pub reader_batch_size: usize,
    /// Number of batches delivered to the ordered sink.
    pub batches_written: usize,
    /// Number of mapped results delivered to the ordered sink.
    pub reads_written: usize,
}

/// Invalid execution settings detected before a pool is started.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WorkerPoolConfigError {
    ZeroWorkers,
    ZeroChunkSize,
    ZeroReaderBatchSize,
}

impl fmt::Display for WorkerPoolConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::ZeroWorkers => "worker pool requires at least one worker",
            Self::ZeroChunkSize => "worker pool chunk_size must be greater than zero",
            Self::ZeroReaderBatchSize => "worker pool reader_batch_size must be greater than zero",
        })
    }
}

impl std::error::Error for WorkerPoolConfigError {}

/// Errors raised by the reader, mapper, ordered sink, or pool protocol.
///
/// The three type parameters intentionally stay separate: an adapter can
/// preserve its reader, aligner, and output error types instead of converting
/// everything to a string at this boundary.
#[derive(Debug)]
pub enum WorkerPoolError<SourceError, MapperError, SinkError> {
    InvalidConfig(WorkerPoolConfigError),
    Source(SourceError),
    Mapper(MapperError),
    Sink(SinkError),
    /// A worker disappeared before all numbered batches were delivered.
    ChannelClosed,
    /// A producer/worker returned the same batch number twice.
    DuplicateBatch(u64),
    /// A reader or mapper thread panicked.
    ThreadPanicked,
}

impl<S, M, W> fmt::Display for WorkerPoolError<S, M, W>
where
    S: fmt::Debug,
    M: fmt::Debug,
    W: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig(error) => write!(f, "invalid worker-pool config: {error}"),
            Self::Source(error) => write!(f, "read source failed: {error:?}"),
            Self::Mapper(error) => write!(f, "mapper failed: {error:?}"),
            Self::Sink(error) => write!(f, "output sink failed: {error:?}"),
            Self::ChannelClosed => f.write_str("worker output channel closed before completion"),
            Self::DuplicateBatch(id) => write!(f, "duplicate mapped batch {id}"),
            Self::ThreadPanicked => f.write_str("worker-pool thread panicked"),
        }
    }
}

impl<S, M, W> std::error::Error for WorkerPoolError<S, M, W>
where
    S: fmt::Debug + 'static,
    M: fmt::Debug + 'static,
    W: fmt::Debug + 'static,
{
}

/// The one RS-LRA scheduling implementation.
///
/// `WorkerPoolConfig` is copied into the runner so a caller may reuse its
/// configuration object for more than one input.  A raw and a mapped queue
/// are both bounded to four times the worker count, with a minimum of 16,
/// matching FlashMap's current LR pipeline.
#[derive(Clone, Debug)]
pub struct WorkerPool {
    config: WorkerPoolConfig,
    queue_capacity: usize,
    reader_batch_size: usize,
}

impl WorkerPool {
    /// Construct a pool without starting any threads.
    pub fn new(config: WorkerPoolConfig) -> Result<Self, WorkerPoolConfigError> {
        if config.workers == 0 {
            return Err(WorkerPoolConfigError::ZeroWorkers);
        }
        if config.chunk_size == 0 {
            return Err(WorkerPoolConfigError::ZeroChunkSize);
        }
        if config.reader_batch_size == Some(0) {
            return Err(WorkerPoolConfigError::ZeroReaderBatchSize);
        }

        let queue_capacity = config.workers.saturating_mul(4).max(16);
        let reader_batch_size = config.reader_batch_size.unwrap_or(config.chunk_size);
        Ok(Self {
            config,
            queue_capacity,
            reader_batch_size,
        })
    }

    pub fn config(&self) -> &WorkerPoolConfig {
        &self.config
    }

    pub fn worker_count(&self) -> usize {
        self.config.workers
    }

    pub fn chunk_size(&self) -> usize {
        self.config.chunk_size
    }

    pub fn reader_batch_size(&self) -> usize {
        self.reader_batch_size
    }

    /// Map a source and deliver completed batches to an ordered sink.
    ///
    /// `source` yields owned reads (or source errors).  The source is consumed
    /// by a dedicated reader thread.  `mapper` is shared by all workers and
    /// must return an owned result; `sink` runs on the calling thread and is
    /// invoked once per batch in increasing `batch_id` order.
    pub fn run<I, R, T, SE, ME, WE, F, W>(
        &self,
        source: I,
        mapper: F,
        sink: W,
    ) -> Result<WorkerPoolStats, WorkerPoolError<SE, ME, WE>>
    where
        I: IntoIterator<Item = Result<R, SE>> + Send,
        R: Send,
        T: Send,
        SE: Send,
        ME: Send,
        F: Fn(R) -> Result<T, ME> + Sync,
        W: FnMut(MappedBatch<T>) -> Result<(), WE>,
    {
        self.run_with_encoder(source, mapper, |_, _| Ok(()), sink)
    }

    /// As [`run`], with a per-batch encoder that runs on the worker.
    ///
    /// The collector has to serialize batch order, so anything it does is
    /// single-threaded no matter how many workers there are. Formatting a
    /// whole-genome run there is hundreds of gigabytes of work on one core.
    /// `encoder` moves that onto the worker that produced the batch, leaving
    /// the collector with a write.
    pub fn run_with_encoder<I, R, T, SE, ME, WE, F, G, W>(
        &self,
        source: I,
        mapper: F,
        encoder: G,
        mut sink: W,
    ) -> Result<WorkerPoolStats, WorkerPoolError<SE, ME, WE>>
    where
        I: IntoIterator<Item = Result<R, SE>> + Send,
        R: Send,
        T: Send,
        SE: Send,
        ME: Send,
        F: Fn(R) -> Result<T, ME> + Sync,
        G: Fn(&[T], &mut Vec<u8>) -> Result<(), ME> + Sync,
        W: FnMut(MappedBatch<T>) -> Result<(), WE>,
    {
        let cancellation = Arc::new(AtomicBool::new(false));
        let raw_capacity = self.queue_capacity;
        let mapped_capacity = self.queue_capacity;
        let chunk_size = self.config.chunk_size;
        let reader_batch_size = self.reader_batch_size;
        let worker_count = self.config.workers;

        thread::scope(|scope| {
            let raw_queue = Arc::new(BoundedQueue::<Result<ReadBatch<R>, SE>>::new(raw_capacity));
            let mapped_queue = Arc::new(BoundedQueue::<
                Result<MappedBatch<T>, InternalFailure<SE, ME>>,
            >::new(mapped_capacity));
            let active_workers = Arc::new(AtomicUsize::new(worker_count));

            let reader_cancel = Arc::clone(&cancellation);
            let reader_raw = Arc::clone(&raw_queue);
            let reader_handle = scope.spawn(move || {
                reader_loop(
                    source,
                    reader_raw,
                    chunk_size,
                    reader_batch_size,
                    reader_cancel,
                );
            });

            let mut worker_handles = Vec::with_capacity(worker_count);
            for _ in 0..worker_count {
                let raw_q = Arc::clone(&raw_queue);
                let mapped_q = Arc::clone(&mapped_queue);
                let worker_cancel = Arc::clone(&cancellation);
                let workers_remaining = Arc::clone(&active_workers);
                let mapper_ref = &mapper;
                let encoder_ref = &encoder;
                worker_handles.push(scope.spawn(move || {
                    worker_loop(
                        raw_q,
                        mapped_q,
                        mapper_ref,
                        encoder_ref,
                        worker_cancel,
                        workers_remaining,
                    );
                }));
            }

            let mut stats = WorkerPoolStats {
                workers: worker_count,
                chunk_size,
                reader_batch_size,
                ..WorkerPoolStats::default()
            };
            let mut pending = BTreeMap::<u64, MappedBatch<T>>::new();
            let mut next_batch_id = 0u64;

            let mut sink_result = Ok(());
            let mut pool_error = None;
            while let Some(message) = mapped_queue.pop() {
                match message {
                    Ok(batch) => {
                        let batch_id = batch.batch_id;
                        if pending.insert(batch_id, batch).is_some() {
                            cancellation.store(true, Ordering::Release);
                            raw_queue.close();
                            mapped_queue.close();
                            pool_error = Some(WorkerPoolError::DuplicateBatch(batch_id));
                            break;
                        }

                        while let Some(batch) = pending.remove(&next_batch_id) {
                            stats.batches_written += 1;
                            stats.reads_written += batch.results.len();
                            if let Err(error) = sink(batch) {
                                cancellation.store(true, Ordering::Release);
                                raw_queue.close();
                                mapped_queue.close();
                                sink_result = Err(error);
                                break;
                            }
                            next_batch_id += 1;
                        }
                        if sink_result.is_err() {
                            break;
                        }
                    }
                    Err(error) => {
                        cancellation.store(true, Ordering::Release);
                        raw_queue.close();
                        mapped_queue.close();
                        pool_error = Some(error.into_public());
                        break;
                    }
                }
            }

            raw_queue.close();
            mapped_queue.close();

            let mut thread_panicked = false;
            for handle in worker_handles {
                thread_panicked |= handle.join().is_err();
            }
            thread_panicked |= reader_handle.join().is_err();

            if let Some(error) = pool_error {
                return Err(error);
            }
            if let Err(error) = sink_result {
                return Err(WorkerPoolError::Sink(error));
            }
            if thread_panicked {
                return Err(WorkerPoolError::ThreadPanicked);
            }
            if !pending.is_empty() {
                return Err(WorkerPoolError::ChannelClosed);
            }
            Ok(stats)
        })
    }

    /// Convenience form that collects the ordered mapped values into one
    /// vector.  Streaming callers should use [`WorkerPool::run`] to avoid
    /// retaining all results in memory.
    pub fn map<I, R, T, SE, ME, F>(
        &self,
        source: I,
        mapper: F,
    ) -> Result<Vec<T>, WorkerPoolError<SE, ME, Infallible>>
    where
        I: IntoIterator<Item = Result<R, SE>> + Send,
        R: Send,
        T: Send,
        SE: Send,
        ME: Send,
        F: Fn(R) -> Result<T, ME> + Sync,
    {
        let mut results = Vec::new();
        self.run(source, mapper, |batch| {
            results.extend(batch.results);
            Ok::<(), Infallible>(())
        })?;
        Ok(results)
    }
}

#[derive(Debug)]
enum InternalFailure<S, M> {
    Source(S),
    Mapper(M),
}

impl<S, M> InternalFailure<S, M> {
    fn into_public<W>(self) -> WorkerPoolError<S, M, W> {
        match self {
            Self::Source(error) => WorkerPoolError::Source(error),
            Self::Mapper(error) => WorkerPoolError::Mapper(error),
        }
    }
}

fn reader_loop<I, R, E>(
    source: I,
    raw_queue: Arc<BoundedQueue<Result<ReadBatch<R>, E>>>,
    chunk_size: usize,
    reader_batch_size: usize,
    cancellation: Arc<AtomicBool>,
) where
    I: IntoIterator<Item = Result<R, E>>,
{
    let mut source = source.into_iter();
    let mut batch_id = 0u64;

    loop {
        if cancellation.load(Ordering::Acquire) {
            raw_queue.close();
            return;
        }

        let mut reader_batch = Vec::with_capacity(reader_batch_size);
        for _ in 0..reader_batch_size {
            if cancellation.load(Ordering::Acquire) {
                raw_queue.close();
                return;
            }
            match source.next() {
                Some(Ok(read)) => reader_batch.push(read),
                Some(Err(error)) => {
                    if !cancellation.load(Ordering::Acquire) {
                        let _ = raw_queue.push(Err(error));
                    }
                    raw_queue.close();
                    return;
                }
                None => break,
            }
        }

        if reader_batch.is_empty() {
            raw_queue.close();
            return;
        }

        let mut reads = reader_batch.into_iter();
        loop {
            let mut chunk = Vec::with_capacity(chunk_size);
            for _ in 0..chunk_size {
                match reads.next() {
                    Some(read) => chunk.push(read),
                    None => break,
                }
            }
            if chunk.is_empty() {
                break;
            }
            if cancellation.load(Ordering::Acquire) {
                raw_queue.close();
                return;
            }
            if raw_queue.push(Ok(ReadBatch::new(batch_id, chunk))).is_err() {
                return;
            }
            batch_id += 1;
        }
    }
}

struct WorkerGuard<T, SE, ME> {
    active: Arc<AtomicUsize>,
    mapped_queue: MappedQueue<T, SE, ME>,
}

impl<T, SE, ME> Drop for WorkerGuard<T, SE, ME> {
    fn drop(&mut self) {
        if self.active.fetch_sub(1, Ordering::AcqRel) == 1 {
            self.mapped_queue.close();
        }
    }
}

fn worker_loop<R, T, SE, ME, F, G>(
    raw_queue: RawQueue<R, SE>,
    mapped_queue: MappedQueue<T, SE, ME>,
    mapper: &F,
    encoder: &G,
    cancellation: Arc<AtomicBool>,
    active_workers: Arc<AtomicUsize>,
) where
    R: Send,
    T: Send,
    SE: Send,
    ME: Send,
    F: Fn(R) -> Result<T, ME> + Sync,
    G: Fn(&[T], &mut Vec<u8>) -> Result<(), ME> + Sync,
{
    let _guard = WorkerGuard {
        active: active_workers,
        mapped_queue: Arc::clone(&mapped_queue),
    };

    while let Some(message) = raw_queue.pop() {
        if cancellation.load(Ordering::Acquire) {
            return;
        }

        let raw_batch = match message {
            Ok(batch) => batch,
            Err(error) => {
                cancellation.store(true, Ordering::Release);
                raw_queue.close();
                let _ = mapped_queue.push(Err(InternalFailure::Source(error)));
                return;
            }
        };

        let batch_id = raw_batch.batch_id;
        let mut results = Vec::with_capacity(raw_batch.reads.len());
        for read in raw_batch.reads {
            if cancellation.load(Ordering::Acquire) {
                return;
            }
            match mapper(read) {
                Ok(result) => results.push(result),
                Err(error) => {
                    if !cancellation.swap(true, Ordering::AcqRel) {
                        raw_queue.close();
                        let _ = mapped_queue.push(Err(InternalFailure::Mapper(error)));
                    }
                    return;
                }
            }
        }

        if cancellation.load(Ordering::Acquire) {
            return;
        }
        // Encode here, on the worker that mapped the batch, so the collector
        // is left with an ordered write rather than the formatting too.
        let mut encoded = Vec::new();
        if let Err(error) = encoder(&results, &mut encoded) {
            if !cancellation.swap(true, Ordering::AcqRel) {
                raw_queue.close();
                let _ = mapped_queue.push(Err(InternalFailure::Mapper(error)));
            }
            return;
        }

        if mapped_queue
            .push(Ok(MappedBatch {
                batch_id,
                results,
                encoded,
            }))
            .is_err()
        {
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::sync::{Arc, Barrier};
    use std::thread;
    use std::time::Duration;

    fn config(
        workers: usize,
        chunk_size: usize,
        reader_batch_size: Option<usize>,
    ) -> WorkerPoolConfig {
        WorkerPoolConfig {
            workers,
            chunk_size,
            reader_batch_size,
        }
    }

    #[test]
    fn ordered_sink_resequences_out_of_order_workers() {
        let pool = WorkerPool::new(config(3, 1, None)).unwrap();
        let source = (0..12).map(Ok::<_, Infallible>);
        let mut observed = Vec::new();

        let stats = pool
            .run(
                source,
                |value| {
                    // Force later batches to finish first without changing
                    // the ordering contract exposed to the sink.
                    thread::sleep(Duration::from_millis((12 - value) as u64));
                    Ok::<_, Infallible>(value * 10)
                },
                |batch| {
                    observed.extend(batch.results);
                    Ok::<_, Infallible>(())
                },
            )
            .unwrap();

        assert_eq!(
            observed,
            (0..12).map(|value| value * 10).collect::<Vec<_>>()
        );
        assert_eq!(stats.batches_written, 12);
        assert_eq!(stats.reads_written, 12);
    }

    #[test]
    fn reader_splits_batches_at_configured_chunk_size() {
        let pool = WorkerPool::new(config(2, 3, Some(5))).unwrap();
        let source = (0..13).map(Ok::<_, Infallible>);
        let mut batches = Vec::new();

        let stats = pool
            .run(source, Ok::<_, Infallible>, |batch| {
                batches.push((batch.batch_id, batch.results));
                Ok::<_, Infallible>(())
            })
            .unwrap();

        assert_eq!(
            batches
                .iter()
                .map(|(batch_id, _)| *batch_id)
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 4]
        );
        assert_eq!(
            batches
                .iter()
                .map(|(_, values)| values.len())
                .collect::<Vec<_>>(),
            vec![3, 2, 3, 2, 3]
        );
        assert_eq!(stats.chunk_size, 3);
        assert_eq!(stats.reader_batch_size, 5);
        assert_eq!(stats.batches_written, 5);
    }

    #[test]
    fn configured_worker_count_is_used() {
        let pool = WorkerPool::new(config(3, 1, None)).unwrap();
        assert_eq!(pool.worker_count(), 3);

        // The barrier makes the first three mapper calls overlap, proving
        // that all configured workers participate rather than merely
        // checking a copied configuration value.
        let barrier = Arc::new(Barrier::new(3));
        let seen = Arc::new(Mutex::new(HashSet::new()));
        let barrier_for_mapper = Arc::clone(&barrier);
        let seen_for_mapper = Arc::clone(&seen);
        let source = (0..3).map(Ok::<_, Infallible>);
        pool.run(
            source,
            move |value| {
                barrier_for_mapper.wait();
                seen_for_mapper
                    .lock()
                    .unwrap()
                    .insert(thread::current().id());
                Ok::<_, Infallible>(value)
            },
            |_| Ok::<_, Infallible>(()),
        )
        .unwrap();

        assert_eq!(seen.lock().unwrap().len(), 3);
    }

    #[test]
    fn source_and_mapper_errors_are_forwarded() {
        let pool = WorkerPool::new(config(2, 1, None)).unwrap();
        let source = vec![Ok(1), Ok(2), Err("source")].into_iter();
        let error = pool.map(source, Ok::<_, &str>).unwrap_err();
        assert!(matches!(error, WorkerPoolError::Source("source")));

        let source = (0..4).map(Ok::<_, Infallible>);
        let error = pool
            .map(
                source,
                |value| {
                    if value == 2 {
                        Err("mapper")
                    } else {
                        Ok(value)
                    }
                },
            )
            .unwrap_err();
        assert!(matches!(error, WorkerPoolError::Mapper("mapper")));
    }
}
