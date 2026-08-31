//! Fast, zero-dependency hasher for non-cryptographic internal maps and sets.
//!
//! Uses the FxHash algorithm (multiplication + rotation) identical to rustc
//! and FlashMap. This eliminates SipHash overhead on small integer and k-mer keys.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
pub struct FxHasher {
    hash: usize,
}

impl Hasher for FxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.hash = self.hash.rotate_left(5) ^ (byte as usize);
            self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.hash = self.hash.rotate_left(5) ^ (i as usize);
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.hash = self.hash.rotate_left(5) ^ (i as usize);
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.hash = self.hash.rotate_left(5) ^ (i as usize);
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.hash = self.hash.rotate_left(5) ^ i;
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn write_i32(&mut self, i: i32) {
        self.hash = self.hash.rotate_left(5) ^ (i as usize);
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn write_i64(&mut self, i: i64) {
        self.hash = self.hash.rotate_left(5) ^ (i as usize);
        self.hash = self.hash.wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

pub type FxBuildHasher = BuildHasherDefault<FxHasher>;
pub type FxHashMap<K, V> = HashMap<K, V, FxBuildHasher>;
pub type FxHashSet<T> = HashSet<T, FxBuildHasher>;

pub trait FxHashMapExt {
    fn new() -> Self;
    #[allow(dead_code)]
    fn with_capacity(capacity: usize) -> Self;
}

impl<K, V> FxHashMapExt for FxHashMap<K, V> {
    #[inline]
    fn new() -> Self {
        HashMap::with_hasher(FxBuildHasher::default())
    }

    #[inline]
    fn with_capacity(capacity: usize) -> Self {
        HashMap::with_capacity_and_hasher(capacity, FxBuildHasher::default())
    }
}

pub trait FxHashSetExt {
    fn new() -> Self;
    #[allow(dead_code)]
    fn with_capacity(capacity: usize) -> Self;
}

impl<T> FxHashSetExt for FxHashSet<T> {
    #[inline]
    fn new() -> Self {
        HashSet::with_hasher(FxBuildHasher::default())
    }

    #[inline]
    fn with_capacity(capacity: usize) -> Self {
        HashSet::with_capacity_and_hasher(capacity, FxBuildHasher::default())
    }
}
