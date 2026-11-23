//! HiGHS source: parallel/HighsCacheAlign.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete. `highs::cache_aligned` over-allocates and stores the original pointer in front
//! of the aligned block (lines 21-42). Rust gets the same 64-byte placement from the type:
//! `CacheAligned<T>` is `#[repr(C, align(64))]`, so `Box<CacheAligned<T>>` is
//! `cache_aligned::unique_ptr<T>` and `Arc<CacheAligned<T>>` is `cache_aligned::shared_ptr<T>`.
//! Alignment only separates cache lines between threads; no result depends on it.

use std::ops::{Deref, DerefMut};

/// HiGHS: HighsCacheAlign.h::cache_aligned::alignment (line 19)
pub const kCacheAlignment: usize = 64;

/// HiGHS: HighsCacheAlign.h::cache_aligned (lines 18-80): a value on its own cache line(s).
#[repr(C, align(64))]
#[derive(Debug, Default)]
pub struct CacheAligned<T>(pub T);

impl<T> CacheAligned<T> {
    /// HiGHS: HighsCacheAlign.h::cache_aligned::make_unique (lines 69-72)
    pub fn new(value: T) -> Self {
        CacheAligned(value)
    }
}

impl<T> Deref for CacheAligned<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T> DerefMut for CacheAligned<T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

/// HiGHS: HighsCacheAlign.h::cache_aligned::unique_ptr<T> (lines 57-58)
pub type CacheAlignedBox<T> = Box<CacheAligned<T>>;

#[cfg(test)]
mod tests {
    // Unit tests.
    use super::*;

    #[test]
    fn cache_aligned_values_start_on_a_64_byte_boundary() {
        let v: Vec<CacheAlignedBox<u8>> = (0..4).map(|i| Box::new(CacheAligned::new(i))).collect();
        for b in &v {
            assert_eq!(
                (&**b as *const CacheAligned<u8> as usize) % kCacheAlignment,
                0
            );
        }
        assert_eq!(std::mem::align_of::<CacheAligned<u8>>(), kCacheAlignment);
        assert_eq!(**v[3], 3);
    }
}
