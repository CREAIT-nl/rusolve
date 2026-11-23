//! Ours (no C++ counterpart): `SharedMut<T>`, the Rust spelling of "a plain C++ pointer to
//! shared mutable state that several MIP workers touch under HiGHS' own protocol".
//!
//! HiGHS MIP workers (`HighsMipWorker`, `HighsMipSolver::runTask`) hand each task
//! references to objects that are shared across threads and protected by HiGHS'
//! conventions rather than by the type system: `parallelLockActive()` /
//! `HighsSpinMutex` sections, disjoint per-worker slots, or "only the main worker after
//! `taskWait`". Wrapping such an object in `SharedMut` makes the sharing explicit at the
//! call site; every access is `unsafe` and each call site carries a `// SAFETY:` comment
//! naming the C++ protocol it follows.
//!
//! Complete; the parallel port only uses it.

use std::cell::UnsafeCell;

/// A `T` reachable mutably from several threads under a protocol the caller upholds.
#[derive(Debug, Default)]
pub struct SharedMut<T> {
    value: UnsafeCell<T>,
}

// SAFETY: `SharedMut` hands out access only through `unsafe` methods whose callers
// guarantee the C++ exclusion protocol; moving `T` between threads requires `T: Send`.
unsafe impl<T: Send> Send for SharedMut<T> {}
// SAFETY: as above; `&SharedMut<T>` is shared only under the caller's protocol.
unsafe impl<T: Send> Sync for SharedMut<T> {}

impl<T> SharedMut<T> {
    pub const fn new(value: T) -> Self {
        Self {
            value: UnsafeCell::new(value),
        }
    }

    /// Shared access.
    ///
    /// # Safety
    /// No thread holds a `&mut T` from [`SharedMut::get_mut_unchecked`] for the life of the
    /// returned reference.
    pub unsafe fn get(&self) -> &T {
        // SAFETY: forwarded to the caller
        unsafe { &*self.value.get() }
    }

    /// Exclusive access through a shared handle.
    ///
    /// # Safety
    /// The caller holds the C++ lock (or slot ownership) that makes this the only access
    /// to the `T` for the life of the returned reference.
    #[allow(clippy::mut_from_ref)]
    pub unsafe fn get_mut_unchecked(&self) -> &mut T {
        // SAFETY: forwarded to the caller
        unsafe { &mut *self.value.get() }
    }

    /// Safe exclusive access when the handle itself is unique.
    pub fn get_mut(&mut self) -> &mut T {
        self.value.get_mut()
    }

    pub fn into_inner(self) -> T {
        self.value.into_inner()
    }

    /// Raw pointer to the value (C++ `T*`).
    pub fn as_ptr(&self) -> *mut T {
        self.value.get()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disjoint_writes_from_two_threads() {
        let s = SharedMut::new([0i32; 2]);
        std::thread::scope(|sc| {
            for i in 0..2 {
                let s = &s;
                sc.spawn(move || {
                    // SAFETY: each thread writes only its own slot (disjoint)
                    let p = s.as_ptr() as *mut i32;
                    unsafe { *p.add(i) = i as i32 + 1 };
                });
            }
        });
        assert_eq!(s.into_inner(), [1, 2]);
    }

    #[test]
    fn get_mut_and_unchecked() {
        let mut s = SharedMut::new(1);
        *s.get_mut() += 1;
        // SAFETY: single thread, no other reference
        unsafe { *s.get_mut_unchecked() += 1 };
        assert_eq!(unsafe { *s.get() }, 3);
    }
}
