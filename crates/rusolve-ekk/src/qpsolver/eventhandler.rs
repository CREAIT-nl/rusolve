//! HiGHS source: qpsolver/eventhandler.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ `Eventhandler<T&>` stores `std::function<void(T&)>`; here a
//! subscriber is `Rc<dyn Fn(&mut T)>` so `Settings` stays `Clone` (C++ copies
//! `Settings` by value in `Runtime`).

use crate::prelude::*;

/// HiGHS: eventhandler.hpp::Eventhandler (lines 14-28)
pub struct Eventhandler<T> {
    subscribers: Vec<Rc<dyn Fn(&mut T)>>,
}

impl<T> Clone for Eventhandler<T> {
    fn clone(&self) -> Self {
        Self {
            subscribers: self.subscribers.clone(),
        }
    }
}

impl<T> Default for Eventhandler<T> {
    fn default() -> Self {
        Self {
            subscribers: Vec::new(),
        }
    }
}

impl<T> std::fmt::Debug for Eventhandler<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Eventhandler")
            .field("subscribers", &self.subscribers.len())
            .finish()
    }
}

impl<T> Eventhandler<T> {
    /// HiGHS: eventhandler.hpp::Eventhandler::subscribe (lines 18-20)
    pub fn subscribe(&mut self, subscriber: Rc<dyn Fn(&mut T)>) {
        self.subscribers.push(subscriber);
    }

    /// HiGHS: eventhandler.hpp::Eventhandler::fire (lines 22-26)
    pub fn fire(&self, args: &mut T) {
        for fun in &self.subscribers {
            fun(args);
        }
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    #[test]
    fn fires_every_subscriber() {
        let mut e: Eventhandler<HighsInt> = Eventhandler::default();
        let seen = Rc::new(std::cell::Cell::new(0));
        let s1 = seen.clone();
        e.subscribe(Rc::new(move |x: &mut HighsInt| s1.set(s1.get() + *x)));
        let s2 = seen.clone();
        e.subscribe(Rc::new(move |x: &mut HighsInt| s2.set(s2.get() + 10 * *x)));
        e.fire(&mut 2);
        assert_eq!(seen.get(), 22);
    }
}
