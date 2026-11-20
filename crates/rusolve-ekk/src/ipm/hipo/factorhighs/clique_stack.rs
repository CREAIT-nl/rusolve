//! HiGHS source: ipm/hipo/factorhighs/CliqueStack.h, ipm/hipo/factorhighs/CliqueStack.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The serial stack of contribution blocks (cliques). The C++ hands out raw
//! `double*` into `stack_`; the Rust hands out `Int64` offsets into `stack_`,
//! which stay valid across a reallocation (the
//! C++ `reallocation` flag is still reported, callers use it exactly as C++).
//! `std::stack<std::pair<Int, Int64>>` is a `Vec` used as a stack.

use super::super::{Int, Int64};

/// HiGHS: CliqueStack.h `class CliqueStack`
#[derive(Clone, Debug)]
pub struct CliqueStack {
    pub stack_: Vec<f64>,
    /// C++ `double* workspace_` (uninitialised until `setup`): offset into `stack_`.
    pub workspace_: Option<Int64>,
    /// No initialiser in C++; `setup` sets it.
    pub worksize_: Int64,
    pub empty_: bool,
    pub sn_pushed_: Vec<(Int, Int64)>,
    pub top_: Int64,
}

impl Default for CliqueStack {
    fn default() -> Self {
        Self {
            stack_: Vec::new(),
            workspace_: None,
            worksize_: 0,
            empty_: true,
            sn_pushed_: Vec::new(),
            top_: 0,
        }
    }
}

impl CliqueStack {
    /// HiGHS: CliqueStack.cpp::CliqueStack::init (lines 8-15)
    pub fn init(&mut self, stack_size: Int64) {
        self.top_ = 0;
        self.workspace_ = None;
        self.worksize_ = 0;
        self.empty_ = false;

        self.stack_.resize(stack_size as usize, 0.0);
    }

    /// HiGHS: CliqueStack.cpp::CliqueStack::setup (lines 17-41)
    /// Returns the workspace offset (`None` = the C++ `nullptr` for an empty clique).
    pub fn setup(&mut self, clique_size: Int64, reallocation: &mut bool) -> Option<Int64> {
        // Clear workspace

        debug_assert!(self.workspace_.is_none() && self.worksize_ == 0);
        *reallocation = false;

        if clique_size > 0 {
            // This should not trigger reallocation, because the resize in init is done
            // with the maximum possible size of the stack.
            if self.top_ + clique_size > self.stack_.len() as Int64 {
                *reallocation = true;
                self.stack_.resize((self.top_ + clique_size) as usize, 0.0);
            }

            // accessing stack[top] is valid only if the clique is not empty,
            // otherwise it may be out of bounds.
            self.workspace_ = Some(self.top_);
            self.worksize_ = clique_size;

            // initialize workspace to zero
            let ws = self.workspace_.unwrap() as usize;
            let ws_len = self.worksize_ as usize;
            self.stack_[ws..ws + ws_len].fill(0.0);
        }

        self.workspace_
    }

    /// HiGHS: CliqueStack.cpp::CliqueStack::empty (line 43; lines 43-43)
    pub fn empty(&self) -> bool {
        self.empty_
    }

    /// HiGHS: CliqueStack.cpp::CliqueStack::getChild (lines 45-54)
    /// Returns the offset of the child clique in `stack_`.
    pub fn getChild(&self, child_sn: &mut Int) -> Int64 {
        // Get the top of the stack, in terms of supernode ID of the child and pointer
        // to its data.

        *child_sn = self.sn_pushed_.last().unwrap().0;
        let child_size = self.sn_pushed_.last().unwrap().1;
        self.top_ - child_size
    }

    /// HiGHS: CliqueStack.cpp::CliqueStack::popChild (lines 56-63)
    pub fn popChild(&mut self) {
        // Remove top child from the stack

        let child_size = self.sn_pushed_.last().unwrap().1;
        self.sn_pushed_.pop();

        self.top_ -= child_size;
    }

    /// HiGHS: CliqueStack.cpp::CliqueStack::pushWork (lines 65-80)
    /// `memmove` -> `copy_within`.
    pub fn pushWork(&mut self, sn: Int) {
        // Put the content of the workspace at the top of the stack

        if self.worksize_ > 0 {
            // stack_[top_] has lower address than workspace, so no need to resize.
            let ws = self.workspace_.unwrap() as usize;
            let ws_len = self.worksize_ as usize;
            let top = self.top_ as usize;
            self.stack_.copy_within(ws..ws + ws_len, top);

            self.top_ += self.worksize_;

            // keep track of supernodes pushed
            self.sn_pushed_.push((sn, self.worksize_));
        }

        self.worksize_ = 0;
        self.workspace_ = None;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn push_and_pop() {
        let mut cs = CliqueStack::default();
        cs.init(16);
        let mut realloc = false;
        let off = cs.setup(4, &mut realloc).unwrap();
        cs.stack_[off as usize] = 1.5;
        cs.pushWork(3);
        assert!(!cs.empty());
        let mut child = -1;
        let c = cs.getChild(&mut child);
        assert_eq!(child, 3);
        assert_eq!(cs.stack_[c as usize], 1.5);
        cs.popChild();
        // C++ CliqueStack.cpp lines 12 and 43: `empty_` is set false by init and
        // never set back, so `empty()` stays false after the pop.
        assert!(!cs.empty());
    }
}
