//! HiGHS source: util/HSet.h, util/HSet.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use crate::prelude::*;

/// HiGHS: HSet.h::HSet
#[derive(Clone, Debug)]
pub struct HSet {
    pub count_: HighsInt,
    pub entry_: Vec<HighsInt>,
    pub setup_: bool,
    pub debug_: bool,
    pub allow_assert_: bool,
    pub output_flag_: bool,
    pub log_file_: Option<HighsLogStream>,
    pub max_entry_: HighsInt,
    pub pointer_: Vec<HighsInt>,
}

impl Default for HSet {
    fn default() -> Self {
        Self {
            count_: 0,
            entry_: Default::default(),
            setup_: false,
            debug_: false,
            allow_assert_: true,
            output_flag_: false,
            log_file_: None,
            max_entry_: 0,
            pointer_: Default::default(),
        }
    }
}

impl HSet {
    /// HiGHS: HSet.cpp::HSet::setup (lines 15-31, 17 lines)
    /// C++ default arguments: output_flag = false; log_file = NULL; debug = false; allow_assert = true
    pub fn setup(
        &mut self,
        size: HighsInt,
        max_entry: HighsInt,
        output_flag: bool,
        log_file: Option<HighsLogStream>,
        debug: bool,
        allow_assert: bool,
    ) -> bool {
        // HiGHS: HSet.h file-scope constants (C++ `min_entry`, `no_pointer`).
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        self.setup_ = false;
        if size <= 0 {
            return false;
        }
        if max_entry < min_entry {
            return false;
        }
        self.max_entry_ = max_entry;
        self.debug_ = debug;
        self.allow_assert_ = allow_assert;
        self.output_flag_ = output_flag;
        self.log_file_ = log_file;
        self.entry_.resize(size as usize, 0);
        self.pointer_ = vec![no_pointer; (self.max_entry_ + 1) as usize];
        self.count_ = 0;
        self.setup_ = true;
        true
    }

    /// HiGHS: HSet.cpp::HSet::clear (lines 33-38, 6 lines)
    pub fn clear(&mut self) {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        if !self.setup_ {
            self.setup(1, 0, false, None, false, true);
        }
        self.pointer_ = vec![no_pointer; (self.max_entry_ + 1) as usize];
        self.count_ = 0;
        if self.debug_ {
            self.debug();
        }
    }

    /// HiGHS: HSet.cpp::HSet::add (lines 40-64, 25 lines)
    pub fn add(&mut self, entry: HighsInt) -> bool {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        if entry < min_entry {
            return false;
        }
        if !self.setup_ {
            self.setup(1, entry, false, None, false, true);
        }
        if entry > self.max_entry_ {
            // Entry exceeds what's allowable so far so can't be in the list
            self.pointer_.resize((entry + 1) as usize, 0);
            for ix in (self.max_entry_ + 1)..entry {
                self.pointer_[ix as usize] = no_pointer;
            }
            self.max_entry_ = entry;
        } else if self.pointer_[entry as usize] > no_pointer {
            // Duplicate
            if self.debug_ {
                self.debug();
            }
            return false;
        }
        // New entry
        let mut size = self.entry_.len() as HighsInt;
        if self.count_ == size {
            size += 1;
            self.entry_.resize(size as usize, 0);
        }
        self.pointer_[entry as usize] = self.count_;
        self.entry_[self.count_ as usize] = entry;
        self.count_ += 1;
        if self.debug_ {
            self.debug();
        }
        true
    }

    /// HiGHS: HSet.cpp::HSet::remove (lines 66-85, 20 lines)
    pub fn remove(&mut self, entry: HighsInt) -> bool {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        if !self.setup_ {
            self.setup(1, 0, false, None, false, true);
            if self.debug_ {
                self.debug();
            }
            return false;
        }
        if entry < min_entry {
            return false;
        }
        if entry > self.max_entry_ {
            return false;
        }
        let pointer = self.pointer_[entry as usize];
        if pointer == no_pointer {
            return false;
        }
        self.pointer_[entry as usize] = no_pointer;
        if pointer < self.count_ - 1 {
            let last_entry = self.entry_[(self.count_ - 1) as usize];
            self.entry_[pointer as usize] = last_entry;
            self.pointer_[last_entry as usize] = pointer;
        }
        self.count_ -= 1;
        if self.debug_ {
            self.debug();
        }
        true
    }

    /// HiGHS: HSet.cpp::HSet::in (lines 87-91, 5 lines)
    pub fn r#in(&self, entry: HighsInt) -> bool {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        if entry < min_entry {
            return false;
        }
        if entry > self.max_entry_ {
            return false;
        }
        self.pointer_[entry as usize] != no_pointer
    }

    /// HiGHS: HSet.cpp::HSet::debug (lines 93-169, 77 lines)
    pub fn debug(&self) -> bool {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        // C++ fprintf(log_file_, ...) is only reached when output_flag_ is set;
        // a null log_file_ would be UB in C++, so it is treated as no output here.
        let log_on = self.output_flag_ && self.log_file_.is_some();
        if !self.setup_ {
            if log_on {
                println!("HSet: ERROR setup_ not called");
            }
            if self.allow_assert_ {
                debug_assert!(self.setup_);
            }
            return false;
        }
        let max_entry_ok = self.max_entry_ >= min_entry;
        if !max_entry_ok {
            if log_on {
                println!(
                    "HSet: ERROR max_entry_ = {} < {}",
                    self.max_entry_, min_entry
                );
                self.print();
            }
            if self.allow_assert_ {
                debug_assert!(max_entry_ok);
            }
            return false;
        }
        let size = self.entry_.len() as HighsInt;
        let size_count_ok = size >= self.count_;
        if !size_count_ok {
            if log_on {
                println!(
                    "HSet: ERROR entry_.size() = {} is less than count_ = {}",
                    size, self.count_
                );
                self.print();
            }
            if self.allow_assert_ {
                debug_assert!(size_count_ok);
            }
            return false;
        }
        // Check pointer_ is consistent with count_ and entry_
        let mut count: HighsInt = 0;
        for ix in 0..=self.max_entry_ {
            let pointer = self.pointer_[ix as usize];
            if pointer == no_pointer {
                continue;
            }
            let pointer_ok = pointer >= 0 && pointer < self.count_;
            if !pointer_ok {
                if log_on {
                    println!(
                        "HSet: ERROR pointer_[{}] = {} is not in [0, {}]",
                        ix, pointer, self.count_
                    );
                    self.print();
                }
                if self.allow_assert_ {
                    debug_assert!(pointer_ok);
                }
                return false;
            }
            count += 1;
            let entry = self.entry_[pointer as usize];
            let entry_ok = entry == ix;
            if !entry_ok {
                if log_on {
                    println!("HSet: ERROR entry_[{}] is {}, not {}", pointer, entry, ix);
                    self.print();
                }
                if self.allow_assert_ {
                    debug_assert!(entry_ok);
                }
                return false;
            }
        }
        let count_ok = count == self.count_;
        if !count_ok {
            if log_on {
                println!(
                    "HSet: ERROR pointer_ has {} pointers, not {}",
                    count, self.count_
                );
                self.print();
            }
            if self.allow_assert_ {
                debug_assert!(count_ok);
            }
            return false;
        }
        true
    }

    /// HiGHS: HSet.cpp::HSet::print (lines 171-197, 27 lines)
    pub fn print(&self) {
        const min_entry: HighsInt = 0;
        const no_pointer: HighsInt = min_entry - 1;
        if !self.setup_ {
            return;
        }
        if self.log_file_.is_none() {
            return;
        }
        let size = self.entry_.len() as HighsInt;
        println!("\nSet({}, {}):", size, self.max_entry_);
        print!("Pointers: Pointers|");
        for ix in 0..=self.max_entry_ {
            if self.pointer_[ix as usize] != no_pointer {
                print!(" {:>4}", self.pointer_[ix as usize]);
            }
        }
        println!();
        print!("          Entries |");
        for ix in 0..=self.max_entry_ {
            if self.pointer_[ix as usize] != no_pointer {
                print!(" {:>4}", ix);
            }
        }
        println!();
        print!("Entries:  Indices |");
        for ix in 0..self.count_ {
            print!(" {:>4}", ix);
        }
        println!();
        print!("          Entries |");
        for ix in 0..self.count_ {
            print!(" {:>4}", self.entry_[ix as usize]);
        }
        println!();
    }

    /// HiGHS: HSet.h::HSet::count (inline getter)
    pub fn count(&self) -> &HighsInt {
        &self.count_
    }

    /// HiGHS: HSet.h::HSet::entry (inline getter)
    pub fn entry(&self) -> &[HighsInt] {
        &self.entry_
    }
}
