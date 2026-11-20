//! HiGHS source: ipm/hipo/factorhighs/DataCollector.h, ipm/hipo/factorhighs/DataCollector.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Factorisation statistics. With `HIPO_TIMING_LEVEL == 0` and
//! `HIPO_COLLECT_EXPENSIVE_DATA` undefined (the HiGHS build) every collecting
//! method is an empty body and `IterData` has no fields, so this file is
//! complete: the methods take `&self` (no mutex) and do nothing.

use super::super::Int;
use super::super::auxiliary::logger::Logger;
use super::*;

/// HiGHS: DataCollector.h `struct IterData` (empty without
/// `HIPO_COLLECT_EXPENSIVE_DATA`).
#[derive(Clone, Copy, Debug, Default)]
pub struct IterData {}

/// HiGHS: DataCollector.h `class DataCollector` (`mutex_` not ported).
#[derive(Clone, Debug)]
pub struct DataCollector {
    pub times_: Vec<f64>,
    pub blas_calls_: Vec<Int>,
    pub iter_data_record_: Vec<IterData>,
}

impl Default for DataCollector {
    fn default() -> Self {
        Self::new()
    }
}

impl DataCollector {
    /// HiGHS: DataCollector.cpp::DataCollector::DataCollector (lines 10-13)
    pub fn new() -> Self {
        Self {
            times_: vec![0.0; kTimeSize as usize],
            blas_calls_: vec![0; (kTimeBlasEnd - kTimeBlasStart + 1) as usize],
            iter_data_record_: Vec::new(),
        }
    }

    /// HiGHS: DataCollector.cpp::DataCollector::append (lines 14-17)
    pub fn append(&mut self) {
        self.iter_data_record_.push(IterData::default());
    }

    /// HiGHS: DataCollector.cpp::DataCollector::back (lines 18-22)
    pub fn back(&mut self) -> &mut IterData {
        if self.iter_data_record_.is_empty() {
            self.append();
        }
        self.iter_data_record_.last_mut().unwrap()
    }

    /// HiGHS: DataCollector.cpp::DataCollector::sumTime (lines 28-38; empty at timing level 0)
    pub fn sumTime(&self, i: Int, t: f64) {}

    /// HiGHS: DataCollector.cpp::DataCollector::setExtremeEntries (lines 39-49; empty)
    pub fn setExtremeEntries(&self, minD: f64, maxD: f64, minoffD: f64, maxoffD: f64) {}

    /// HiGHS: DataCollector.cpp::DataCollector::countRegPiv (lines 50-56; empty)
    pub fn countRegPiv(&self) {}

    /// HiGHS: DataCollector.cpp::DataCollector::countSwap (lines 57-62; empty)
    pub fn countSwap(&self) {}

    /// HiGHS: DataCollector.cpp::DataCollector::count2x2 (lines 63-68; empty)
    pub fn count2x2(&self) {}

    /// HiGHS: DataCollector.cpp::DataCollector::setWrongSign (lines 69-75; empty)
    pub fn setWrongSign(&self, p: f64) {}

    /// HiGHS: DataCollector.cpp::DataCollector::setMaxReg (lines 76-82; empty)
    pub fn setMaxReg(&self, new_reg: f64) {}

    /// HiGHS: DataCollector.cpp::DataCollector::setNorms (lines 84-90; empty)
    pub fn setNorms(&self, norm1: f64, maxdiag: f64) {}

    /// HiGHS: DataCollector.cpp::DataCollector::printTimes (lines 96-294; empty at timing level 0)
    pub fn printTimes(&self, logger: &Logger) {}

    /// HiGHS: DataCollector.cpp::DataCollector::printIter (lines 296-324; empty)
    pub fn printIter(&self, logger: &Logger) {}
}
