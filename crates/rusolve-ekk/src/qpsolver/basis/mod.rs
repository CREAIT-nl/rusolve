//! HiGHS source: qpsolver/basis.hpp, qpsolver/basis.cpp (lines 13-220; see `part2.rs` for 224-425)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//!
//!
//! PORT NOTES (ownership, lib.rs rule (a)):
//! * C++ `Basis` stores `Runtime& runtime`, read only in the constructor
//!   (`runtime.instance.num_var/num_con`, `rt.instance.A.t()`). It is a
//!   constructor parameter here and not stored. `new` takes `&mut Runtime`
//!   because `Matrix::t()` lazily builds and caches the transpose.
//! * `activate` / `updatebasis` take `rt: &Runtime` in place of
//!   `const Settings& settings` (every C++ caller passes `runtime.settings`;
//!   read `rt.settings.reinvertfrequency`), because the pricing callback
//!   (`Pricing::update_weights`) needs `rt` too.
//! * `Pricing* pricing` -> `pricing: &mut dyn Pricing`. In `updatebasis`
//!   compute `hvec2vec(col_aq)` / `hvec2vec(row_ep)` into locals first, then
//!   call `pricing.update_weights(rt, self, &aq, &ep, droppedcon, newactivecon)`.
//! * `HFactor::setup(..., baseindex.data())` -> `setup_raw` (the basic
//!   index is not stored; pass `&mut self.baseindex` to `build` / `update`,
//!   lib.rs rule (c)). `basisfactor.build()` -> `build(&mut baseindex, None)`.
//! * `std::map<int, BasisStatus> basisstatus` -> `BTreeMap`; C++
//!   `basisstatus[conid]` INSERTS a default (`kInactive`) when missing, so
//!   `getstatus` is `*self.basisstatus.entry(conid).or_insert(kInactive)`.

pub mod part2;

use crate::prelude::*;
use crate::qpsolver::*;
use crate::util::hfactor_const::*;

/// HiGHS: basis.hpp::Basis (lines 24-157)
#[derive(Clone, Debug)]
pub struct Basis {
    pub(crate) buffer_vec2hvec: HVector,
    pub(crate) Ztprod_res: QpVector,
    pub(crate) buffer_Zprod: QpVector,

    pub(crate) basisfactor: HFactor,
    pub(crate) updatessinceinvert: HighsInt,

    pub(crate) Atran: MatrixBase,

    /// indices of active constraints in basis
    pub(crate) active_constraint_index: Vec<HighsInt>,
    /// ids of constraints that are in the basis but not active
    pub(crate) non_active_constraint_index: Vec<HighsInt>,
    /// ids of constraints that are in the basis
    pub(crate) baseindex: Vec<HighsInt>,

    pub(crate) basisstatus: BTreeMap<i32, BasisStatus>,

    /// -1 if constraint not in basis, [0, num_var] if in basis
    pub(crate) constraintindexinbasisfactor: Vec<HighsInt>,

    pub(crate) buffer_column_aq: QpVector,
    pub(crate) buffer_row_ep: QpVector,

    pub(crate) buffered_q: HighsInt,
    pub(crate) buffered_p: HighsInt,
    pub(crate) row_ep: HVector,
    pub(crate) col_aq: HVector,

    pub(crate) reinversion_hint: bool,
    /// rusolve HiGHS 1.16 #3354 (h16_hfactor_pivot): `rt.settings.h16_hfactor_pivot`, read in
    /// the constructor and applied to every `basisfactor` that `build` creates.
    pub(crate) h16_hfactor_pivot: bool,
}

impl Basis {
    /// HiGHS: basis.hpp::Basis::vec2hvec (lines 29-38)
    pub(crate) fn vec2hvec(&mut self, vec: &QpVector) -> &mut HVector {
        self.buffer_vec2hvec.clear();
        for i in 0..vec.num_nz {
            self.buffer_vec2hvec.index[i as usize] = vec.index[i as usize];
            self.buffer_vec2hvec.array[vec.index[i as usize] as usize] =
                vec.value[vec.index[i as usize] as usize];
        }
        self.buffer_vec2hvec.count = vec.num_nz;
        self.buffer_vec2hvec.pack_flag = true;
        &mut self.buffer_vec2hvec
    }

    /// HiGHS: basis.hpp::Basis::hvec2vec(const HVector&, QpVector&) (lines 40-52)
    pub(crate) fn hvec2vec_into<'a>(hvec: &HVector, target: &'a mut QpVector) -> &'a mut QpVector {
        target.reset();
        for i in 0..hvec.count {
            target.index[i as usize] = hvec.index[i as usize];
            target.value[target.index[i as usize] as usize] =
                hvec.array[hvec.index[i as usize] as usize];
        }
        // C++ lines 46-49: commented-out loop, dropped.
        target.num_nz = hvec.count;
        target
    }

    /// HiGHS: basis.hpp::Basis::hvec2vec(const HVector&) (lines 54-58)
    pub(crate) fn hvec2vec(hvec: &HVector) -> QpVector {
        let mut vec = QpVector::new(hvec.size);
        Self::hvec2vec_into(hvec, &mut vec);
        vec
    }

    /// HiGHS: basis.hpp::Basis::getreinversionhint (line 100)
    pub fn getreinversionhint(&self) -> bool {
        self.reinversion_hint
    }

    /// HiGHS: basis.hpp::Basis::getnupdatessinceinvert (line 104)
    pub fn getnupdatessinceinvert(&self) -> HighsInt {
        self.updatessinceinvert
    }

    /// HiGHS: basis.hpp::Basis::getnumactive (line 106)
    pub fn getnumactive(&self) -> HighsInt {
        self.active_constraint_index.len() as HighsInt
    }

    /// HiGHS: basis.hpp::Basis::getnuminactive (lines 108-110)
    pub fn getnuminactive(&self) -> HighsInt {
        self.non_active_constraint_index.len() as HighsInt
    }

    /// HiGHS: basis.hpp::Basis::getactive (lines 112-114)
    pub fn getactive(&self) -> &Vec<HighsInt> {
        &self.active_constraint_index
    }

    /// HiGHS: basis.hpp::Basis::getinactive (lines 116-118)
    pub fn getinactive(&self) -> &Vec<HighsInt> {
        &self.non_active_constraint_index
    }

    /// HiGHS: basis.hpp::Basis::getindexinfactor (lines 120-122)
    pub fn getindexinfactor(&self) -> &Vec<HighsInt> {
        &self.constraintindexinbasisfactor
    }

    /// HiGHS: basis.hpp::Basis::getstatus (line 124)
    pub fn getstatus(&mut self, conid: HighsInt) -> BasisStatus {
        *self
            .basisstatus
            .entry(conid)
            .or_insert(BasisStatus::kInactive)
    }

    /// HiGHS: basis.cpp::Basis::Basis (lines 13-42)
    pub fn new(
        rt: &mut Runtime,
        active: Vec<HighsInt>,
        status: Vec<BasisStatus>,
        inactive: Vec<HighsInt>,
    ) -> Self {
        let num_var = rt.instance.num_var;

        let mut basis = Basis {
            buffer_vec2hvec: HVector::default(),
            Ztprod_res: QpVector::new(num_var),
            buffer_Zprod: QpVector::new(num_var),

            basisfactor: HFactor::default(),
            updatessinceinvert: 0,

            Atran: MatrixBase::default(),

            active_constraint_index: Vec::new(),
            non_active_constraint_index: Vec::new(),
            baseindex: Vec::new(),

            basisstatus: BTreeMap::new(),

            constraintindexinbasisfactor: Vec::new(),

            buffer_column_aq: QpVector::new(num_var),
            buffer_row_ep: QpVector::new(num_var),

            buffered_q: -1,
            buffered_p: -1,
            row_ep: HVector::default(),
            col_aq: HVector::default(),

            reinversion_hint: false,
            h16_hfactor_pivot: rt.settings.h16_hfactor_pivot,
        };

        basis.buffer_vec2hvec.setup(num_var);

        for i in 0..(rt.instance.num_var + rt.instance.num_con) {
            basis.basisstatus.insert(i, BasisStatus::kInactive);
        }

        for i in 0..active.len() {
            basis.active_constraint_index.push(active[i]);
            basis
                .basisstatus
                .insert(basis.active_constraint_index[i], status[i]);
        }
        for i in 0..inactive.len() {
            basis.non_active_constraint_index.push(inactive[i]);
            basis.basisstatus.insert(
                basis.non_active_constraint_index[i],
                BasisStatus::kInactiveInBasis,
            );
        }

        basis.Atran = rt.instance.A.t().clone();

        basis.col_aq.setup(num_var);
        basis.row_ep.setup(num_var);

        basis.build();

        basis
    }

    /// HiGHS: basis.cpp::Basis::build (lines 44-76)
    pub(crate) fn build(&mut self) {
        // C++ line 45: report(); commented out.

        self.updatessinceinvert = 0;

        self.baseindex.resize(
            self.active_constraint_index.len() + self.non_active_constraint_index.len(),
            0,
        );
        self.constraintindexinbasisfactor.clear();

        self.basisfactor = HFactor::default();
        // HiGHS 1.16 #3354 (h16_hfactor_pivot): C++ `HFactor basisfactor;` is the 1.16.0 kernel
        self.basisfactor.h16_hfactor_pivot = self.h16_hfactor_pivot;

        self.constraintindexinbasisfactor.clear();
        self.constraintindexinbasisfactor
            .resize((self.Atran.num_row + self.Atran.num_col) as usize, -1);
        debug_assert!(
            (self.non_active_constraint_index.len() + self.active_constraint_index.len())
                as HighsInt
                == self.Atran.num_row
        );

        let mut counter: HighsInt = 0;
        for i in self.non_active_constraint_index.clone() {
            self.baseindex[counter as usize] = i;
            counter += 1;
        }
        for i in self.active_constraint_index.clone() {
            self.baseindex[counter as usize] = i;
            counter += 1;
        }

        self.basisfactor.setup_raw(
            self.Atran.num_col,
            self.Atran.num_row,
            &self.Atran.start,
            &self.Atran.index,
            &self.Atran.value,
            kDefaultPivotThreshold,
            kDefaultPivotTolerance,
            kHighsDebugLevelMin,
            None,
            true,
            kUpdateMethodFt,
        );
        self.basisfactor.build(&mut self.baseindex, None);

        for i in 0..(self.active_constraint_index.len() + self.non_active_constraint_index.len()) {
            self.constraintindexinbasisfactor[self.baseindex[i] as usize] = i as HighsInt;
        }
    }

    /// HiGHS: basis.cpp::Basis::rebuild (lines 78-96)
    pub fn rebuild(&mut self) {
        // C++ line 79: report(); commented out.

        self.updatessinceinvert = 0;
        self.constraintindexinbasisfactor.clear();

        self.constraintindexinbasisfactor.clear();
        self.constraintindexinbasisfactor
            .resize((self.Atran.num_row + self.Atran.num_col) as usize, -1);
        debug_assert!(
            (self.non_active_constraint_index.len() + self.active_constraint_index.len())
                as HighsInt
                == self.Atran.num_row
        );

        self.basisfactor.build(&mut self.baseindex, None);

        for i in 0..(self.active_constraint_index.len() + self.non_active_constraint_index.len()) {
            self.constraintindexinbasisfactor[self.baseindex[i] as usize] = i as HighsInt;
        }
        self.reinversion_hint = false;
    }

    /// HiGHS: basis.cpp::Basis::report (lines 98-220)
    pub fn report(&mut self) {
        // C++ lines 99-130: comment block, dropped.

        let qp_num_var: HighsInt = self.Atran.num_row;
        let qp_num_con: HighsInt = self.Atran.num_col;
        let num_active_in_basis: HighsInt = self.active_constraint_index.len() as HighsInt;
        let num_inactive_in_basis: HighsInt = self.non_active_constraint_index.len() as HighsInt;

        let mut num_var_inactive: HighsInt = 0;
        let mut num_var_active_at_lower: HighsInt = 0;
        let mut num_var_active_at_upper: HighsInt = 0;
        let mut num_var_inactive_in_basis: HighsInt = 0;
        let mut num_con_inactive: HighsInt = 0;
        let mut num_con_active_at_lower: HighsInt = 0;
        let mut num_con_active_at_upper: HighsInt = 0;
        let mut num_con_inactive_in_basis: HighsInt = 0;

        for i in 0..qp_num_var {
            match self.getstatus(qp_num_con + i) {
                BasisStatus::kInactive => {
                    num_var_inactive += 1;
                    continue;
                }
                BasisStatus::kActiveAtLower => {
                    num_var_active_at_lower += 1;
                    continue;
                }
                BasisStatus::kActiveAtUpper => {
                    num_var_active_at_upper += 1;
                    continue;
                }
                BasisStatus::kInactiveInBasis => {
                    num_var_inactive_in_basis += 1;
                    continue;
                }
            }
        }

        for i in 0..qp_num_con {
            match self.getstatus(i) {
                BasisStatus::kInactive => {
                    num_con_inactive += 1;
                    continue;
                }
                BasisStatus::kActiveAtLower => {
                    num_con_active_at_lower += 1;
                    continue;
                }
                BasisStatus::kActiveAtUpper => {
                    num_con_active_at_upper += 1;
                    continue;
                }
                BasisStatus::kInactiveInBasis => {
                    num_con_inactive_in_basis += 1;
                    continue;
                }
            }
        }

        let print_basis = num_inactive_in_basis + num_active_in_basis < 100;
        if print_basis {
            print!("basis: ");
            for a_ in self.active_constraint_index.clone() {
                if a_ < qp_num_con {
                    print!("c{:<3} ", a_);
                } else {
                    print!("v{:<3} ", a_ - qp_num_con);
                }
            }
            print!(" - ");
            for n_ in self.non_active_constraint_index.clone() {
                if n_ < qp_num_con {
                    print!("c{:<3} ", n_);
                } else {
                    print!("v{:<3} ", n_ - qp_num_con);
                }
            }
            println!();
        }

        print!(
            "Basis::report: QP({:>6} [inact {:>6}; act {:>6}], {:>6})",
            qp_num_var, num_inactive_in_basis, num_active_in_basis, qp_num_con
        );
        println!(
            " (inact / lo / up / basis) for var ({:>6} / {:>6} / {:>6} / {:>6}) and con ({:>6} / {:>6} / {:>6} / {:>6})",
            num_var_inactive,
            num_var_active_at_lower,
            num_var_active_at_upper,
            num_var_inactive_in_basis,
            num_con_inactive,
            num_con_active_at_lower,
            num_con_active_at_upper,
            num_con_inactive_in_basis
        );
        debug_assert!(qp_num_var == num_inactive_in_basis + num_active_in_basis);
        debug_assert!(qp_num_con == num_var_inactive + num_con_inactive);
        debug_assert!(
            num_inactive_in_basis == num_var_inactive_in_basis + num_con_inactive_in_basis
        );
        debug_assert!(
            num_active_in_basis
                == num_var_active_at_lower
                    + num_var_active_at_upper
                    + num_con_active_at_lower
                    + num_con_active_at_upper
        );
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

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10 and, when
    /// `with_row`, x0 + x1 <= 10. Optimum x = (1, 1), objective -3.
    fn box_qp(with_row: bool) -> Instance {
        let nc = if with_row { 1 } else { 0 };
        let mut inst = Instance::new(2, nc);
        inst.c = QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![-2.0, -4.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 1],
            value: vec![2.0, 4.0],
        };
        inst.A.mat = if with_row {
            MatrixBase {
                num_row: 1,
                num_col: 2,
                start: vec![0, 1, 2],
                index: vec![0, 0],
                value: vec![1.0, 1.0],
            }
        } else {
            MatrixBase {
                num_row: 0,
                num_col: 2,
                start: vec![0, 0, 0],
                index: vec![],
                value: vec![],
            }
        };
        inst.con_lo = vec![-kHighsInf; nc as usize];
        inst.con_up = vec![10.0; nc as usize];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variables active at their lower bound (constraint ids
    /// `num_con + i`); the row (id 0) is not in the basis.
    fn at_lower(rt: &mut Runtime) -> Basis {
        // C++ basis.cpp:27-30 reads status[i] for active[i], so `status`
        // has one entry per active index, not per constraint id.
        let status = vec![BasisStatus::kActiveAtLower, BasisStatus::kActiveAtLower];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kInactiveInBasis,
            BasisStatus::kInactiveInBasis,
        ];
        Basis::new(rt, vec![], status, vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn construct_counts_and_status() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut b = all_free(&mut rt);
        assert_eq!((b.getnumactive(), b.getnuminactive()), (0, 2));
        assert_eq!(b.getstatus(0), BasisStatus::kInactive);
        assert_eq!(b.getstatus(2), BasisStatus::kInactiveInBasis);
        let mut b = at_lower(&mut rt);
        assert_eq!((b.getnumactive(), b.getnuminactive()), (2, 0));
        assert_eq!(b.getstatus(1), BasisStatus::kActiveAtLower);
        assert!(!b.getreinversionhint());
    }
}
