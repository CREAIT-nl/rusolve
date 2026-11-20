//! HiGHS source: qpsolver/basis.cpp (lines 224-425)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! See `mod.rs` for the ownership notes. `Basis::write`
//! (declared in basis.hpp line 150) is commented out in basis.cpp
//! (lines 427-443) and not ported.

use crate::prelude::*;
use crate::qpsolver::*;

impl Basis {
    /// HiGHS: basis.cpp::Basis::deactivate (lines 224-230)
    pub fn deactivate(&mut self, conid: HighsInt) {
        // C++ line 225: printf commented out, dropped.
        debug_assert!(crate::qpsolver::snippets::contains(
            &self.active_constraint_index,
            &conid
        ));
        self.basisstatus
            .insert(conid, BasisStatus::kInactiveInBasis);
        crate::qpsolver::snippets::remove(&mut self.active_constraint_index, &conid);
        self.non_active_constraint_index.push(conid);
    }

    /// HiGHS: basis.cpp::Basis::activate (lines 232-259)
    pub fn activate(
        &mut self,
        rt: &Runtime,
        conid: HighsInt,
        newstatus: BasisStatus,
        nonactivetoremove: HighsInt,
        pricing: &mut dyn Pricing,
    ) -> QpSolverStatus {
        // C++ line 235: printf commented out, dropped.
        if !crate::qpsolver::snippets::contains(&self.active_constraint_index, &(conid as HighsInt))
        {
            self.basisstatus
                .insert(nonactivetoremove, BasisStatus::kInactive);
            self.basisstatus.insert(conid, newstatus);
            self.active_constraint_index.push(conid);
        } else {
            println!("Degeneracy? constraint {} already in basis", conid);
            return QpSolverStatus::DEGENERATE;
        }

        // C++ line 246: printf commented out, dropped.
        // C++ line 247: comment.
        let rowtoremove = self.constraintindexinbasisfactor[nonactivetoremove as usize];

        self.baseindex[rowtoremove as usize] = conid;
        crate::qpsolver::snippets::remove(
            &mut self.non_active_constraint_index,
            &nonactivetoremove,
        );
        self.updatebasis(rt, conid, nonactivetoremove, pricing);

        if self.updatessinceinvert != 0 {
            self.constraintindexinbasisfactor[nonactivetoremove as usize] = -1;
            self.constraintindexinbasisfactor[conid as usize] = rowtoremove;
        }
        QpSolverStatus::OK
    }

    /// HiGHS: basis.cpp::Basis::updatebasis (lines 261-294)
    pub fn updatebasis(
        &mut self,
        rt: &Runtime,
        newactivecon: HighsInt,
        droppedcon: HighsInt,
        pricing: &mut dyn Pricing,
    ) {
        if newactivecon == droppedcon {
            return;
        }

        const kHintNotChanged: HighsInt = 99999;
        let mut hint: HighsInt = kHintNotChanged;

        let droppedcon_rowindex = self.constraintindexinbasisfactor[droppedcon as usize];
        if self.buffered_p != droppedcon {
            self.row_ep.clear();
            self.row_ep.pack_flag = true;
            self.row_ep.index[0] = droppedcon_rowindex;
            self.row_ep.array[droppedcon_rowindex as usize] = 1.0;
            self.row_ep.count = 1;
            self.basisfactor.btran_call(&mut self.row_ep, 1.0, None);
        }

        // Compute hvec2vec(col_aq) and hvec2vec(row_ep) into locals first,
        // then call pricing.update_weights (see mod.rs port notes).
        let aq = Self::hvec2vec(&self.col_aq);
        let ep = Self::hvec2vec(&self.row_ep);
        pricing.update_weights(rt, self, &aq, &ep, droppedcon, newactivecon);

        let mut row_out = droppedcon_rowindex;

        // C++ HFactor holds a pointer to baseindex; passed here (lib.rs rule (c)).
        self.basisfactor.update(
            &mut self.baseindex,
            &mut self.col_aq,
            &mut self.row_ep,
            &mut row_out,
            &mut hint,
        );

        self.updatessinceinvert += 1;
        if self.updatessinceinvert >= rt.settings.reinvertfrequency || hint != kHintNotChanged {
            self.reinversion_hint = true;
        }
        // since basis changed, buffered values are no longer valid
        self.buffered_p = -1;
        self.buffered_q = -1;
    }

    /// HiGHS: basis.cpp::Basis::btran(const QpVector&, QpVector&, bool, HighsInt) (lines 296-311)
    ///
    /// C++ default arguments: buffer = false, p = -1
    pub fn btran_into<'a>(
        &mut self,
        rhs: &QpVector,
        target: &'a mut QpVector,
        buffer: bool,
        p: HighsInt,
    ) -> &'a mut QpVector {
        // C++ `HVector rhs_hvec = vec2hvec(rhs);` copies the member buffer,
        // so `buffer_vec2hvec` keeps its pre-btran contents.
        let mut rhs_hvec = self.vec2hvec(rhs).clone();
        self.basisfactor.btran_call(&mut rhs_hvec, 1.0, None);
        if buffer {
            self.row_ep.copy(&rhs_hvec);
            for i in 0..rhs_hvec.pack_count {
                self.row_ep.pack_index[i as usize] = rhs_hvec.pack_index[i as usize];
                self.row_ep.pack_value[i as usize] = rhs_hvec.pack_value[i as usize];
            }
            self.row_ep.pack_count = rhs_hvec.pack_count;
            self.row_ep.pack_flag = rhs_hvec.pack_flag;
            self.buffered_q = p;
        }
        Self::hvec2vec_into(&rhs_hvec, target)
    }

    /// HiGHS: basis.cpp::Basis::btran(const QpVector&, bool, HighsInt) (lines 313-327)
    ///
    /// C++ default arguments: buffer = false, p = -1
    pub fn btran(&mut self, rhs: &QpVector, buffer: bool, p: HighsInt) -> QpVector {
        // C++ copies the member buffer (see `btran_into`).
        let mut rhs_hvec = self.vec2hvec(rhs).clone();
        self.basisfactor.btran_call(&mut rhs_hvec, 1.0, None);
        if buffer {
            self.row_ep.copy(&rhs_hvec);
            for i in 0..rhs_hvec.pack_count {
                self.row_ep.pack_index[i as usize] = rhs_hvec.pack_index[i as usize];
                self.row_ep.pack_value[i as usize] = rhs_hvec.pack_value[i as usize];
            }
            self.row_ep.pack_count = rhs_hvec.pack_count;
            self.row_ep.pack_flag = rhs_hvec.pack_flag;
            self.buffered_q = p;
        }
        Self::hvec2vec(&rhs_hvec)
    }

    /// HiGHS: basis.cpp::Basis::ftran(const QpVector&, QpVector&, bool, HighsInt) (lines 329-345)
    ///
    /// C++ default arguments: buffer = false, q = -1
    pub fn ftran_into<'a>(
        &mut self,
        rhs: &QpVector,
        target: &'a mut QpVector,
        buffer: bool,
        q: HighsInt,
    ) -> &'a mut QpVector {
        // C++ copies the member buffer (see `btran_into`).
        let mut rhs_hvec = self.vec2hvec(rhs).clone();
        self.basisfactor.ftran_call(&mut rhs_hvec, 1.0, None);
        if buffer {
            self.col_aq.copy(&rhs_hvec);
            for i in 0..rhs_hvec.pack_count {
                self.col_aq.pack_index[i as usize] = rhs_hvec.pack_index[i as usize];
                self.col_aq.pack_value[i as usize] = rhs_hvec.pack_value[i as usize];
            }
            self.col_aq.pack_count = rhs_hvec.pack_count;
            self.col_aq.pack_flag = rhs_hvec.pack_flag;
            self.buffered_q = q;
        }

        Self::hvec2vec_into(&rhs_hvec, target)
    }

    /// HiGHS: basis.cpp::Basis::ftran(const QpVector&, bool, HighsInt) (lines 347-361)
    ///
    /// C++ default arguments: buffer = false, q = -1
    pub fn ftran(&mut self, rhs: &QpVector, buffer: bool, q: HighsInt) -> QpVector {
        // C++ copies the member buffer (see `btran_into`).
        let mut rhs_hvec = self.vec2hvec(rhs).clone();
        self.basisfactor.ftran_call(&mut rhs_hvec, 1.0, None);
        if buffer {
            self.col_aq.copy(&rhs_hvec);
            for i in 0..rhs_hvec.pack_count {
                self.col_aq.pack_index[i as usize] = rhs_hvec.pack_index[i as usize];
                self.col_aq.pack_value[i as usize] = rhs_hvec.pack_value[i as usize];
            }
            self.col_aq.pack_count = rhs_hvec.pack_count;
            self.col_aq.pack_flag = rhs_hvec.pack_flag;
            self.buffered_q = q;
        }
        Self::hvec2vec(&rhs_hvec)
    }

    /// HiGHS: basis.cpp::Basis::recomputex (lines 363-396)
    pub fn recomputex(&mut self, inst: &Instance) -> QpVector {
        debug_assert!(self.active_constraint_index.len() as HighsInt == inst.num_var);
        let mut rhs = QpVector::new(inst.num_var);

        for i in 0..inst.num_var {
            let con = self.active_constraint_index[i as usize];
            if self.constraintindexinbasisfactor[con as usize] == -1 {
                println!("error");
            }
            if self.getstatus(con) == BasisStatus::kActiveAtLower {
                if con < inst.num_con {
                    rhs.value[self.constraintindexinbasisfactor[con as usize] as usize] =
                        inst.con_lo[con as usize];
                } else {
                    rhs.value[self.constraintindexinbasisfactor[con as usize] as usize] =
                        inst.var_lo[(con - inst.num_con) as usize];
                }
            } else {
                if con < inst.num_con {
                    rhs.value[self.constraintindexinbasisfactor[con as usize] as usize] =
                        inst.con_up[con as usize];
                    // C++ line 382: commented-out, dropped.
                } else {
                    rhs.value[self.constraintindexinbasisfactor[con as usize] as usize] =
                        inst.var_up[(con - inst.num_con) as usize];
                    // C++ line 386: commented-out, dropped.
                }
            }

            rhs.index[i as usize] = i;
            rhs.num_nz += 1;
        }
        // C++ copies the member buffer (see `btran_into`).
        let mut rhs_hvec = self.vec2hvec(&rhs).clone();
        self.basisfactor.btran_call(&mut rhs_hvec, 1.0, None);
        Self::hvec2vec(&rhs_hvec)
    }

    /// HiGHS: basis.cpp::Basis::Ztprod (lines 398-411)
    ///
    /// C++ default arguments: buffer = false, q = -1
    pub fn Ztprod<'a>(
        &mut self,
        rhs: &QpVector,
        target: &'a mut QpVector,
        buffer: bool,
        q: HighsInt,
    ) -> &'a mut QpVector {
        // C++ line 400: ftran(rhs, Ztprod_res, buffer, q);
        // Ztprod_res is a member. Use std::mem::take (crate rule (f)).
        let mut ztprod_res = std::mem::take(&mut self.Ztprod_res);
        self.ftran_into(rhs, &mut ztprod_res, buffer, q);
        self.Ztprod_res = ztprod_res;

        target.reset();
        for i in 0..self.non_active_constraint_index.len() {
            let nonactive = self.non_active_constraint_index[i];
            let idx = self.constraintindexinbasisfactor[nonactive as usize];
            target.index[i] = i as HighsInt;
            target.value[i] = self.Ztprod_res.value[idx as usize];
        }
        target.resparsify();
        target
    }

    /// HiGHS: basis.cpp::Basis::Zprod (lines 413-425)
    pub fn Zprod<'a>(&mut self, rhs: &QpVector, target: &'a mut QpVector) -> &'a mut QpVector {
        // C++ line 414: buffer_Zprod.reset();
        // buffer_Zprod is a member. Use std::mem::take (crate rule (f)).
        let mut buffer_zprod = std::mem::take(&mut self.buffer_Zprod);
        buffer_zprod.reset();
        buffer_zprod.dim = target.dim;
        for i in 0..rhs.num_nz {
            let nz = rhs.index[i as usize];
            let nonactive = self.non_active_constraint_index[nz as usize];
            let idx = self.constraintindexinbasisfactor[nonactive as usize];
            buffer_zprod.index[i as usize] = idx;
            buffer_zprod.value[idx as usize] = rhs.value[nz as usize];
        }
        buffer_zprod.num_nz = rhs.num_nz;

        // C++ line 424: return btran(buffer_Zprod, target);
        // btran here is the overload that takes (rhs, target, buffer, p).
        // buffer defaults to false, p defaults to -1.
        let result = self.btran_into(&buffer_zprod, target, false, -1);

        self.buffer_Zprod = buffer_zprod;

        result
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
        // C++ Basis::Basis indexes `status` by position in `active`
        // (basis.cpp lines 27-29): one status per active constraint.
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
    fn recompute_x_and_deactivate() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let inst = rt.instance.clone();
        let mut b = at_lower(&mut rt);
        assert!(close(&b.recomputex(&inst).value, &[0.0, 0.0]));
        b.deactivate(2);
        assert_eq!((b.getnumactive(), b.getnuminactive()), (1, 1));
        assert_eq!(b.getstatus(2), BasisStatus::kInactiveInBasis);
    }

    /// Pricing stub: counts `update_weights` calls; never prices.
    struct CountPricing(HighsInt);
    impl Pricing for CountPricing {
        fn price(
            &mut self,
            _rt: &Runtime,
            _basis: &mut Basis,
            _redcosts: &mut ReducedCosts,
            _gradient: &mut Gradient,
            _x: &QpVector,
        ) -> HighsInt {
            unreachable!()
        }
        fn update_weights(
            &mut self,
            _rt: &Runtime,
            _basis: &mut Basis,
            _aq: &QpVector,
            _ep: &QpVector,
            _p: HighsInt,
            _q: HighsInt,
        ) {
            self.0 += 1;
        }
        fn recompute(&mut self, _rt: &Runtime, _basis: &mut Basis) {
            unreachable!()
        }
    }

    fn vec2(a: f64, b: f64) -> QpVector {
        QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![a, b],
        }
    }

    /// Slack-only basis (ids 1, 2 are the column bounds): B = I, so
    /// ftran / btran are the identity; `buffer` stores into col_aq / row_ep
    /// and sets `buffered_q` (basis.cpp lines 308, 324, 341, 358).
    #[test]
    fn ftran_btran_identity_and_buffering() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut b = all_free(&mut rt);
        let v = vec2(3.0, -5.0);
        assert!(close(&b.ftran(&v, false, 7).value, &[3.0, -5.0]));
        assert_eq!(b.buffered_q, -1);
        assert!(close(&b.ftran(&v, true, 7).value, &[3.0, -5.0]));
        assert_eq!(b.buffered_q, 7);
        assert_eq!(b.col_aq.count, 2);
        assert!(close(&b.col_aq.array[..2], &[3.0, -5.0]));
        let mut t = QpVector::new(2);
        assert!(close(
            &b.btran_into(&v, &mut t, true, 4).value,
            &[3.0, -5.0]
        ));
        assert_eq!(b.buffered_q, 4);
        assert!(close(&b.row_ep.array[..2], &[3.0, -5.0]));
        let mut t = QpVector::new(2);
        assert!(close(
            &b.ftran_into(&v, &mut t, false, -1).value,
            &[3.0, -5.0]
        ));
        assert!(close(&b.btran(&v, false, -1).value, &[3.0, -5.0]));
        // C++ copies vec2hvec's buffer before the solve: the member keeps the input.
        assert_eq!(b.buffer_vec2hvec.count, 2);
        assert!(close(&b.buffer_vec2hvec.array[..2], &[3.0, -5.0]));
    }

    /// Z = I for the all-free basis: Ztprod and Zprod return the input.
    #[test]
    fn z_products_identity() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut b = all_free(&mut rt);
        let v = vec2(2.0, 0.5);
        let mut t = QpVector::new(2);
        let r = b.Ztprod(&v, &mut t, false, -1);
        assert_eq!(r.num_nz, 2);
        assert!(close(&r.value, &[2.0, 0.5]));
        let mut t = QpVector::new(2);
        assert!(close(&b.Zprod(&v, &mut t).value, &[2.0, 0.5]));
    }

    /// activate on an already-active constraint returns DEGENERATE before
    /// touching the basis (basis.cpp lines 240-244); updatebasis returns
    /// early when the entering and leaving ids match (lines 263-265).
    #[test]
    fn activate_degenerate_and_updatebasis_noop() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        let mut b = at_lower(&mut rt);
        let mut pricing = CountPricing(0);
        let st = b.activate(&rt, 1, BasisStatus::kActiveAtUpper, 2, &mut pricing);
        assert_eq!(st, QpSolverStatus::DEGENERATE);
        assert_eq!(b.getstatus(1), BasisStatus::kActiveAtLower);
        assert_eq!(b.getstatus(2), BasisStatus::kActiveAtLower);
        b.updatebasis(&rt, 2, 2, &mut pricing);
        assert_eq!((pricing.0, b.getnupdatessinceinvert()), (0, 0));
    }
}
