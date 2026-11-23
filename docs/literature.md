# Core literature — implementing a modern LP/MIP solver

✔ = free copy available · ◻ = paywalled (author/tech-report
versions often findable) · [FREE] = publisher open access.

## Overviews
- Huangfu & Hall (2018), *Parallelizing the dual revised simplex method* — the HiGHS simplex paper. ✔ arXiv:1503.01889
- Bestuzheva et al. (2025), *The HiGHS class of linear and mixed integer optimization solvers*, Math. Prog. Comp. ◻ (Springer s12532)
- Bestuzheva et al. (2026), *Presolve heuristics in HiGHS*. ✔ arXiv:2609.22938

## Simplex
- Gold & Price (1983), steepest edge ◻
- Forsgren, Gill & Wright (2015), non-exact **Devex** pricing, SIAM J. Optim. ◻ (free tech report on Gill's Stanford page)
- Kober, *Implementing Devex*, MSc thesis Freiburg ◻ (free university PDF)
- Forrest & Gold (2016), **branch-dual** pivot rule, Math. Prog. Comp. ◻
- Bixby (1991), *Empirical behavior of the primal simplex algorithm* — classic on degenerate LPs ◻
- Network simplex: Arbenz & Geissmann tutorial ◻; Kovács (2015) EJOR ◻

## Presolve
- **Bestuzheva, Kochenderfer, Achterberg et al., *Presolve reductions in mixed integer programming*** ✔ https://opus4.kobv.de/opus4-zib/files/6037/Presolve.pdf (MIPC 2023)
- Andersen & Gondzio (1999), *Sparse LP pre-processing* ◻
- Gleixner et al., *Multi-pass high-level presolving*, IJCAI'15 ✔ https://www.ijcai.org/Proceedings/15/Papers/055.pdf
- *Exploiting Variable Implications in Presolve for MIP* (2026) ✔ arXiv:2607.04313

## Interior point / first-order / crossover
- Wright (1997), *Primal-Dual Interior-Point Methods* (SIAM) ✔ https://www.mcs.anl.gov/research/projects/otc/InteriorPoint/
- Ma, O'Donoghue, Diamond, Boyd, **PDLP** ✔ arXiv:2106.04756; full paper ✔ arXiv:2501.07018
- GPU first-order methods survey ✔ arXiv:2506.02174; cuPDLP.jl ✔ arXiv:2311.12180
- Inexact IPM + incomplete LDLᵀ preconditioning ✔ arXiv:1708.04298
- O'Donoghue (2021), operator splitting for homogeneous LCP embedding (crossover-adjacent) ◻

## Branch & bound / cuts / heuristics
- Achterberg (2009), *SCIP: solving constraint integer programs*, MIPC ◻ [likely free via Springer]
- Achterberg & Berthold (2007), *Conflict analysis in MIP* ✔ ZIB report https://opus4.kobv.de/opus4-zib/files/853/ZR-05-19.pdf
- Achterberg, *Status of MIP cuts and separation* — recent slides ✔ https://www.mixedinteger.org/2025/slides/2025-06-04%20-%20Talks%20-%20Tobias%20Achterberg.pdf
- Gamrath (2010), *Generic cutting plane framework with parallel cut selection at the root node*, ZIB 10-32 ◻ (opus4 open-access)
- *Implementing cutting plane management and selection* (Mops) ✔ https://optimization-online.org/wp-content/uploads/2012-12-3714.pdf
- Achterberg & Berthold, *Rank & Select* ◻ (opus4)
- Fischetti & Lodi (2003), *Feasibility pump* ◻ (author page dei.unipd.it/~fisch/papers/)
- Berthold (2020), *Heuristic algorithms in commercial MIP solvers*, MIPC ◻ [likely OA]; RENS/RINS/undercover papers ◻
- Wolsey & Sammlung, *Cutting planes: a guided tour* ◻ (Wolsey's page); Wolsey, *Integer Programming* (author-distributed)
- MIPLIB 2017 (Amon et al., MIPC 2021) [FREE ✔]

## Parallelization
- Huangfu & Hall ✔ (above) — dual-simplex task design
- Rothberg, *MIP solver internals* (implicit parallelism in branch-and-bound) ✔ slides https://arpa-e.energy.gov/sites/default/files/5_Rothberg_presentation.pdf
- Boland & Dahal, *Parallel solvers for MILP* survey ✔ Lehigh copy http://coral.ie.lehigh.edu/~ted/files/papers/ParallelMILPSurvey16.pdf
- HiGHS `parallel/` task executor (~1.7k LOC) — the practical spec for our Rust task graph

## Permissively licensed solvers and libraries
- **HiGHS** (MIT) ✔ https://github.com/ERGO-Code/HiGHS — the solver rusolve ports subsystem by
  subsystem (simplex, IPX, HiPO, PDLP, presolve, MIP); attribution in `THIRD_PARTY_NOTICES.md`.
- **SCIP** (Apache-2.0) ✔ https://github.com/scipopt/scip — source of the ported conflict
  resolution, cut selection and LNS heuristics.
- **Clarabel** (Apache-2.0) ✔ https://github.com/oxfordcontrol/Clarabel.rs — Rust conic
  interior-point solver, linked for SOCP and convex QCQP (Goulart & Chen 2024, arXiv:2405.12762).
- **mip-heuristics** (MIT) ✔ https://github.com/spoorendonk/mip-heuristics — FPR and LocalMIP.
- METIS 5.2.1 (Apache-2.0) and SuiteSparse AMD (BSD-3) — the fill-reducing orderings HiPO uses.
- `highs`/`highs-sys` crates (rust-or, MIT) — in-process HiGHS bindings, useful as a
  differential-testing oracle; `sprs` (sparse linear algebra in Rust).
