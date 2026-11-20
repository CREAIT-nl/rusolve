# Third-party notices

Portions of rusolve are translated from HiGHS (https://github.com/ERGO-Code/HiGHS),
including its bundled interior point solver IPX and the BASICLU LU factorization
(Copyright (C) 2016-2019 ERGO-Code, Lukas Schork), all distributed with HiGHS and
used under the MIT License:

```
MIT License

Copyright (c) 2026 HiGHS

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

`crates/rusolve-ekk/src/util/highs_pdqsort.rs` is an altered (translated to Rust)
version of pdqsort, vendored by HiGHS at `extern/pdqsort/pdqsort.h`, used under the
zlib License:

```
pdqsort.h - Pattern-defeating quicksort.

Copyright (c) 2021 Orson Peters

This software is provided 'as-is', without any express or implied warranty. In no event will the
authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial
applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the
   original software. If you use this software in a product, an acknowledgment in the product
   documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as
   being the original software.

3. This notice may not be removed or altered from any source distribution.
```

`crates/rusolve-ekk/src/util/cpp_std.rs` is translated from the libc++ algorithm headers
(`<__algorithm/push_heap.h`, `pop_heap.h`, `make_heap.h`, `sift_down.h`, `sort_heap.h`,
`partial_sort.h`, `partition.h`; libc++ 210106, from the macOS SDK). The `Mt19937` type in
`crates/rusolve-ekk/src/mip/mod.rs` is translated from libc++ `<random>`
(`mersenne_twister_engine` seeding and twist, and `generate_canonical<double, 53>` behind
`uniform_real_distribution`) from the same SDK, under the same license. libc++ is part of the
LLVM Project, used under the Apache License v2.0 with LLVM Exceptions (SPDX:
`Apache-2.0 WITH LLVM-exception`). The full license text is at https://llvm.org/LICENSE.txt:

```
Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
See https://llvm.org/LICENSE.txt for license information.
SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
```

## HiPO orderings

The fill-reducing orderings under `crates/rusolve-ekk/src/ipm/hipo/orderings/` are
translated (to Rust, trimmed to the routines HiPO reaches) from the copies HiGHS vendors
in `extern/`. Each translated file carries a header naming its origin.

`crates/rusolve-ekk/src/ipm/hipo/orderings/metis/` is translated from METIS 5.2.1 and the
parts of GKlib it uses (`extern/metis`), used under the Apache License 2.0
(https://www.apache.org/licenses/LICENSE-2.0). The translation is a modified version: the
files were translated to Rust and trimmed to the `METIS_NodeND` call graph.

```
Copyright & License Notice
---------------------------

Copyright 1995-2013, Regents of the University of Minnesota

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or 
implied. See the License for the specific language governing 
permissions and limitations under the License.
```

GKlib's `gk_mksort.h` (the `GK_MKQSORT` template behind `isortd` / `ikvsorti`) states that
it is adopted from GNU glibc `stdlib/qsort.c` under the GNU Lesser General Public License
2.1 or later. It is not translated and no part of it is included in the crates: the sort in
`orderings/metis/gklib.rs` is an independent implementation written from a description of
the sort's behaviour, by an implementer who did not see the LGPL source. `tests/t_gksort.rs` checks
it against reference outputs (data only) produced by `tools/highs-oracle/gksort`, a test
tool that compiles the vendored header and is not part of any crate.

`crates/rusolve-ekk/src/ipm/hipo/orderings/amd/` is translated from SuiteSparse AMD
(`extern/amd`), used under the BSD 3-clause license:

```
AMD, Copyright (c), 1996-2022, Timothy A. Davis,
Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.

Availability:

    http://www.suitesparse.com

-------------------------------------------------------------------------------
AMD License: BSD 3-clause:
-------------------------------------------------------------------------------

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of the organizations to which the authors are
          affiliated, nor the names of its contributors may be used to endorse
          or promote products derived from this software without specific prior
          written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.
```

`crates/rusolve-ekk/src/ipm/hipo/orderings/rcm/` is translated from the reverse
Cuthill-McKee code vendored by HiGHS (`extern/rcm/rcm.cpp`), used under the MIT License:

```
Copyright (c) 2005 Alan George, Joseph Liu, John Burkardt.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## SCIP conflict resolution

`crates/rusolve-ekk/src/mip/conflict_resolution/` is translated from SCIP 11.0.0
(commit 32dd425; `src/scip/conflict_resolution.c`, `def.h`, `set.c`), used under the
Apache License 2.0 (https://www.apache.org/licenses/LICENSE-2.0). The translation is a
modified version: the generalized-resolution conflict analysis is ported to Rust for HiGHS
bound-change stacks and runs only behind the default-off `mip_conflict_resolution` option.

```
Copyright & License Notice
---------------------------

Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing
permissions and limitations under the License.
```

## SCIP cut selection

`crates/rusolve-ekk/src/mip/cut_selection/` is translated from SCIP 11.0.0
(commit 32dd425; `src/scip/cutsel_hybrid.c`, `cutsel.c`, `sepastore.c`, `lp.c`, `misc.c`,
`set.c`, `def.h`), used under the Apache License 2.0 (https://www.apache.org/licenses/LICENSE-2.0).
The translation is a modified version: SCIP's hybrid cut selector and sepastore apply logic are
ported to Rust on HiGHS cut pools and run only behind the default-off `mip_cut_selection_hybrid`
option. `pool_s4.rs` also copies HiGHS `HighsCutPool.cpp` (MIT License, see the HiGHS section above).

```
Copyright & License Notice
---------------------------

Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing
permissions and limitations under the License.
```

## SCIP LNS heuristics

`crates/rusolve-ekk/src/mip/scip_heuristics/` is translated from SCIP 11.0.0
(commit 32dd425; `src/scip/heur_oneopt.c`, `heur_crossover.c`, `primal.c`, `sol.c`,
`misc.c`, `sorttpl.c`, `heur.c`, `set.c`, `scip_randnumgen.c`, `scip_copy.c`, `def.h`),
used under the Apache License 2.0 (https://www.apache.org/licenses/LICENSE-2.0). The
translation is a modified version: the oneopt and crossover primal heuristics are ported
to Rust for HiGHS model and solution structures and run only behind the default-off
`mip_heuristic_run_oneopt` and `mip_heuristic_run_crossover` options. The
`sub_mip_s6.rs` module also copies HiGHS `HighsPrimalHeuristics.cpp` (MIT License).

`proximity.rs` and `local_branching.rs` in the same directory are translated from SCIP
(commit 63417ae; `src/scip/heur_proximity.c`, `heur_localbranching.c`; same Apache License 2.0
notice below) and run only behind the default-off `mip_proximity` and `mip_local_branching`
options. They are modified: the sub-problem is a HiGHS sub-MIP; SCIP's LP-iteration limit is
expressed as a work limit, and the sub-SCIP settings are mapped to the sub-MIP options that exist.

`guided_dive.rs` is translated from SCIP (`src/scip/heur_guideddiving.c`, `heur_dive.c`; same Apache License 2.0
notice below) and runs only behind the default-off `mip_guided_dive` option. It is modified: it runs on
the HiGHS heuristic search and LP copy, and all candidates are scored as unroundable.

```
Copyright & License Notice
---------------------------

Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing
permissions and limitations under the License.
```

`scip_heuristics/shiftprop.rs` (the default-off `mip_shift_and_propagate` option) is translated
from SCIP `src/scip/heur_shiftandpropagate.c` and `sorttpl.c` (scipopt/scip commit 63417ae, same
Apache License 2.0 and copyright notice as above). The translation is a modified version: the
SCIP probing mode, bound-change events and probing LP are replaced by a local `HighsDomain` copy
synchronised from its change stack and by the port's `try_rounded_point`; propagation runs to
the fixpoint of `HighsDomain::propagate` instead of 10 rounds.

`mip/knapsack_cover_separator.rs` (the default-off `mip_knapsack_cover_sepa` option) is translated
from SCIP `src/scip/cons_knapsack.c` and the row selection of `cons_linear.c::separateCons`
(scipopt/scip commit 96b1a0d, same Apache License 2.0 and copyright notice as above):
`SCIPseparateRelaxedKnapsack`, `SCIPseparateKnapsackCuts`, `SCIPsolveKnapsackApproximately`,
`getCover`, `makeCoverMinimal`, `getFeasibleSet`, `sequentialUpAndDownLifting` and the
`_mcseq` / `_ewseq` separation functions. It is a modified version: SCIP's variable bounds are
read from the HiGHS implication tables, GUB lifting and superadditive lifting (both off by default
in SCIP) are not ported, and the partial weighted-median sort is a full sort.

## mip-heuristics FPR / LocalMIP

`crates/rusolve-ekk/src/mip/fpr_local_mip/` is translated from mip-heuristics
(https://github.com/spoorendonk/mip-heuristics, commit 36f7d9d), used under the MIT License.
The translation is a modified version: the sequential FPR and LocalMIP presolve heuristics are
ported to Rust and run only behind the default-off `mip_heuristic_run_fpr` and
`mip_heuristic_run_local_mip` options. `context.rs` and `sink.rs` also read HiGHS
`HighsMipSolver` state (MIT License, see the HiGHS section above).

Files translated (under `src/`): `clique_cover.cpp/.h`, `continuous_loop.h`, `deadline.h`,
`effort_ledger.cpp/.h`, `fpr.cpp/.h`, `fpr_core.cpp/.h`, `fpr_strategies.h`,
`fpr_val_select.cpp/.h`, `fpr_var_order.cpp/.h`, `heuristic_common.h`,
`heuristic_context.cpp/.h`, `highs_context.h`, `incumbent_sink.cpp/.h`, `least_violated.h`,
`local_mip.cpp/.h`, `local_mip_caches.h`, `local_mip_construction.cpp/.h`,
`local_mip_core.cpp/.h`, `local_mip_search.cpp`, `local_mip_worker.cpp/.h`,
`mode_dispatch.cpp/.h`, `opportunistic_runner.h`, `prop_engine.cpp/.h`,
`repair_search.cpp/.h`, `repair_walk.cpp/.h`, `rng.h`, `solution_pool.cpp/.h`,
`solution_sink.cpp/.h`, `walksat.cpp/.h`, `worker_base.h`. The test fixture
`crates/rusolve-ekk/tests/fixtures/s5/rng_oracle.cpp` includes `rng.h`.

```
MIT License

Copyright (c) 2026 Simon Spoorendonk

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

`rng.rs` re-expresses the algorithms of libstdc++ 13.3 (GCC 13.3.0)
`uniform_int_distribution`, `generate_canonical`, `uniform_real_distribution`,
`bernoulli_distribution`, `std::shuffle` and the `std::sort` introsort, so that draws match the
GCC-built upstream bit for bit. libstdc++ is licensed under the GPL-3.0 with the GCC Runtime
Library Exception. Only the algorithms are re-expressed; no libstdc++ source text is copied.

`crates/rusolve-ekk/src/util/cpp_std.rs` (module `gnu`) and `crates/rusolve-ekk/src/mip/fpr_local_mip/lm_search.rs`
(the heap functions near the top of the file) re-express the control flow of the libstdc++ 13.3
`bits/stl_heap.h` heap algorithms (`__adjust_heap`, `__pop_heap`, `__make_heap`, `__sort_heap`,
`__heap_select`, `__push_heap`) so that tie order matches a GCC-built upstream. Same license and
terms as above (GPL-3.0 with the GCC Runtime Library Exception). `cpp_std.rs::gnu` was written from the
author's recollection of the algorithm and checked against the output of a small g++ 13.3 test program; it mirrors the structure of those functions, so treat it
as a re-expression of the algorithm, not as independent work. No libstdc++ source text is copied.

## SCIP parity presolve and orbitope fixing

`crates/rusolve-ekk/src/presolve/hpresolve/xor_gf2.rs` and `orbitope_root_fix` in
`mip/mip_solver_data/part1.rs` re-implement the algorithms of SCIP 10 `cons_xor.c`
(`linconsUpgdXor`, `checkSystemGF2`) and `cons_orbitope.c` (Apache-2.0, https://scipopt.org).
They are written from the algorithm descriptions, not translated line by line. SCIP is
Copyright (c) Zuse Institute Berlin, licensed under the Apache License 2.0.

## SCIP full-orbitope propagation (reference only)

`mip_orbitope_full_fix` applies lexicographic orbitopal fixing at the
root. The routine itself is the HiGHS port (`orbital_fixing_for_full_orbitope`, MIT, above); the
idea of running it over all rows follows SCIP `cons_orbitope.c` (Apache-2.0, Zuse Institute Berlin)
and Kaibel, Peinhardt, Pfetsch, "Orbitopal fixing", Discrete Optimization 8 (2011). No SCIP code
was copied.

## MathOptFormat example (test fixture)

`crates/rusolve-io/tests/mof/milp.mof.json` is `examples/milp.mof.json` from
https://github.com/jump-dev/MathOptFormat, copied verbatim, used under the MIT License:

```
MIT License

Copyright (c) 2020 Oscar Dowson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Clarabel (linked dependency)

`crates/rusolve-cli` links the `clarabel` crate (0.11.1, https://github.com/oxfordcontrol/Clarabel.rs),
an interior-point solver for conic programs (Goulart and Chen 2024, arXiv:2405.12762), used
unmodified under the Apache License 2.0 (https://www.apache.org/licenses/LICENSE-2.0). It solves the
continuous SOCP and convex QCQP models routed by `crates/rusolve-cli/src/conic_solve.rs`.
No Clarabel source is copied into this repository.
