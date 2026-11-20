# benches/data — instance attribution

All instance data here is fetched, never committed — with one exception:
`netlib/` is committed (~1.9 MB). `miplib`, `lp47`, `lpopt`, `mip50` are
fetch-only (`tools/bench/fetch_corpus.py` fetches the MIPLIB instances; `lpopt/`
and `mip50/` are downloaded from the sources below, `lp47/` is derived from them).

| Corpus | Source | Attribution / license |
|---|---|---|
| `netlib/` | Netlib LP test set via github.com/ozy4dm/lp-data-netlib (emps-expanded) | D. Gay / J. Hall collection, https://netlib.org/lp/data/ |
| `lpopt/` | https://plato.asu.edu/ftp/lptestset/ (+ misc/ nug/ pds/ rail/ fome/) | Robert E. Bixby's LP test set, maintained and published by Hans D. Mittelmann (ASU); MPC files expanded with netlib `emps` (emp.c, netlib, https://www.netlib.org/lp/data/) |
| `mip50/`, `miplib/` | https://miplib.zib.de/WebData/instances/ | MIPLIB 2017, Gleixner, Gomeshari, Salani, Wiggers, "MIPLIB 2017: data-driven compilation of the 6th mixed-integer programming library", Math. Prog. Comp. 13 (2021); instance files CC BY-SA 4.0 by default — each instance page names its original contributor and any instance-specific license |
| `lp47/` | LP relaxations of the 47 pinned MIPLIB 2017 instances | as `miplib/` above (derived works) |
