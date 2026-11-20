#!/usr/bin/env python3
"""Fetch the benchmark corpora: Netlib LP (default) and the pinned MIPLIB 2017 slice.

    python tools/bench/fetch_corpus.py                  # Netlib -> benches/data/netlib/
    python tools/bench/fetch_corpus.py --set miplib     # MIPLIB slice -> benches/data/miplib/
    python tools/bench/fetch_corpus.py --set miplib-ref # all 47 pinned MIPLIB instances

Provenance
----------
Netlib LP test set (https://netlib.org/lp/data/, D. Gay / J. Hall collection).
The canonical site stores the instances in a *compressed* MPS format that
requires the `emps` program to expand; the historical `stand-inized/`
subdirectory is no longer served. We therefore fetch the emps-expanded plain
MPS files from https://github.com/ozy4dm/lp-data-netlib ("Netlib LP problem
data converted into readily usable MPS format", expanded from the compressed
netlib files — numerically identical instances, standard fixed/free MPS).

MIPLIB 2017 (https://miplib.zib.de/, Gleixner et al. 2021, "MIPLIB 2017:
data-driven compilation of the 6th mixed-integer programming library", Math.
Prog. Comp. 13): 47 pinned benchmark-set instances, each with the sha256 of its
.mps.gz and the published optimal objective from its instance page. The
`miplib` slice is a quick subset: instances HiGHS 1.15 single-thread
solves to optimality within 60 s with m (n + 3m) <= 4e6. The files are large and licensed per instance, so
benches/data/miplib/ is gitignored: fetch, never commit. A sha256 mismatch
is a hard failure (the published objective no longer applies).

Idempotent: existing files with a plausible MPS header (Netlib) or a
matching .mps.gz hash (MIPLIB) are kept.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import ssl
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO = "ozy4dm/lp-data-netlib"
BRANCH = "master"
URL = f"https://raw.githubusercontent.com/{REPO}/{BRANCH}/mps_files/{{name}}.mps"

# Small-to-medium Netlib instances, chosen for structural diversity:
# degenerate classics (afiro, sc50a, sc105), blending (blend), set-covering
# (sctap1), standard-form with free vars (scsd1, scsd6), ranged rows
# (share1b/share2b, boeing2), network-ish (capri, bore3d, stair), larger
# medium (scagr25, standmps, finnis).
INSTANCES = [
    "afiro", "adlittle", "blend", "boeing2", "bore3d", "brandy", "capri",
    "e226", "finnis", "grow7", "kb2", "lotfi", "sc105", "sc50a", "scagr7",
    "scagr25", "scsd1", "scsd6", "scorpion", "scrs8", "sctap1", "share1b",
    "share2b", "stair", "standata", "standmps", "tuff", "vtp.base",
]

DEST = Path(__file__).resolve().parents[2] / "benches" / "data" / "netlib"

MIPLIB_URL = "https://miplib.zib.de/WebData/instances/{name}.mps.gz"
MIPLIB_DEST = Path(__file__).resolve().parents[2] / "benches" / "data" / "miplib"

# (name, published optimal objective, rows, cols, in slice, why, sha256 of
# <name>.mps.gz). Objectives from the MIPLIB 2017 instance pages (all 47 are
# status "easy": optimum known). "HiGHS" columns: highspy 1.15.1, threads=1,
# mip_rel_gap=0, 60 s limit. "dense" = m (n + 3m), the doubles
# our dense LP allocates per node (to_standard + basis).
MIPLIB = [
    ("30n20b8", 302.0, 576, 18380, False, "cut: HiGHS >60 s, dense 11.6M > 4M",
     "61681c962ef830c2f958cb39bc03d6a1799f44158e5efeb97a1e6c3592f2ac40"),
    ("air05", 26374.0, 426, 7195, True, "HiGHS 16.4 s, 57 nodes; dense 3.61M",
     "8e03d2b90834dc1fe15c908247fdd8baf76e208aa310508caa0734e5245692cc"),
    ("assign1-5-8", 211.999999999998, 161, 156, False, "cut: HiGHS >60 s",
     "305c26c9eee8da7242f2f42504cb0988c9503f65972859ad5f26d9615086cde3"),
    ("beasleyC3", 753.9999999999128, 1750, 2500, False, "cut: dense 13.6M > 4M",
     "37350bc46bdb84a6373e9044eb3d0d65202427fe2bcba64ed282f894699d75dd"),
    ("binkar10_1", 6741.38002393972, 1026, 2298, False, "cut: dense 5.5M > 4M",
     "94282c38edd6c83bb6c2d788d823fd82f218218e8e07c98d1d1ab42c1e91039e"),
    ("bppc4-08", 53.0, 111, 1456, False, "cut: HiGHS >60 s",
     "f7441940af00d8a1a68b889f95dcabbe0185b68f53598cb70595fab23eca1cec"),
    ("cod105", -12.0, 1024, 1024, False, "cut: HiGHS >60 s, dense 4.2M > 4M",
     "f0007f47d356d70d1843e371ac75c7dbfdef4a2947785d47059e61821097e15a"),
    ("dano3_3", 576.34463303, 3202, 13873, False, "cut: dense 75.2M > 4M",
     "dfae83681f99a910db2f05018df647c21e8df6e09fd1d160b0a6dc29c76d15e7"),
    ("drayage-25-23", 101282.647018, 4630, 11090, False, "cut: dense 115.7M > 4M",
     "207c00284a4edbe92f1cc11705e47f03156830a804fbaf748e19f7c84103e131"),
    ("enlight_hard", 37.0, 100, 200, True, "HiGHS 6.6 s, 1986 nodes; dense 0.05M",
     "942168c2126a2a91ae3ec1ededea59bc1af0cad55f94223edf4c03d20e831f66"),
    ("fiball", 138.0, 3707, 34219, False, "cut: HiGHS >60 s, dense 168.1M > 4M",
     "ee5835347bf31af5a0541d1b9ca6db8458872be2892c67eb137056b843277d1f"),
    ("gen-ip002", -4783.733392, 24, 41, False, "cut: HiGHS >60 s",
     "cb37e723fa2187947981215bf9c8ba40d452984d7b9be3467e9f9d18f0343f78"),
    ("gen-ip054", 6840.96564179, 27, 30, False, "cut: HiGHS >60 s",
     "b5089ebc97c43bd583fcadf30955b8800861044d632d747c0b7ac73f19b9e6f3"),
    ("glass4", 1200012599.972384, 396, 322, False, "cut: HiGHS >60 s",
     "0d6e6ceaf921e5d5489b9926fdfed2ce0aa25fd9acfa9deb6794669b56b59591"),
    ("gmu-35-40", -2406733.3688, 424, 1205, False, "cut: HiGHS >60 s",
     "0c91e71117164becea8109a414127a2223be82215498f243767a7f5f51b1f11d"),
    ("graph20-20-1rand", -9.0, 5587, 2183, False, "cut: HiGHS >60 s, dense 105.8M > 4M",
     "8aa8510022fd286b61306bbfa165fe9ea47708aa8fd5d7d1c6a55f45750f5bd3"),
    ("h80x6320d", 6382.099048245999, 6558, 12640, False, "cut: dense 211.9M > 4M",
     "da72f2491ab4c7c058074b71a140e12af3c2e87e52e68161f88b2aca8a5f9a23"),
    ("icir97_tension", 6375.0, 1203, 2494, False, "cut: HiGHS >60 s, dense 7.3M > 4M",
     "b900d1eacbb02a57253f4bac839e0c1f3464a18a4686f8a621610d5d735e92ff"),
    ("istanbul-no-cutoff", 204.08170701, 20346, 5282, False, "cut: HiGHS >60 s, dense 1349.3M > 4M",
     "2b9171a1b1e4236c33ea444b24a98892fb0228ac9309495b4fac280e2ad5653a"),
    ("lotsize", 1480195.0, 1920, 2985, False, "cut: HiGHS >60 s, dense 16.8M > 4M",
     "14702c1793474b78a30b91512a12917790ba860f7abe765b01b296fb86835902"),
    ("markshare_4_0", 1.0, 4, 34, False, "cut: HiGHS >60 s",
     "1a1feb04f5637db9f69be53bd6224fcd9fe5b8d750870d9b2002b535d25f4353"),
    ("mas74", 11801.18572, 13, 151, False, "cut: HiGHS >60 s",
     "019a9ca2dac800342a887cd07e326526b463dca9646ba16ad955a16fa4b6ecc1"),
    ("mas76", 40005.05398999999, 12, 151, False, "cut: HiGHS >60 s",
     "91f624659e181538020e0c627910f02ffa6b355d6df708ef324e697a69485a74"),
    ("mc11", 11688.99999999966, 1920, 3040, False, "cut: dense 16.9M > 4M",
     "55f4d9c9d2022cffeb3e07cd7398753f5aab720f16fa49ab17223375bab83874"),
    ("mik-250-20-75-4", -52301.0, 195, 270, True, "HiGHS 9.8 s, 5043 nodes; dense 0.17M",
     "1fa6126895fa59d1c613787ad19bb10e6080a39207224b018a43b95440dca4ef"),
    ("milo-v12-6-r2-40-1", 326481.14282799, 5628, 2688, False, "cut: HiGHS >60 s, dense 110.2M > 4M",
     "a9fdc6d2cbcc4dcd8596f09e45c5ee729ed29a3d0ed93cce4d1f9d683ac8f1b7"),
    ("n5-3", 8104.999999993999, 1062, 2550, False, "cut: dense 6.1M > 4M",
     "48efedf21c4057851a4ef7efc79fbbcbd292b68ee5f6b6a16e904d11ab11bacc"),
    ("neos-1445765", -17783.0, 2147, 20617, False, "cut: dense 58.1M > 4M",
     "195639d1186c8f70b6cf563b0f67f6090b55b3936ade34a32ed40e5dbb67d6fe"),
    ("neos-3004026-krka", 0.0, 12545, 17030, False, "cut: HiGHS >60 s, dense 685.8M > 4M",
     "98e00480fca9f452b1185089c154e14e3af9e79113e1bc0237a2fd444a271516"),
    ("neos-3754480-nidda", 12939.7540104743, 402, 253, False, "cut: HiGHS >60 s",
     "63f322dd9be49acf31f7579d8f04e3bb578bdd7ca7d3525683384190c9a8bf7f"),
    ("neos-5107597-kakapo", 3644.99999999952, 6498, 3114, False, "cut: HiGHS >60 s, dense 146.9M > 4M",
     "53441c15ede747db5fb98670d82633e261054fc36c046b8184e552855e2f5023"),
    ("neos-860300", 3200.999999999998, 850, 1385, True, "HiGHS 9.0 s, 1 nodes; dense 3.34M",
     "d1b2a0a0500834ef6e3a97a6dbfc17c48226c6619056edfedf990e4d83c01504"),
    ("neos-911970", 54.76, 107, 888, True, "HiGHS 5.2 s, 1638 nodes; dense 0.13M",
     "cf9ad35ade0015d140f33ce83a593d76672c09fd8cc30c0539213be7e4a65ebb"),
    ("neos5", 15.0, 63, 63, True, "HiGHS 56.8 s, 139100 nodes; dense 0.02M",
     "2888e9aeb10a4cf16d3b9ad221302fd3de815ce1454b31e660ebd060589c1072"),
    ("ns1830653", 20622.0, 2932, 1629, False, "cut: HiGHS >60 s, dense 30.6M > 4M",
     "a4e3c130f7ff6ebab4d6781a7e95fadc6c4c043e3e65d95b90f3e2e5e46f1017"),
    ("nu25-pr12", 53904.99999999999, 2313, 5868, False, "cut: dense 29.6M > 4M",
     "161c4c89b351d8671ae4e4d7af0fc222e1910e14564419053f22600460ef3236"),
    ("p200x1188c", 15078.0, 1388, 2376, False, "cut: dense 9.1M > 4M",
     "edf0cc0c502379119a65005b5c9a59ad242211598b1741e7fd9f9376b11bc30f"),
    ("pg", -8674.34260712, 125, 2700, True, "HiGHS 3.3 s, 60 nodes; dense 0.38M",
     "b4344d7361aca3d0ff039e5de05d2ad38cfe7c9b9b8300fab4573b35821457b2"),
    ("pg5_34", -14339.35345, 225, 2600, False, "cut: HiGHS >60 s",
     "d169bbd600ddc0be5625d95d52092f63e86d8bdfc78f86662c05429d26519567"),
    ("pk1", 11.0, 45, 86, False, "cut: HiGHS >60 s",
     "6c64530538b254a8ae648f8f8f188847366ae52da560fdec0ad290c642973ba0"),
    ("qap10", 339.9999999983871, 1820, 4150, False, "cut: dense 17.5M > 4M",
     "86c1f147d75688fe154f9dd402032f0c3b8ce164a70e4d0e350fec6c00d20c66"),
    ("ran14x18-disj-8", 3712.0, 447, 504, False, "cut: HiGHS >60 s",
     "2a063d948b69e080e8712c449000bf8f26b05ff4f619552a35b06ca5b6eb38c8"),
    ("rmatr100-p10", 423.0, 7260, 7359, False, "cut: dense 211.5M > 4M",
     "86e4b03e3db2298c29e3ce5863e76b319e95252da2a834d7441d9346dfa4a423"),
    ("sp150x300d", 69.0, 450, 600, True, "HiGHS 0.0 s, 1 nodes; dense 0.88M",
     "47051527397e27d1a62d0462ded35042d5eb5f4c4a3b5e3257c82279b97f1e42"),
    ("supportcase26", 1745.123813, 870, 436, False, "cut: HiGHS >60 s",
     "e5769bd12cc2b5e711165d986ec324d68d336d6314e2b9db4156d5003b9f4ae4"),
    ("timtab1", 764771.99999978, 171, 397, True, "HiGHS 35.3 s, 9692 nodes; dense 0.16M",
     "5867d1a9c1067a849211c1a61bdb8839fb0025255dde09989d0ed333e0674d7a"),
    ("tr12-30", 130595.9999999999, 750, 1080, False, "cut: HiGHS >60 s",
     "94f8cf2fd2fc920ed6ea2452b0755aeec110eb8404716554485c25e7a1e06a3d"),
]
MIPLIB_SLICE = [r[0] for r in MIPLIB if r[4]]
MIPLIB_OBJ = {r[0]: r[1] for r in MIPLIB}


def _download(url: str) -> bytes:
    """urllib first; system curl when this Python has no CA bundle (python.org
    macOS builds without certifi). The sha256 pin checks the bytes either way."""
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            return r.read()
    except urllib.error.URLError as e:
        if not isinstance(e.reason, ssl.SSLError):
            raise
    return subprocess.run(["curl", "-fsSL", "--max-time", "300", url], check=True,
                          capture_output=True).stdout


def fetch_miplib(names) -> int:
    MIPLIB_DEST.mkdir(parents=True, exist_ok=True)
    pins = {r[0]: r[6] for r in MIPLIB}
    failures = []
    for name in names:
        gz, mps = MIPLIB_DEST / f"{name}.mps.gz", MIPLIB_DEST / f"{name}.mps"
        if gz.exists() and hashlib.sha256(gz.read_bytes()).hexdigest() == pins[name] \
                and mps.exists():
            print(f"skip  {name} (already present, hash ok)")
            continue
        try:
            data = _download(MIPLIB_URL.format(name=name))
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"FAIL  {name}: {e}")
            failures.append(name)
            continue
        h = hashlib.sha256(data).hexdigest()
        if h != pins[name]:
            print(f"FAIL  {name}: sha256 {h} != pinned {pins[name]}")
            failures.append(name)
            continue
        gz.write_bytes(data)
        mps.write_bytes(gzip.decompress(data))
        print(f"got   {name}  ({len(data) / 1e6:.2f} MB gz)")
    print(f"miplib: {len(names) - len(failures)}/{len(names)} instances in {MIPLIB_DEST}")
    return 1 if failures else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--set", choices=("netlib", "miplib", "miplib-ref"), default="netlib",
                    help="miplib = the measurable slice, miplib-ref = all 47 pinned instances")
    args = ap.parse_args()
    if args.set == "miplib":
        return fetch_miplib(MIPLIB_SLICE)
    if args.set == "miplib-ref":
        return fetch_miplib([r[0] for r in MIPLIB])
    DEST.mkdir(parents=True, exist_ok=True)
    failures = []
    for name in INSTANCES:
        target = DEST / f"{name}.mps"
        if target.exists() and target.read_text(errors="ignore").lstrip().upper().startswith("NAME"):
            print(f"skip  {name} (already present)")
            continue
        url = URL.format(name=name)
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"FAIL  {name}: {e}")
            failures.append(name)
            continue
        if not data.lstrip().upper().startswith(b"NAME"):
            print(f"FAIL  {name}: body is not an MPS file")
            failures.append(name)
            continue
        target.write_bytes(data)
        print(f"got   {name}  ({len(data) / 1e6:.2f} MB)")
    total = sum(p.stat().st_size for p in DEST.glob("*.mps"))
    print(f"corpus: {len(INSTANCES) - len(failures)}/{len(INSTANCES)} files, "
          f"{total / 1e6:.2f} MB total in {DEST}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
