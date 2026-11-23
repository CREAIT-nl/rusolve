//! AMPL `.nl` reader against independent evidence.
//!
//! * `tests/nl/features.{nl,col,row,lp,mps}`: one Pyomo 6.10 model written
//!   by Pyomo's own .nl, LP and MPS writers (`tests/nl/make_nl.py`); the
//!   `.nl` read must be the same optimisation problem as the LP and MPS
//!   reads (rows compared as single-sided named-term inequalities, so
//!   Pyomo's LP/MPS split of ranged rows and its `ONE_VAR_CONSTANT` column
//!   for the objective constant are normalised away).
//! * `tests/nl/features_mp_{binary,text}.nl`: the same model written by
//!   AMPL's own NL writer (`tests/nl/make_nl_mp.py`), binary and text; Pyomo
//!   writes text only, so this is the real-producer binary file.
//! * Every text file re-encoded as binary in both byte orders by an encoder
//!   in this test (big-endian has no producer on this host); the binary read
//!   must equal the text read exactly.

mod common;

use common::canon;
use rusolve_io::{NlNames, ObjSense, read_model_file, read_nl};

const INF: f64 = f64::INFINITY;

fn path(f: &str) -> String {
    format!("{}/tests/{f}", env!("CARGO_MANIFEST_DIR"))
}

fn names(stem: &str) -> NlNames {
    let lines = |ext: &str| {
        std::fs::read_to_string(path(&format!("{stem}.{ext}")))
            .ok()
            .map(|s| s.lines().map(str::to_string).collect())
    };
    NlNames {
        col: lines("col"),
        row: lines("row"),
    }
}

/// Re-encode a linear text (`g`) .nl file as binary (`b`), written from the
/// segment layout of Gay, "Writing .nl Files" (binary: "int (4 byte) and
/// double (8 byte) values ... a string consists of an int length followed
/// by that many characters"; segment keys and r/b line types one byte) and
/// independent of the reader: it tokenises the text itself. The header
/// stays text with `b` first and `arith` (line 6, third field) set to 1
/// (little-endian) or 2 (big-endian). Integer constants are written as `s`
/// (2-byte) in C segments and `l` (4-byte) in O segments so every numeric
/// node kind is exercised.
fn encode_binary(text: &str, big_endian: bool) -> Vec<u8> {
    let lines: Vec<&str> = text.split('\n').collect();
    let mut head: Vec<String> = lines[..10].iter().map(|l| l.to_string()).collect();
    assert!(head[0].starts_with('g'));
    head[0].replace_range(0..1, "b");
    let mut l6: Vec<String> = head[5]
        .split('#')
        .next()
        .unwrap()
        .split_whitespace()
        .map(str::to_string)
        .collect();
    l6.resize(4.max(l6.len()), "0".into());
    l6[2] = if big_endian { "2" } else { "1" }.into();
    head[5] = format!(" {}", l6.join(" "));
    let ints = |l: &str| -> Vec<usize> {
        l.split('#')
            .next()
            .unwrap()
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect()
    };
    let (n_var, n_con) = (ints(&head[1])[0], ints(&head[1])[1]);
    let toks: Vec<&str> = lines[10..]
        .iter()
        .flat_map(|l| l.split('#').next().unwrap().split_whitespace())
        .collect();
    let mut out = Vec::new();
    let int = |o: &mut Vec<u8>, v: &str| {
        let v: i32 = v.parse().unwrap();
        o.extend(if big_endian {
            v.to_be_bytes()
        } else {
            v.to_le_bytes()
        });
    };
    let real = |o: &mut Vec<u8>, v: &str| {
        let v: f64 = v.parse().unwrap();
        o.extend(if big_endian {
            v.to_be_bytes()
        } else {
            v.to_le_bytes()
        });
    };
    let mut it = toks.into_iter();
    while let Some(t) = it.next() {
        let (key, rest) = t.split_at(1);
        out.push(key.as_bytes()[0]);
        match key {
            "C" | "O" => {
                int(&mut out, rest);
                if key == "O" {
                    int(&mut out, it.next().unwrap());
                }
                let e = it.next().unwrap();
                assert!(e.starts_with('n'), "linear models only: {e}");
                let v: f64 = e[1..].parse().unwrap();
                if v.fract() == 0.0 && v.abs() < 32768.0 && key == "C" {
                    out.push(b's');
                    let v = v as i16;
                    out.extend(if big_endian {
                        v.to_be_bytes()
                    } else {
                        v.to_le_bytes()
                    });
                } else if v.fract() == 0.0 && v.abs() < 2e9 && key == "O" {
                    out.push(b'l');
                    int(&mut out, &(v as i32).to_string());
                } else {
                    out.push(b'n');
                    real(&mut out, &e[1..]);
                }
            }
            "x" | "d" => {
                int(&mut out, rest);
                for _ in 0..rest.parse::<usize>().unwrap() {
                    int(&mut out, it.next().unwrap());
                    real(&mut out, it.next().unwrap());
                }
            }
            "r" | "b" => {
                for _ in 0..if key == "r" { n_con } else { n_var } {
                    let ty = it.next().unwrap();
                    out.push(ty.as_bytes()[0]);
                    for _ in 0..[2, 1, 1, 0, 1][ty.parse::<usize>().unwrap()] {
                        real(&mut out, it.next().unwrap());
                    }
                }
            }
            "k" | "K" => {
                int(&mut out, rest);
                for _ in 0..rest.parse::<usize>().unwrap() {
                    int(&mut out, it.next().unwrap());
                }
            }
            "J" | "G" => {
                int(&mut out, rest);
                let m = it.next().unwrap();
                int(&mut out, m);
                for _ in 0..m.parse::<usize>().unwrap() {
                    int(&mut out, it.next().unwrap());
                    real(&mut out, it.next().unwrap());
                }
            }
            "S" => {
                int(&mut out, rest);
                let kind: i32 = rest.parse().unwrap();
                let m = it.next().unwrap();
                int(&mut out, m);
                let name = it.next().unwrap();
                int(&mut out, &name.len().to_string());
                out.extend(name.as_bytes());
                for _ in 0..m.parse::<usize>().unwrap() {
                    int(&mut out, it.next().unwrap());
                    if kind & 4 != 0 {
                        real(&mut out, it.next().unwrap());
                    } else {
                        int(&mut out, it.next().unwrap());
                    }
                }
            }
            _ => panic!("segment {t} not in the linear subset"),
        }
    }
    let mut bytes = (head.join("\n") + "\n").into_bytes();
    bytes.extend(out);
    bytes
}

#[test]
fn binary_nl_reads_exactly_like_text_in_both_byte_orders() {
    for stem in ["dialects/pyomo", "nl/features", "nl/features_mp_text"] {
        let n = names(stem);
        let text = std::fs::read_to_string(path(&format!("{stem}.nl"))).unwrap();
        let want = read_nl(text.as_bytes(), &n).unwrap();
        for big_endian in [false, true] {
            let b = encode_binary(&text, big_endian);
            let got = read_nl(&b, &n).unwrap_or_else(|e| panic!("{stem} BE={big_endian}: {e}"));
            assert_eq!(got, want, "{stem}, big-endian {big_endian}");
            // The byte order really comes from `arith`: the same bytes
            // labelled with the other order do not read to this model.
            let (from, to) = if big_endian {
                (" 2 ", " 1 ")
            } else {
                (" 1 ", " 2 ")
            };
            let hl: Vec<usize> = b
                .iter()
                .enumerate()
                .filter(|x| *x.1 == b'\n')
                .map(|x| x.0)
                .take(6)
                .collect();
            let line6 = String::from_utf8(b[hl[4] + 1..hl[5]].to_vec()).unwrap();
            let mut lied = b[..hl[4] + 1].to_vec();
            lied.extend(line6.replacen(from, to, 1).into_bytes());
            lied.extend(&b[hl[5]..]);
            assert_ne!(
                read_nl(&lied, &n).ok(),
                Some(want.clone()),
                "{stem}: arith ignored"
            );
        }
    }
}

/// AMPL's own NL writer (MP NLWriter2 via `nlwpy`, `tests/nl/make_nl_mp.py`)
/// wrote the features model as binary and as text on a little-endian host:
/// both read to the same problem as Pyomo's .nl, LP and MPS files of it.
#[test]
fn binary_nl_from_ampls_writer() {
    let bin = std::fs::read(path("nl/features_mp_binary.nl")).unwrap();
    assert_eq!(bin[0], b'b');
    let m = read_model_file(&path("nl/features_mp_binary.nl")).unwrap();
    let t = read_model_file(&path("nl/features_mp_text.nl")).unwrap();
    assert_eq!(m.matrix, t.matrix);
    assert_eq!((&m.col_names, &m.row_names), (&t.col_names, &t.row_names));
    assert_eq!(m.obj_name.as_deref(), Some("obj"));
    let c = canon(&m);
    assert_eq!(canon(&t), c);
    for other in ["nl/features.nl", "nl/features.lp", "nl/features.mps"] {
        assert_eq!(canon(&read_model_file(&path(other)).unwrap()), c, "{other}");
    }
}

#[test]
fn nl_is_the_same_problem_as_pyomos_lp_and_mps_exports() {
    let nl = read_model_file(&path("nl/features.nl")).unwrap();
    assert_eq!(nl.row_names, ["eq", "rng", "le", "ge", "cst", "dup", "box"]);
    assert_eq!(nl.obj_name.as_deref(), Some("obj"));
    assert_eq!((nl.obj_sense, nl.obj_offset), (ObjSense::Minimize, -3.5));
    let c = canon(&nl);
    for other in ["nl/features.lp", "nl/features.mps"] {
        let m = read_model_file(&path(other)).unwrap();
        assert_eq!(canon(&m), c, "{other}");
    }
    // Spot-check the feature mix the model was built for.
    let col = |n: &str| nl.bounds[nl.col_names.iter().position(|c| c == n).unwrap()];
    assert_eq!(
        (col("i").lower, col("i").upper, col("i").integer),
        (-5.0, 5.0, true)
    );
    assert_eq!(
        (col("b").lower, col("b").upper, col("b").integer),
        (0.0, 1.0, true)
    );
    assert_eq!((col("f").lower, col("f").upper), (-INF, INF));
    assert_eq!((col("u").lower, col("u").upper), (-INF, 4.0));
    assert_eq!((col("w").lower, col("w").upper), (1.5, 1.5));
}
