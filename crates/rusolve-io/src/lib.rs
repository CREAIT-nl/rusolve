//! File-format readers/writers for rusolve: MPS, CPLEX LP, AMPL `.nl` and
//! MathOptFormat JSON (readers), MPS (writer).
//!
//! Format-agnostic file entry points ([`read_model_file`] and friends)
//! decompress transparently and pick the reader by extension, falling back
//! to a content sniff for opaque paths:
//!
//! * Compression is detected by the leading magic bytes, never by the name,
//!   so a mislabelled or extension-less file still decodes: gzip `1f 8b`
//!   (RFC 1952 §2.3.1; concatenated members decoded, as `gzip -d` does),
//!   bzip2 `BZh` (the bzip2 1.0 stream header), xz `fd 37 7a 58 5a 00`
//!   (The .xz File Format 1.2.1, §2.1.1.1), zstd `28 b5 2f fd`
//!   (RFC 8878 §3.1.1; several frames and skippable frames handled). All
//!   decoders are pure Rust: `flate2` (zlib-rs backend), `bzip2` (libbz2-rs-sys),
//!   `lzma-rs`, `ruzstd`.
//! * The inner format comes from the extension with the compression suffix
//!   stripped (`.mps`, `.lp`, `.nl`, `.mof.json`), else from the
//!   content: see [`sniff_format`].

use std::io::{BufRead as _, Read as _};

pub mod conic;
pub mod cp;
pub mod cp_lin;
pub mod diag;
pub mod gc;
pub mod gc_lin;
pub mod lp;
pub mod lp_write;
pub mod model;
pub mod mof;
pub mod mps;
pub mod nl;

pub use cp::{
    AllDifferent, CpError, CpErrorKind, CpModel, CpTerm, CpTol, CpViolation, Cumulative, Element,
    Interval, NoOverlap, Table,
};
pub use gc::{
    FuncCon, FuncKind, GcModel, GcTol, GcViolation, Indicator, Pwl, SemiVar, SosKind, SosSet,
};
pub use lp::{LpError, read_lp, read_lp_conic};
pub use lp_write::{
    LpNames, LpWriteReport, cplex_read, lp_names, write_lp, write_lp_report, write_model_file,
    write_model_file_report,
};
pub use model::{Model, ObjSense, RowSense, VariableBound};
pub use mof::{MofError, read_mof};
pub use mps::{
    MpsError, MpsFileOrder, read_mps, read_mps_conic, read_mps_conic_order, read_mps_file_order,
    write_mps,
};
pub use nl::{NlError, NlNames, read_nl};

/// A file read failure: I/O (decompression included), or one of the format
/// parse errors.
#[derive(Debug)]
pub enum ReadError {
    Io(std::io::Error),
    Mps(MpsError),
    Lp(LpError),
    Nl(NlError),
    Mof(MofError),
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Io(e) => write!(f, "read error: {e}"),
            ReadError::Mps(e) => e.fmt(f),
            ReadError::Lp(e) => e.fmt(f),
            ReadError::Nl(e) => e.fmt(f),
            ReadError::Mof(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for ReadError {}

/// Compression wrapper of a model file, by magic bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Bzip2,
    Xz,
    Zstd,
}

impl Compression {
    /// The wrapper the leading bytes announce (see the module doc for the
    /// magic numbers and their sources).
    pub fn sniff(head: &[u8]) -> Compression {
        if head.starts_with(&[0x1f, 0x8b]) {
            Compression::Gzip
        } else if head.len() >= 4 && head.starts_with(b"BZh") && (b'1'..=b'9').contains(&head[3]) {
            Compression::Bzip2
        } else if head.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]) {
            Compression::Xz
        } else if head.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
            Compression::Zstd
        } else {
            Compression::None
        }
    }
}

/// File-name suffixes of the compression wrappers, stripped before the
/// inner format is read from the extension.
const COMPRESSION_SUFFIXES: [&str; 6] = [".gz", ".bz2", ".xz", ".zst", ".zstd", ".gzip"];

fn io_err(what: &str, e: impl std::fmt::Display) -> ReadError {
    ReadError::Io(std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("{what}: {e}"),
    ))
}

/// Whole file, decompressed according to its magic bytes.
pub fn read_file_bytes(path: &str) -> Result<Vec<u8>, ReadError> {
    let file = std::fs::File::open(path).map_err(ReadError::Io)?;
    let file_len = file.metadata().map(|m| m.len()).unwrap_or(0);
    let mut r = std::io::BufReader::with_capacity(1 << 16, file);
    let comp = Compression::sniff(r.fill_buf().map_err(ReadError::Io)?);
    let mut out = Vec::new();
    match comp {
        Compression::None => {
            out.reserve_exact(file_len as usize);
            r.read_to_end(&mut out).map_err(ReadError::Io)?;
        }
        Compression::Gzip => {
            // Pre-size from the compressed length (text models compress
            // ~5-10x); `read_to_end` grows from there.
            out.reserve((file_len as usize).saturating_mul(6));
            flate2::bufread::MultiGzDecoder::new(r)
                .read_to_end(&mut out)
                .map_err(|e| io_err("gzip", e))?;
        }
        Compression::Bzip2 => {
            out.reserve((file_len as usize).saturating_mul(6));
            bzip2::bufread::MultiBzDecoder::new(r)
                .read_to_end(&mut out)
                .map_err(|e| io_err("bzip2", e))?;
        }
        Compression::Xz => {
            out.reserve((file_len as usize).saturating_mul(6));
            lzma_rs::xz_decompress(&mut r, &mut out).map_err(|e| io_err("xz", e))?;
        }
        Compression::Zstd => {
            let mut input = Vec::with_capacity(file_len as usize);
            r.read_to_end(&mut input).map_err(ReadError::Io)?;
            out.reserve(input.len().saturating_mul(6));
            zstd_decode_all(&input, &mut out)?;
        }
    }
    Ok(out)
}

/// Every zstd frame in `input`, skippable frames skipped (RFC 8878 §3.1.2).
fn zstd_decode_all(mut input: &[u8], out: &mut Vec<u8>) -> Result<(), ReadError> {
    use ruzstd::decoding::errors::{FrameDecoderError, ReadFrameHeaderError};
    use ruzstd::decoding::{BlockDecodingStrategy, FrameDecoder};
    let mut dec = FrameDecoder::new();
    while !input.is_empty() {
        match dec.init(&mut input) {
            Ok(()) => {}
            Err(FrameDecoderError::ReadFrameHeaderError(ReadFrameHeaderError::SkipFrame {
                length,
                ..
            })) => {
                input = input
                    .get(length as usize..)
                    .ok_or_else(|| io_err("zstd", "truncated skippable frame"))?;
                continue;
            }
            Err(e) => return Err(io_err("zstd", e)),
        }
        loop {
            dec.decode_blocks(&mut input, BlockDecodingStrategy::UptoBytes(1 << 20))
                .map_err(|e| io_err("zstd", e))?;
            dec.collect_to_writer(&mut *out)
                .map_err(|e| io_err("zstd", e))?;
            if dec.is_finished() {
                break;
            }
        }
        dec.collect_to_writer(&mut *out)
            .map_err(|e| io_err("zstd", e))?;
    }
    Ok(())
}

fn into_text(bytes: Vec<u8>) -> Result<String, ReadError> {
    String::from_utf8(bytes).map_err(|e| {
        ReadError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("stream did not contain valid UTF-8: {e}"),
        ))
    })
}

/// Whole file as text, decompressed according to its magic bytes.
fn read_text(path: &str) -> Result<String, ReadError> {
    into_text(read_file_bytes(path)?)
}

/// Read an MPS file (compressed or not).
pub fn read_mps_file(path: &str) -> Result<Model, ReadError> {
    read_mps(&read_text(path)?).map_err(ReadError::Mps)
}

/// [`read_mps_file`], also returning the COLUMNS entries in file order.
pub fn read_mps_file_file_order(path: &str) -> Result<(Model, MpsFileOrder), ReadError> {
    read_mps_file_order(&read_text(path)?).map_err(ReadError::Mps)
}

/// Read a CPLEX LP-format file (compressed or not).
pub fn read_lp_file(path: &str) -> Result<Model, ReadError> {
    read_lp(&read_text(path)?).map_err(ReadError::Lp)
}

/// A model file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Free (and fixed) MPS.
    Mps,
    /// CPLEX LP format.
    Lp,
    /// AMPL `.nl`, text (`g`) or binary (`b`).
    Nl,
    /// MathOptFormat JSON (`.mof.json`).
    Mof,
}

/// The format a path names by extension, compression suffix stripped
/// (case-insensitive); `None` for anything else.
pub fn format_from_extension(path: &str) -> Option<Format> {
    let lower = path.to_ascii_lowercase();
    let mut stem = lower.as_str();
    if let Some(s) = COMPRESSION_SUFFIXES
        .iter()
        .find_map(|x| stem.strip_suffix(x))
    {
        stem = s;
    }
    if stem.ends_with(".mps") {
        Some(Format::Mps)
    } else if stem.ends_with(".lp") {
        Some(Format::Lp)
    } else if stem.ends_with(".nl") {
        Some(Format::Nl)
    } else if stem.ends_with(".json") {
        Some(Format::Mof)
    } else {
        None
    }
}

/// The format of decompressed content, for paths whose extension names
/// none:
/// * `.nl`: the first line is `g` or `b` followed by digits (the header of
///   Gay, "Writing .nl Files", §2) and the second line is all integers;
/// * MathOptFormat: the first non-blank byte is `{` (a JSON object);
/// * LP: the first non-blank line is a `\` comment (LP-format comments;
///   MPS comments start with `*`) or starts with an objective keyword
///   (`min`/`max` and their long forms, the LP file's first section);
/// * MPS otherwise.
pub fn sniff_format(bytes: &[u8]) -> Format {
    let mut lines = bytes.split(|&b| b == b'\n');
    let first = lines.next().unwrap_or_default();
    if matches!(first.first(), Some(b'g' | b'b'))
        && first.len() > 1
        && first[1].is_ascii_digit()
        && lines.next().is_some_and(|l| {
            let l = l.split(|&b| b == b'#').next().unwrap_or_default();
            let mut toks = l
                .split(|b| b.is_ascii_whitespace())
                .filter(|t| !t.is_empty());
            toks.clone().count() >= 3 && toks.all(|t| t.iter().all(u8::is_ascii_digit))
        })
    {
        return Format::Nl;
    }
    if bytes.iter().find(|b| !b.is_ascii_whitespace()) == Some(&b'{') {
        return Format::Mof;
    }
    let start = bytes.iter().position(|b| !b.is_ascii_whitespace());
    let Some(start) = start else {
        return Format::Mps;
    };
    let line = bytes[start..]
        .split(|&b| b == b'\n')
        .next()
        .unwrap_or_default();
    if line.first() == Some(&b'\\') {
        return Format::Lp;
    }
    let word: Vec<u8> = line
        .iter()
        .take_while(|b| b.is_ascii_alphabetic())
        .map(u8::to_ascii_lowercase)
        .collect();
    if matches!(
        word.as_slice(),
        b"min" | b"max" | b"minimize" | b"maximize" | b"minimum" | b"maximum"
    ) {
        return Format::Lp;
    }
    Format::Mps
}

/// The `.nl` sidecar name files next to `path` (`stem.col`, `stem.row`,
/// AMPL's `option auxfiles rc`), when present; a compression suffix and the
/// `.nl` extension are stripped to find the stem.
fn nl_sidecar_names(path: &str) -> Result<NlNames, ReadError> {
    let lower = path.to_ascii_lowercase();
    let mut cut = path.len();
    if let Some(s) = COMPRESSION_SUFFIXES
        .iter()
        .find_map(|x| lower.strip_suffix(x))
    {
        cut = s.len();
    }
    if lower[..cut].ends_with(".nl") {
        cut -= 3;
    }
    let stem = &path[..cut];
    let read = |ext: &str| -> Result<Option<Vec<String>>, ReadError> {
        let p = format!("{stem}.{ext}");
        match std::fs::read_to_string(&p) {
            Ok(s) => Ok(Some(s.lines().map(|l| l.trim_end().to_string()).collect())),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(ReadError::Io(std::io::Error::new(
                e.kind(),
                format!("{p}: {e}"),
            ))),
        }
    };
    Ok(NlNames {
        col: read("col")?,
        row: read("row")?,
    })
}

fn read_format(fmt: Format, bytes: Vec<u8>, path: &str) -> Result<Model, ReadError> {
    match fmt {
        Format::Mps => read_mps(&into_text(bytes)?).map_err(ReadError::Mps),
        Format::Lp => read_lp(&into_text(bytes)?).map_err(ReadError::Lp),
        Format::Nl => read_nl(&bytes, &nl_sidecar_names(path)?).map_err(ReadError::Nl),
        Format::Mof => read_mof(&bytes).map_err(ReadError::Mof),
    }
}

/// Read a model file, auto-detecting compression (magic bytes) and format
/// (extension, else [`sniff_format`]).
pub fn read_model_file(path: &str) -> Result<Model, ReadError> {
    let bytes = read_file_bytes(path)?;
    let fmt = format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes));
    read_format(fmt, bytes, path)
}

/// [`read_model_file`], but an MPS or LP file with quadratic rows or cones
/// is kept (in [`Model::conic`]) instead of rejected. Only the conic route
/// calls it; every other format reads as [`read_model_file`] does.
pub fn read_model_file_conic(path: &str) -> Result<Model, ReadError> {
    let bytes = read_file_bytes(path)?;
    let fmt = format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes));
    match fmt {
        Format::Mps => read_mps_conic(&into_text(bytes)?).map_err(ReadError::Mps),
        Format::Lp => read_lp_conic(&into_text(bytes)?).map_err(ReadError::Lp),
        other => read_format(other, bytes, path),
    }
}

/// Whether the MPS or LP file at `path` may carry quadratic rows or cones: an MPS line that
/// starts a `QCMATRIX` or `CSECTION` section, or any `[` in an LP file (which also opens an
/// objective quadratic, so `true` is "look closer", never "conic"). Any other format, and any
/// read or decode failure, is `false`: the legacy path then reads the file and reports its own
/// error. One pass over the (decompressed) bytes, no parse. The CLI uses it to pick the conic
/// route before a legacy reader can print a rejection.
pub fn may_have_conic_content(path: &str) -> bool {
    let Ok(bytes) = read_file_bytes(path) else {
        return false;
    };
    match format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes)) {
        Format::Mps => bytes
            .split(|&b| b == b'\n')
            .any(|l| l.starts_with(b"QCMATRIX") || l.starts_with(b"CSECTION")),
        Format::Lp => bytes.contains(&b'['),
        _ => false,
    }
}

/// Whether the file at `path` may carry anything beyond a plain LP or MIP: a quadratic term, a
/// quadratic row, a cone, semi-continuous columns or a general constraint. A liberal byte scan,
/// never a parse: `true` means "read it and classify" (a plain model can still say `true`, e.g. an
/// MPS column named `SC`), `false` means plain, so the router leaves such a file to the LP/MIP
/// path without reading it twice. `.nl` and MathOptFormat files are always `true`; a read or
/// decode failure is `false`, so the legacy path reports its own error.
pub fn may_have_special_content(path: &str) -> bool {
    let Ok(bytes) = read_file_bytes(path) else {
        return false;
    };
    match format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes)) {
        Format::Mps => bytes.split(|&b| b == b'\n').any(|l| {
            let first = l.first().copied().unwrap_or(b' ');
            if first == b'*' || first == b'\r' {
                return false;
            }
            let word = l
                .split(|b| b.is_ascii_whitespace())
                .find(|w| !w.is_empty())
                .unwrap_or_default();
            if first.is_ascii_whitespace() {
                // a bound type of semi-continuous / semi-integer columns
                word.eq_ignore_ascii_case(b"SC") || word.eq_ignore_ascii_case(b"SI")
            } else {
                ![
                    &b"NAME"[..],
                    b"ROWS",
                    b"COLUMNS",
                    b"RHS",
                    b"RANGES",
                    b"BOUNDS",
                    b"ENDATA",
                    b"OBJSENSE",
                    b"OBJNAME",
                ]
                .iter()
                .any(|k| word.eq_ignore_ascii_case(k))
            }
        }),
        Format::Lp => {
            let lower = bytes.to_ascii_lowercase();
            let has = |needle: &[u8]| lower.windows(needle.len()).any(|w| w == needle);
            has(b"[") || has(b"->") || has(b"semi") || has(b"sos") || has(b"constraints")
        }
        _ => true,
    }
}

/// [`read_model_file`], also returning the constraint entries in file
/// order: the MPS reader's COLUMNS order for MPS input, the CSC column
/// order (and costs in column order) for every other format, which has no
/// COLUMNS section to keep.
pub fn read_model_file_order(path: &str) -> Result<(Model, MpsFileOrder), ReadError> {
    let bytes = read_file_bytes(path)?;
    let fmt = format_from_extension(path).unwrap_or_else(|| sniff_format(&bytes));
    if fmt == Format::Mps {
        return read_mps_file_order(&into_text(bytes)?).map_err(ReadError::Mps);
    }
    let m = read_format(fmt, bytes, path)?;
    let order = MpsFileOrder::of_csc(&m);
    Ok((m, order))
}

#[cfg(test)]
mod tests {
    use super::*;

    const MPS: &str = "\
NAME GZ
ROWS
  N  OF
  L  R1
COLUMNS
  X  OF  1  R1  1
RHS
  rhs R1 3
ENDATA
";

    const LP: &str = "\
\\ LP format for model gz
Minimize
 obj: x
Subject To
 c: x >= 2
End
";

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("rusolve-io-{}-{}", tag, std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn gzip_bytes(text: &str) -> Vec<u8> {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write as _;
        let mut e = GzEncoder::new(Vec::new(), Compression::default());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn gz_and_plain_paths_parse_identically() {
        let d = tmpdir("gz");
        std::fs::write(d.join("m.mps"), MPS).unwrap();
        std::fs::write(d.join("m.mps.gz"), gzip_bytes(MPS)).unwrap();
        std::fs::write(d.join("m.lp"), LP).unwrap();
        std::fs::write(d.join("m.lp.gz"), gzip_bytes(LP)).unwrap();

        let m_plain = read_mps_file(&d.join("m.mps").display().to_string()).unwrap();
        let m_gz = read_mps_file(&d.join("m.mps.gz").display().to_string()).unwrap();
        assert_eq!(m_gz, m_plain);
        let l_plain = read_lp_file(&d.join("m.lp").display().to_string()).unwrap();
        let l_gz = read_lp_file(&d.join("m.lp.gz").display().to_string()).unwrap();
        assert_eq!(l_gz, l_plain);
    }

    #[test]
    fn read_model_file_dispatches_by_extension_and_sniff() {
        let d = tmpdir("dispatch");
        std::fs::write(d.join("m.mps.gz"), gzip_bytes(MPS)).unwrap();
        std::fs::write(d.join("m.lp"), LP).unwrap();
        // opaque extensions: sniff the first non-empty line
        std::fs::write(d.join("opaque.gz"), gzip_bytes(LP)).unwrap();
        std::fs::write(d.join("opaque2.gz"), gzip_bytes(MPS)).unwrap();

        let m = read_model_file(&d.join("m.mps.gz").display().to_string()).unwrap();
        assert_eq!(m.name, "GZ");
        let l = read_model_file(&d.join("m.lp").display().to_string()).unwrap();
        assert_eq!(l.name, "gz");
        let l2 = read_model_file(&d.join("opaque.gz").display().to_string()).unwrap();
        assert_eq!(l2.name, "gz", "LP marker sniffed through bare .gz");
        let m2 = read_model_file(&d.join("opaque2.gz").display().to_string()).unwrap();
        assert_eq!(m2.name, "GZ", "no marker ⇒ MPS");
    }
}
