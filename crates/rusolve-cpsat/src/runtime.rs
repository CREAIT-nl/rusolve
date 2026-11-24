//! Run-time loading of the OR-Tools shared library.
//!
//! Nothing of OR-Tools is needed to build rusolve: the library is opened with the system loader
//! (`dlopen` / `LoadLibrary`, through `libloading`, ISC) the first time CP-SAT is selected, and
//! the one function the engine calls, `SolveCpModelWithParameters` of
//! `ortools/sat/c_api/cp_solver_c.h`, is looked up by name. The library is never unloaded.
//!
//! ## Search order (THE single source of truth: [`candidates`] builds it, the docs link here)
//!
//! The first entry that loads wins.
//!
//! 1. `RUSOLVE_ORTOOLS_LIB`: the path of the library file. When it is set, nothing else is
//!    tried: a variable that names no loadable library is an error, not a reason to use another.
//! 2. `<dir of the executable>/../lib/rusolve/ortools/<libname>`: where `install.sh --with-cpsat`
//!    puts the files of the `ortools` wheel's `.libs/`.
//! 3. The `ortools` Python package: the directory the Python layer reports through [`set_hint`]
//!    (the package found next to `rusolve` in the same environment), then
//!    `<dir of the executable>/../../ortools/.libs/<libname>` (the CLI of the `rusolve` wheel is
//!    `site-packages/rusolve/_bin/rusolve-cli`, `ortools` is `site-packages/ortools`).
//! 4. The system loader's own search for the bare `<libname>` (`LD_LIBRARY_PATH`, `ld.so.cache`,
//!    `DYLD_LIBRARY_PATH`, `PATH`).
//!
//! `<libname>` is `libortools.so.9` (Linux), `libortools.9.dylib` (macOS) or `ortools.dll`
//! (Windows): the file names the OR-Tools 9 wheels carry. The request and the response are
//! protobuf messages of `cp_model.proto` (vendored, 9.15), so any 9.x library speaks them.
//!
//! Version: the library does not export its version through this API, so it is read from the
//! files around it: `VERSION` next to the library (written by `install.sh`), else the
//! `ortools-<version>.dist-info` of the wheel it came from, else the version in the library's
//! real file name (`libortools.so.9.15.6755`), else `unknown`.

use libloading::Library;
use std::ffi::{c_int, c_void};
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// The environment variable of step 1.
pub const ENV_LIB: &str = "RUSOLVE_ORTOOLS_LIB";

/// The file name of the library on this platform (see the module docs).
pub const LIB_NAME: &str = if cfg!(target_os = "windows") {
    "ortools.dll"
} else if cfg!(target_os = "macos") {
    "libortools.9.dylib"
} else {
    "libortools.so.9"
};

/// `SolveCpModelWithParameters` as declared in `ortools/sat/c_api/cp_solver_c.h` (OR-Tools 9.15).
pub type SolveFn = unsafe extern "C" fn(
    creq: *const c_void,
    creq_len: c_int,
    cparams: *const c_void,
    cparams_len: c_int,
    cres: *mut *mut c_void,
    cres_len: *mut c_int,
);

const SYMBOL: &[u8] = b"SolveCpModelWithParameters\0";

/// One place of the search order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Candidate {
    /// Which step: `env`, `install`, `ortools-package`, `ortools-wheel` or `system`.
    pub origin: &'static str,
    /// A file path, or the bare library name for the system loader.
    pub path: PathBuf,
}

impl fmt::Display for Candidate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.path.display(), self.origin)
    }
}

static HINT: OnceLock<PathBuf> = OnceLock::new();

/// The Python layer calls this once with the directory of the `ortools` package it found
/// (`ortools/`, holding `.libs/`), because inside Python the executable is the interpreter and
/// the relative steps do not apply. A second call is ignored.
pub fn set_hint(ortools_package_dir: PathBuf) {
    let _ = HINT.set(ortools_package_dir);
}

/// `p` without `.` and `..` components, decided lexically (the display of a search path).
fn clean(p: &Path) -> PathBuf {
    use std::path::Component;
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            Component::ParentDir if out.file_name().is_some() => {
                out.pop();
            }
            Component::CurDir => {}
            c => out.push(c.as_os_str()),
        }
    }
    out
}

/// The search order for an executable at `exe` and the value of `RUSOLVE_ORTOOLS_LIB`
/// (the pure part of the search, tested without touching the process).
pub fn candidates_for(
    env_lib: Option<PathBuf>,
    exe: Option<&Path>,
    hint: Option<&Path>,
) -> Vec<Candidate> {
    if let Some(p) = env_lib {
        return vec![Candidate {
            origin: "env",
            path: p,
        }];
    }
    let mut v = Vec::new();
    let exe_dir = exe.and_then(Path::parent);
    if let Some(d) = exe_dir {
        v.push(Candidate {
            origin: "install",
            path: clean(&d.join("../lib/rusolve/ortools").join(LIB_NAME)),
        });
    }
    if let Some(h) = hint {
        v.push(Candidate {
            origin: "ortools-package",
            path: h.join(".libs").join(LIB_NAME),
        });
    }
    if let Some(d) = exe_dir {
        v.push(Candidate {
            origin: "ortools-wheel",
            path: clean(&d.join("../../ortools/.libs").join(LIB_NAME)),
        });
    }
    v.push(Candidate {
        origin: "system",
        path: PathBuf::from(LIB_NAME),
    });
    v
}

/// [`candidates_for`] for this process.
pub fn candidates() -> Vec<Candidate> {
    let env_lib = std::env::var_os(ENV_LIB)
        .filter(|v| !v.is_empty())
        .map(PathBuf::from);
    candidates_for(
        env_lib,
        std::env::current_exe().ok().as_deref(),
        HINT.get().map(PathBuf::as_path),
    )
}

/// Why the library is not available.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LoadError {
    /// Every place tried, in order.
    pub searched: Vec<String>,
    /// What the loader said about each file that exists but did not load, or lacks the symbol.
    pub problems: Vec<String>,
}

impl fmt::Display for LoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "the OR-Tools library ({LIB_NAME}) was not found. Searched, in order: {}.",
            self.searched.join("; ")
        )?;
        if !self.problems.is_empty() {
            write!(f, " Problems: {}.", self.problems.join("; "))?;
        }
        write!(
            f,
            " Install it with `install.sh --with-cpsat`, or `pip install 'rusolve[cpsat]'`, \
             or point {ENV_LIB} at the library file (docs/design/cpsat-engine.md)"
        )
    }
}

impl std::error::Error for LoadError {}

/// The loaded library.
pub struct Loaded {
    /// Where it came from: a file path, or the bare name when only the system loader knows.
    pub path: PathBuf,
    /// Which step of the search found it.
    pub origin: &'static str,
    /// `9.15.6755`, or `unknown`.
    pub version: String,
    pub(crate) solve: SolveFn,
    // keeps the library mapped: `solve` points into it. The static below is never dropped.
    _lib: Library,
}

static LOADED: OnceLock<Result<Loaded, LoadError>> = OnceLock::new();

/// Opens the library on first use (once per process) and returns it, or why it is missing.
pub fn library() -> Result<&'static Loaded, &'static LoadError> {
    LOADED.get_or_init(load).as_ref()
}

fn open(path: &Path) -> Result<Library, libloading::Error> {
    // SAFETY: opening a shared library runs its initializers. The file is one the user or the
    // installer put in the documented places (or named in RUSOLVE_ORTOOLS_LIB); there is no
    // safe way to open a library at all, and the loader is the same one the program started with.
    #[cfg(windows)]
    unsafe {
        // LOAD_WITH_ALTERED_SEARCH_PATH: the DLLs next to ortools.dll are found there
        libloading::os::windows::Library::load_with_flags(path, 0x0000_0008).map(Library::from)
    }
    #[cfg(not(windows))]
    unsafe {
        Library::new(path)
    }
}

fn load() -> Result<Loaded, LoadError> {
    let cands = candidates();
    let mut problems = Vec::new();
    for c in &cands {
        // a path with a directory part must exist; the bare name is the system loader's to find
        let bare = c.origin == "system";
        if !bare && !c.path.is_file() {
            continue;
        }
        let lib = match open(&c.path) {
            Ok(l) => l,
            Err(e) => {
                problems.push(format!("{}: {e}", c.path.display()));
                if c.origin == "env" {
                    break;
                }
                continue;
            }
        };
        // SAFETY: the symbol has the signature of `SolveFn` (cp_solver_c.h, OR-Tools 9.x); the
        // pointer is copied out and `lib` is stored next to it in `Loaded`, so it stays valid.
        let solve: SolveFn = match unsafe { lib.get::<SolveFn>(SYMBOL) } {
            Ok(s) => *s,
            Err(e) => {
                problems.push(format!(
                    "{} has no SolveCpModelWithParameters ({e}): not an OR-Tools 9 library",
                    c.path.display()
                ));
                if c.origin == "env" {
                    break;
                }
                continue;
            }
        };
        let path = if bare {
            loaded_path().unwrap_or_else(|| c.path.clone())
        } else {
            c.path.clone()
        };
        return Ok(Loaded {
            version: version_near(&path),
            path,
            origin: c.origin,
            solve,
            _lib: lib,
        });
    }
    Err(LoadError {
        searched: cands.iter().map(Candidate::to_string).collect(),
        problems,
    })
}

/// Where the system loader found `libortools`, on platforms where the process map says so.
fn loaded_path() -> Option<PathBuf> {
    let maps = std::fs::read_to_string("/proc/self/maps").ok()?;
    maps.lines()
        .filter_map(|l| l.split_whitespace().nth(5))
        .find(|p| {
            Path::new(p)
                .file_name()
                .is_some_and(|n| n.to_string_lossy().starts_with("libortools"))
        })
        .map(PathBuf::from)
}

/// The version of the library at `lib`, from the files around it (see the module docs).
pub fn version_near(lib: &Path) -> String {
    let dir = lib.parent().unwrap_or(Path::new("."));
    if let Some(v) = std::fs::read_to_string(dir.join("VERSION"))
        .ok()
        .and_then(|t| t.lines().next().map(|l| l.trim().to_string()))
        .filter(|v| !v.is_empty())
    {
        return v;
    }
    // site-packages/ortools/.libs/libortools.so.9 -> site-packages/ortools-<v>.dist-info
    if let Some(site) = dir.parent().and_then(Path::parent)
        && let Ok(rd) = std::fs::read_dir(site)
    {
        for e in rd.filter_map(Result::ok) {
            let n = e.file_name().to_string_lossy().into_owned();
            if let Some(v) = n
                .strip_prefix("ortools-")
                .and_then(|r| r.strip_suffix(".dist-info"))
            {
                return v.to_string();
            }
        }
    }
    // libortools.so.9.15.6755, followed through the symlink
    if let Ok(real) = std::fs::canonicalize(lib)
        && let Some(n) = real.file_name().map(|n| n.to_string_lossy().into_owned())
        && let Some(v) = n.strip_prefix("libortools.so.").or_else(|| {
            n.strip_prefix("libortools.")
                .and_then(|r| r.strip_suffix(".dylib"))
        })
        && v.matches('.').count() == 2
    {
        return v.to_string();
    }
    "unknown".to_string()
}

/// What an agent or a person wants to know without solving.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Availability {
    pub available: bool,
    /// The loaded library's path.
    pub library: Option<PathBuf>,
    pub version: Option<String>,
    /// The search order, `path (origin)`.
    pub searched: Vec<String>,
    /// Why it is unavailable.
    pub error: Option<String>,
}

/// Loads the library if that has not happened yet and says what came of it.
pub fn availability() -> Availability {
    match library() {
        Ok(l) => Availability {
            available: true,
            library: Some(l.path.clone()),
            version: Some(l.version.clone()),
            searched: candidates().iter().map(Candidate::to_string).collect(),
            error: None,
        },
        Err(e) => Availability {
            available: false,
            library: None,
            version: None,
            searched: e.searched.clone(),
            error: Some(e.to_string()),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn env_alone_when_set() {
        let c = candidates_for(Some("/x/lib.so".into()), Some(Path::new("/a/bin/e")), None);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].origin, "env");
    }

    #[test]
    fn order_without_env() {
        let c = candidates_for(
            None,
            Some(Path::new("/a/bin/e")),
            Some(Path::new("/site/ortools")),
        );
        let o: Vec<_> = c.iter().map(|c| c.origin).collect();
        assert_eq!(o, ["install", "ortools-package", "ortools-wheel", "system"]);
        assert_eq!(
            c[0].path,
            Path::new("/a/lib/rusolve/ortools").join(LIB_NAME)
        );
        assert_eq!(c[1].path, Path::new("/site/ortools/.libs").join(LIB_NAME));
        assert_eq!(c[2].path, Path::new("/ortools/.libs").join(LIB_NAME));
        assert_eq!(c[3].path, PathBuf::from(LIB_NAME));
    }

    #[test]
    fn version_from_file_then_dist_info() {
        let d = std::env::temp_dir().join(format!("rusolve-rt-{}", std::process::id()));
        let libs = d.join("site/ortools/.libs");
        std::fs::create_dir_all(&libs).unwrap();
        std::fs::create_dir_all(d.join("site/ortools-9.15.6755.dist-info")).unwrap();
        let lib = libs.join("libortools.so.9");
        assert_eq!(version_near(&lib), "9.15.6755");
        std::fs::write(libs.join("VERSION"), "9.99.1\n").unwrap();
        assert_eq!(version_near(&lib), "9.99.1");
        let _ = std::fs::remove_dir_all(&d);
    }
}
