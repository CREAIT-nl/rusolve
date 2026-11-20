// extension-module cdylibs must not link libpython; on macOS that needs
// `-undefined dynamic_lookup`, which this adds so a plain `cargo build` works too
// (maturin adds the same flags itself).
fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
