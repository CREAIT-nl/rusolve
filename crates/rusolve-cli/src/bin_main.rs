//! The `rusolve-cli` executable: the crate root (`main.rs`) is the library, so the Python
//! bindings can reuse the subcommand cores (`rusolve_cli::api`).

fn main() -> std::process::ExitCode {
    rusolve_cli::main_entry()
}
