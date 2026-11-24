//! Generates the protobuf types of CP-SAT from the vendored 9.15 `.proto` files. OR-Tools itself
//! is not needed here: the library is opened at run time (`src/runtime.rs`).

fn main() {
    println!("cargo:rerun-if-changed=proto");
    let fds = protox::compile(["cp_model.proto", "sat_parameters.proto"], ["proto"])
        .expect("compile the OR-Tools .proto files");
    prost_build::Config::new()
        .compile_fds(fds)
        .expect("generate the CP-SAT protobuf types");
}
