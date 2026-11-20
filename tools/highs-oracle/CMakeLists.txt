# HiGHS internal-class oracle for the rusolve-ekk faithful port.
#
# Builds HiGHS (a patched *copy* of reference/HiGHS, see gen_fixtures.sh and
# highs-oracle-hooks.patch) as a static, serial library and links the oracle
# executable against it, with HiGHS' internal headers on the include path.
#
#   cmake -S tools/highs-oracle -B target/highs-oracle/build \
#         -DHIGHS_SRC=target/highs-oracle/highs-src -DCMAKE_BUILD_TYPE=Release
cmake_minimum_required(VERSION 3.18)
project(highs_oracle LANGUAGES C CXX)

if(NOT HIGHS_SRC)
  message(FATAL_ERROR "Set -DHIGHS_SRC=<patched copy of reference/HiGHS>")
endif()
get_filename_component(HIGHS_SRC "${HIGHS_SRC}" ABSOLUTE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Bit-exactness: forbid floating-point contraction (a*b+c -> fma). Apple clang
# and recent GCC/clang default to -ffp-contract=on, which fuses multiply-adds
# inside expressions on arm64 / x86-64-v3. Rust never contracts implicitly, so
# the oracle must not either, or the fixtures would not be reproducible by a
# faithful (non-fma) Rust translation. See README.md "Floating point".
# ORACLE_FP_CONTRACT=on exists only for the contraction-sensitivity check in
# gen_fixtures.sh --fp-contract-check; fixtures are always generated with off.
set(ORACLE_FP_CONTRACT "off" CACHE STRING "-ffp-contract value (off|on|fast)")
set(ORACLE_FP_CONTRACT_EFFECTIVE "${ORACLE_FP_CONTRACT}")
# GCC before 14 accepts -ffp-contract=on but implements it as off, so the
# contraction check would compare two identical builds; fast is the GCC
# setting that actually fuses (within and across statements).
if(ORACLE_FP_CONTRACT STREQUAL "on" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
  message(STATUS "GCC ${CMAKE_CXX_COMPILER_VERSION}: -ffp-contract=on is a no-op, using fast")
  set(ORACLE_FP_CONTRACT_EFFECTIVE "fast")
endif()
set(ORACLE_FP_FLAGS "-ffp-contract=${ORACLE_FP_CONTRACT_EFFECTIVE}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ORACLE_FP_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ORACLE_FP_FLAGS}")

# HiGHS build switches: static, serial, no extras.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_CXX_EXE OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_EXTRAS_LIB OFF CACHE BOOL "" FORCE)
set(ZLIB OFF CACHE BOOL "" FORCE)
set(FAST_BUILD ON CACHE BOOL "" FORCE)
set(HIGHS_NO_DEFAULT_THREADS ON CACHE BOOL "" FORCE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF CACHE BOOL "" FORCE)
# PDLP (W3 fixtures): cuPDLP-C and HiPDLP are compiled for the CPU only
# (HConfig.h: CUPDLP_CPU defined, CUPDLP_GPU undefined; no CUDA toolkit).
set(CUPDLP_GPU OFF CACHE BOOL "" FORCE)
set(HIGHS_GPU_LIB OFF CACHE BOOL "" FORCE)

# HiPO (W4 fixtures): a second oracle
# configuration, in its own build dir (gen_fixtures.sh: build-hipo), with
# HIPO=ON and the BLAS pinned to refblas/refblas.c, the netlib-order
# reference CBLAS that crates/rusolve-ekk/src/ipm/hipo/blas.rs mirrors. The
# default (OFF) configuration, and every fixture it writes, is unchanged.
#
# The pin presets FindHipoDeps' cache-free state instead of passing
# -DBLAS_LIBRARIES: highs_configure_blas() returns early when
# HIGHS_BLAS_CONFIGURED is set, and highs_link_blas(highs_extras) then links
# HIGHS_BLAS_TARGET, so CMake's FindBLAS (which clears BLAS_LIBRARIES and
# searches for OpenBLAS on Linux) never runs and no other BLAS can be picked
# up. METIS, AMD and RCM are always the extern/ sources compiled into
# highs_extras.
option(ORACLE_HIPO "Build HiGHS with HIPO=ON, linked against refblas" OFF)
if(ORACLE_HIPO)
  add_library(refblas STATIC refblas/refblas.c refblas/refblas.h)
  set_target_properties(refblas PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON)
  # Pinned whatever ORACLE_FP_CONTRACT / CMAKE_BUILD_TYPE say: these flags
  # are part of refblas' definition (the Rust mirror does no contraction).
  # The later -O2 overrides the Release -O3.
  target_compile_options(refblas PRIVATE -O2 -ffp-contract=off)
  target_include_directories(refblas PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/refblas>)
  # highs_extras is installed with EXPORT highs-targets and links refblas
  # PUBLIC, so refblas must be in that export set for the generate step.
  install(TARGETS refblas EXPORT highs-targets ARCHIVE DESTINATION lib)

  add_executable(test_refblas EXCLUDE_FROM_ALL refblas/test_refblas.c)
  set_target_properties(test_refblas PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
  target_compile_options(test_refblas PRIVATE -O2 -ffp-contract=off)
  target_link_libraries(test_refblas PRIVATE refblas)

  set(HIPO ON CACHE BOOL "" FORCE)
  set(BUILD_OPENBLAS OFF CACHE BOOL "" FORCE)
  set(HIGHS_BLAS_CONFIGURED TRUE)
  set(HIGHS_BLAS_TARGET refblas)
  set(HIGHS_BLAS_INCLUDE_DIRS "")
  set(HIGHS_BLAS_COMPILE_DEFINITION "")
  set(HIGHS_BLAS_VENDOR "rusolve refblas")
  set(HIGHS_BLAS_VERSION "netlib-order")
  set(HIGHS_BLAS_LICENSE "MIT")
  message(STATUS "highs-oracle: HIPO=ON, BLAS pinned to refblas (${CMAKE_CURRENT_SOURCE_DIR}/refblas/refblas.c)")
else()
  set(HIPO OFF CACHE BOOL "" FORCE)
endif()

add_subdirectory(${HIGHS_SRC} highs-build EXCLUDE_FROM_ALL)

# FAST_BUILD batches the cuPDLP-C sources into one unity C file, where
# cupdlp_proj.c includes HighsOracleW3Dump.h before cupdlp_solver.c defines
# HIGHS_ORACLE_W3_DUMP_IMPL, so the include guard drops the writer
# definitions and highs-oracle fails to link (undefined highsOracleW3File).
# Compiling cupdlp_solver.c as its own translation unit keeps the patch-4
# hash unchanged; without contraction the generated code computes the same
# values.
if(EXISTS ${HIGHS_SRC}/highs/pdlp/HighsOracleW3Dump.h)
  set_source_files_properties(${HIGHS_SRC}/highs/pdlp/cupdlp/cupdlp_solver.c
    DIRECTORY ${HIGHS_SRC}/highs PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif()

add_executable(highs-oracle
  src/main.cpp
  src/oracle.h
  src/emit.cpp
  src/unit.cpp
  src/factor.cpp
  src/basiclu.cpp
  src/ipx.cpp
  src/simplex.cpp
  src/w1.cpp
  src/w3.cpp
  src/w4.cpp
  src/w5.cpp
  src/g1.cpp
  src/clqn.cpp)
set_target_properties(highs-oracle PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)
target_include_directories(highs-oracle PRIVATE
  ${HIGHS_SRC}/highs
  ${HIGHS_SRC}/extern
  ${HIGHS_BINARY_DIR})
target_link_libraries(highs-oracle PRIVATE highs)

# Stand-alone drivers on the public Highs API, linked against the same
# patched, serial, non-contracted library: mip_trace (the MIP trace) and
# api_seq (the model-modification sequence). Built on request, e.g.
#   cmake --build target/highs-oracle/build --target mip_trace api_seq
foreach(drv mip_trace api_seq)
  add_executable(${drv} EXCLUDE_FROM_ALL ${drv}.cpp)
  set_target_properties(${drv} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)
  target_include_directories(${drv} PRIVATE
    ${HIGHS_SRC}/highs
    ${HIGHS_SRC}/extern
    ${HIGHS_BINARY_DIR})
  target_link_libraries(${drv} PRIVATE highs)
endforeach()
