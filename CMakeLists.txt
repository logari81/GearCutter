cmake_minimum_required(VERSION 3.10)
project(gear_cutter)

#C++
enable_language(CXX)
add_library(gear_cutter_cpp_lib ${CMAKE_CURRENT_SOURCE_DIR}/cpp/src/gear_cutter_impl.cpp)
target_compile_options(gear_cutter_cpp_lib PRIVATE -fPIC)
set_target_properties(gear_cutter_cpp_lib PROPERTIES OUTPUT_NAME "gear_cutter_cpp")
target_link_libraries(gear_cutter_cpp_lib PUBLIC)
add_executable(gear_cutter_cpp_exe ${CMAKE_CURRENT_SOURCE_DIR}/cpp/src/gear_cutter_main.cpp)
set_target_properties(gear_cutter_cpp_exe PROPERTIES OUTPUT_NAME "gear_cutter_cpp")
target_link_libraries(gear_cutter_cpp_exe gear_cutter_cpp_lib)


#Fortran
enable_language(Fortran)
add_executable(gear_cutter_fortran_exe ${CMAKE_CURRENT_SOURCE_DIR}/fortran/src/gear_cutter_main.f90)
set_target_properties(gear_cutter_fortran_exe PROPERTIES OUTPUT_NAME "gear_cutter_fortran")

#Rust
# Build the Rust executable (gear_cutter_rs)
add_custom_target(gear_cutter_rust_exe ALL
  COMMAND ${CMAKE_COMMAND} -E
    env CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR}/rust
    cargo build --release --bin gear_cutter_rs
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rust  COMMENT "Building Rust executable: gear_cutter_rs"
)
# Build the static Rust library (libgear_cutter.rlib)
add_custom_target(gear_cutter_rust_lib ALL
  COMMAND ${CMAKE_COMMAND} -E
    env CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR}/rust
    cargo build --release --lib
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rust
  COMMENT "Building static Rust library: libgear_cutter.rlib"
)
# Build the Python extension (gear_cutter_rs.so)
add_custom_target(pygear_cutter_rust ALL
  COMMAND ${CMAKE_COMMAND} -E
    env CARGO_TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR}/rust
    cargo build --manifest-path ${CMAKE_CURRENT_SOURCE_DIR}/rust/python/Cargo.toml
    --release --lib
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rust
  COMMENT "Building Python extension: gear_cutter_rs.so"
)
# Copy the generated files to appropriate locations
add_custom_command(TARGET gear_cutter_rust_exe POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/rust/release/gear_cutter_rs
    ${CMAKE_CURRENT_BINARY_DIR}/gear_cutter_rs
  COMMENT "Copying gear_cutter_rs executable"
)
add_custom_command(TARGET gear_cutter_rust_lib POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/rust/release/libgear_cutter.rlib
    ${CMAKE_CURRENT_BINARY_DIR}/libgear_cutter.rlib
  COMMENT "Copying libgear_cutter.rlib"
)
add_custom_command(TARGET pygear_cutter_rust POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/rust/release/libpython_bindings.so
    ${CMAKE_CURRENT_BINARY_DIR}/gear_cutter_rs.so
  COMMENT "Copying Python extension gear_cutter_rs.so"
)

#Python cpp wrapper
find_package(Python3 COMPONENTS Interpreter Development)
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_library(pygear_cutter_cpp MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp/python/gear_cutter_python.cpp
)
target_include_directories(pygear_cutter_cpp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp/src
  ${Python3_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIR}
)
set_target_properties(pygear_cutter_cpp PROPERTIES PREFIX "")
set_target_properties(pygear_cutter_cpp PROPERTIES OUTPUT_NAME "gear_cutter")
target_link_libraries(pygear_cutter_cpp PRIVATE gear_cutter_cpp_lib)