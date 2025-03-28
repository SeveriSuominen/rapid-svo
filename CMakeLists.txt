cmake_minimum_required(VERSION 3.8)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# adding Raylib
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # dont build the supplied examples
set(BUILD_GAMES    OFF CACHE BOOL "" FORCE) # dont build the supplied example games

FetchContent_Declare(
    nanobench
    GIT_REPOSITORY "https://github.com/martinus/nanobench"
    GIT_TAG "v4.3.11"
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(nanobench)

FetchContent_Declare(
    glm
    GIT_REPOSITORY "https://github.com/g-truc/glm"
    GIT_TAG "1.0.1"
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(glm)

FetchContent_Declare(
    libmorton
    GIT_REPOSITORY "https://github.com/Forceflow/libmorton"
    GIT_TAG "v0.2.12"
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(libmorton)

project(rapid-svo LANGUAGES CXX C)

set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE) 

# define runtime intermediate build directory, so that debugging the Visual Studio solution works as expected.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

# set export directory for CMake config.
set(CMAKE_INSTALL_BINARY_DIR    "bin/")
set(CMAKE_INSTALL_LIBRARY_DIR   "bin/")
set(CMAKE_INSTALL_INCLUDE_DIR   "include/")
set(CMAKE_INSTALL_EXPORT_DIR    "cmake/")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)

add_compile_options(-fdiagnostics-color)

add_executable(benchmark "benchmark.cpp")
add_dependencies(benchmark glm)
add_dependencies(benchmark libmorton)
add_dependencies(benchmark nanobench)

message(INFO ${libmorton_SOURCE_DIR}/include)

target_include_directories(benchmark
PUBLIC
    ${glm_SOURCE_DIR}
    ${libmorton_SOURCE_DIR}/include
    ${nanobench_SOURCE_DIR}/src/include
)
