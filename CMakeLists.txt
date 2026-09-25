cmake_minimum_required(VERSION 3.16)
project(BFV_SlotInversion VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(SEAL 4.1 QUIET)
if(NOT SEAL_FOUND)
  find_package(SEAL 4.0 QUIET)
endif()
if(NOT SEAL_FOUND)
  find_package(SEAL 3.7 REQUIRED)
endif()
message(STATUS "Found SEAL: ${SEAL_VERSION}")

# ── Locate source files (src/ or project root) ───────────────────────────
function(find_src VAR FILENAME)
  if(EXISTS "${CMAKE_SOURCE_DIR}/src/${FILENAME}")
    set(${VAR} "${CMAKE_SOURCE_DIR}/src/${FILENAME}" PARENT_SCOPE)
  elseif(EXISTS "${CMAKE_SOURCE_DIR}/${FILENAME}")
    set(${VAR} "${CMAKE_SOURCE_DIR}/${FILENAME}" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Cannot find ${FILENAME}")
  endif()
endfunction()

find_src(SRC_PARAMS  bfv_params.cpp)
find_src(SRC_GALOIS  galois_ops.cpp)
find_src(SRC_BINEXP  binary_exp.cpp)
find_src(SRC_FMC     fmc.cpp)
find_src(SRC_NT      norm_tower.cpp)
find_src(SRC_MAIN    main.cpp)
find_src(SRC_BENCH   bench.cpp)

if(EXISTS "${CMAKE_SOURCE_DIR}/include/bfv_inversion.h")
  set(HDR_DIR "${CMAKE_SOURCE_DIR}/include")
else()
  set(HDR_DIR "${CMAKE_SOURCE_DIR}")
endif()

# ── Threads (for parallel FMC Phase 1) ───────────────────────────────────
find_package(Threads REQUIRED)

# ── Library ──────────────────────────────────────────────────────────────
add_library(fmc_lib STATIC
  ${SRC_PARAMS} ${SRC_GALOIS} ${SRC_BINEXP} ${SRC_FMC} ${SRC_NT})
target_include_directories(fmc_lib PUBLIC ${HDR_DIR})
target_link_libraries(fmc_lib PUBLIC SEAL::seal Threads::Threads)
target_compile_options(fmc_lib PRIVATE -Wall -Wextra
  $<$<CONFIG:Release>:-O3 -march=native>)

# ── Executables ───────────────────────────────────────────────────────────
add_executable(bfv_inversion ${SRC_MAIN})
target_link_libraries(bfv_inversion PRIVATE fmc_lib)
target_compile_options(bfv_inversion PRIVATE
  $<$<CONFIG:Release>:-O3 -march=native>)

add_executable(bfv_bench ${SRC_BENCH})
target_link_libraries(bfv_bench PRIVATE fmc_lib)
target_compile_options(bfv_bench PRIVATE
  $<$<CONFIG:Release>:-O3 -march=native>)
