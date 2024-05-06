# CMake file for fastjet codes used for benchmarks and validation
#
# (C) CERN, 2024
# MIT License

cmake_minimum_required(VERSION 3.12)
project(FastJetBenchmarks)

# Add package utilities to CMake path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

find_package(HepMC3 REQUIRED)
find_package(FastJet REQUIRED)

set(CMAKE_CXX_STANDARD 17)

# fastjet-inclusive runs reconstruction then finds inclusive jets
add_executable(fastjet-inclusive
    src/fastjet-inclusive.cc
    src/fastjet-utils.cc
)

target_include_directories(fastjet-inclusive PRIVATE
    ${FASTJET_INCLUDE_DIRS}
)

target_link_libraries(fastjet-inclusive 
    ${HEPMC3_LIBRARIES}
    ${FASTJET_LIBRARIES}
)

# fastjet-exlusive runs reconstruction then finds exclusive jets
add_executable(fastjet-exclusive
    src/fastjet-exclusive.cc
    src/fastjet-utils.cc
)

target_include_directories(fastjet-exclusive PRIVATE
    ${FASTJET_INCLUDE_DIRS}
)

target_link_libraries(fastjet-exclusive 
    ${HEPMC3_LIBRARIES}
    ${FASTJET_LIBRARIES}
)

