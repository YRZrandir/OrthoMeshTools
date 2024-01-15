cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project ("gumTrimLine" C CXX)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(assimp REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(pybind11 CONFIG)

if(pybind11_FOUND)
    pybind11_add_module(gumTrimLine "PyBind.cpp" "print.cpp" "MeshFix/MeshFix.cpp" "Polyhedron.cpp" "Polyhedron.h"
     "GumTrimLine/GumTrimLine.cpp")
    target_link_libraries(gumTrimLine PRIVATE CGAL::CGAL OpenMP::OpenMP_CXX assimp::assimp nlohmann_json::nlohmann_json)
    target_compile_definitions(gumTrimLine PUBLIC FOUND_PYBIND11)
endif()