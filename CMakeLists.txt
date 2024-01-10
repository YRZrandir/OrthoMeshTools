cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project ("OrthoMeshTools" C CXX)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(assimp REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(pybind11 CONFIG)

add_executable(MeshFix "MeshFix/MeshFix.cpp" "MeshFix/MeshFix.h" "MeshFix/MeshFixApp.cpp" "Polyhedron.cpp" "Polyhedron.h")
target_compile_definitions(MeshFix PRIVATE _USE_MATH_DEFINES)
target_link_libraries(MeshFix CGAL::CGAL OpenMP::OpenMP_CXX assimp::assimp)

add_executable(SegClean "SegClean/SegClean.cpp" "SegClean/SegClean.h" "Polyhedron.cpp" "Polyhedron.h")
target_link_libraries(SegClean CGAL::CGAL nlohmann_json::nlohmann_json assimp::assimp)

add_executable(ReSegment "ReSegment/ReSegment.cpp" "ReSegment/ReSegment.h" "Polyhedron.cpp" "Polyhedron.h")
target_link_libraries(ReSegment PRIVATE CGAL::CGAL nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX assimp::assimp)

add_executable(HoleMerge "HoleMerge/HoleMerge.cpp" "HoleMerge/HoleMerge.h")
target_link_libraries(HoleMerge PRIVATE CGAL::CGAL assimp::assimp)

add_executable(GumTrimLine "GumTrimLine/GumTrimLine.cpp" "GumTrimLine/GumTrimLine.h" "MeshFix/MeshFix.cpp" "Polyhedron.cpp")
target_link_libraries(GumTrimLine PRIVATE CGAL::CGAL assimp::assimp OpenMP::OpenMP_CXX)

add_executable(ColorMeshByLabel "ColorMeshByLabel/ColorMeshByLabel.cpp" "ColorMeshByLabel/ColorMeshByLabel.h" "Polyhedron.cpp")
target_link_libraries(ColorMeshByLabel PRIVATE CGAL::CGAL assimp::assimp)

add_executable(OrthoScanDeform "OrthoScanDeform/OrthoScanDeformApp.cpp" "OrthoScanDeform/OrthoScanDeform.h" "MeshFix/MeshFix.cpp" "Polyhedron.cpp")
target_link_libraries(OrthoScanDeform PRIVATE CGAL::CGAL assimp::assimp)

if(pybind11_FOUND)
    pybind11_add_module(OrthoMeshTools "PyBind.cpp" "MeshFix/MeshFix.cpp" "SegClean/SegClean.cpp"
     "ReSegment/ReSegment.cpp" "HoleMerge/HoleMerge.cpp" "Polyhedron.cpp" "Polyhedron.h"
     "GumTrimLine/GumTrimLine.cpp" "ColorMeshByLabel/ColorMeshByLabel.cpp")
    target_link_libraries(OrthoMeshTools PRIVATE CGAL::CGAL OpenMP::OpenMP_CXX assimp::assimp nlohmann_json::nlohmann_json)
    target_compile_definitions(OrthoMeshTools PUBLIC FOUND_PYBIND11)
endif()