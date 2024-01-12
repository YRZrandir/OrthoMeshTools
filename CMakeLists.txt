cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project ("OrthoMeshTools" C CXX)

add_compile_definitions(NOMINMAX)

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(assimp REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(gflags REQUIRED)
find_package(pybind11 CONFIG)
find_package(Ceres)
include(CGAL_Eigen3_support)
include(CGAL_Ceres_support)

add_library(Ortho "Polyhedron.cpp" "Polyhedron.h" "Ortho.h")
target_link_libraries(Ortho PUBLIC CGAL::CGAL CGAL::Eigen3_support nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX assimp::assimp gflags)

add_executable(MeshFix "MeshFix/MeshFix.cpp" "MeshFix/MeshFix.h" "MeshFix/MeshFixApp.cpp")
target_compile_definitions(MeshFix PRIVATE _USE_MATH_DEFINES)
target_link_libraries(MeshFix PUBLIC Ortho)

add_executable(SegClean "SegClean/SegClean.cpp" "SegClean/SegClean.h")
target_link_libraries(SegClean PUBLIC Ortho)

add_executable(ReSegment "ReSegment/ReSegment.cpp" "ReSegment/ReSegment.h")
target_link_libraries(ReSegment PUBLIC Ortho)

add_executable(HoleMerge "HoleMerge/HoleMerge.cpp" "HoleMerge/HoleMerge.h")
target_link_libraries(HoleMerge PUBLIC Ortho)

add_executable(GumTrimLine "GumTrimLine/GumTrimLine.cpp"  "GumTrimLine/GumTrimLineApp.cpp" "GumTrimLine/GumTrimLine.h" "MeshFix/MeshFix.cpp")
target_link_libraries(GumTrimLine PUBLIC Ortho)

add_executable(ColorMeshByLabel "ColorMeshByLabel/ColorMeshByLabel.cpp" "ColorMeshByLabel/ColorMeshByLabel.h")
target_link_libraries(ColorMeshByLabel PUBLIC Ortho)

add_executable(OrthoScanDeform "OrthoScanDeform/OrthoScanDeformApp.cpp" "OrthoScanDeform/OrthoScanDeform.h" "MeshFix/MeshFix.cpp")
target_link_libraries(OrthoScanDeform PUBLIC Ortho)
if(TARGET CGAL::Ceres_support)
    target_link_libraries(OrthoScanDeform PUBLIC CGAL::Ceres_support)
endif()

if(pybind11_FOUND)
    pybind11_add_module(OrthoMeshTools "PyBind.cpp" "MeshFix/MeshFix.cpp" "SegClean/SegClean.cpp"
     "ReSegment/ReSegment.cpp" "HoleMerge/HoleMerge.cpp" "GumTrimLine/GumTrimLine.cpp" "ColorMeshByLabel/ColorMeshByLabel.cpp")
    target_link_libraries(OrthoMeshTools PUBLIC Ortho)
    target_compile_definitions(OrthoMeshTools PUBLIC FOUND_PYBIND11)
endif()