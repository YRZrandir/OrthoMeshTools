cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project("OrthoMeshTools" C CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-fpermissive -fPIC -fopenmp")
option(BUILD_EXECUTABLES "" ON)
option(BUILD_PYTHON_LIB "" OFF)

add_compile_definitions(NOMINMAX)
add_compile_options("/EHsc")
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(assimp REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Ceres)
include(CGAL_Eigen3_support)
include(CGAL_Ceres_support)

include_directories("./")

add_library(Ortho "Polyhedron.cpp" "Ortho.cpp" "Polyhedron.h" "Ortho.h")
target_link_libraries(Ortho PUBLIC CGAL::CGAL CGAL::Eigen3_support nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX assimp::assimp)

if (BUILD_PYTHON_LIB)
    # find pybind11 dependency when needed
    find_package(pybind11 CONFIG)
    if (pybind11_FOUND)
        add_compile_definitions(FOUND_PYBIND11)
        pybind11_add_module(OrthoMeshTools SHARED "PyBind.cpp" "MeshFix/MeshFix.cpp" "SegClean/SegClean.cpp"
                "ReSegment/ReSegment.cpp" "HoleMerge/HoleMerge.cpp" "GumTrimLine/GumTrimLine.cpp" "ColorMeshByLabel/ColorMeshByLabel.cpp"
                "FakeToothRoot/FakeToothRoot.cpp" "OrthoScanBase/OrthoScanBase.cpp"
        )
        target_link_libraries(OrthoMeshTools PUBLIC Ortho pybind11::module)
        target_compile_definitions(OrthoMeshTools PUBLIC FOUND_PYBIND11)
    else ()
        message(WARNING "Failed to find pybind11.")
    endif ()
endif ()

if (BUILD_EXECUTABLES)
    find_package(argparse REQUIRED)

    add_executable(MeshFix "MeshFix/MeshFix.cpp" "MeshFix/MeshFix.h" "MeshFix/MeshFixApp.cpp")
    target_compile_definitions(MeshFix PRIVATE _USE_MATH_DEFINES)
    target_link_libraries(MeshFix PUBLIC Ortho argparse::argparse)

    add_executable(SegClean "SegClean/SegClean.cpp" "SegClean/SegClean.h")
    target_link_libraries(SegClean PUBLIC Ortho argparse::argparse)

    add_executable(ReSegment "ReSegment/ReSegment.cpp" "ReSegment/ReSegment.h")
    target_link_libraries(ReSegment PUBLIC Ortho argparse::argparse)

    add_executable(HoleMerge "HoleMerge/HoleMerge.cpp" "HoleMerge/HoleMerge.h")
    target_link_libraries(HoleMerge PUBLIC Ortho argparse::argparse)

    add_executable(GumTrimLine "GumTrimLine/GumTrimLine.cpp" "GumTrimLine/GumTrimLineApp.cpp" "GumTrimLine/GumTrimLine.h" "MeshFix/MeshFix.cpp")
    target_link_libraries(GumTrimLine PUBLIC Ortho argparse::argparse)

    add_executable(ColorMeshByLabel "ColorMeshByLabel/ColorMeshByLabel.cpp" "ColorMeshByLabel/ColorMeshByLabel.h")
    target_link_libraries(ColorMeshByLabel PUBLIC Ortho argparse::argparse)

    add_executable(OrthoScanDeform "OrthoScanDeform/OrthoScanDeformApp.cpp" "OrthoScanDeform/OrthoScanDeform.h" "MeshFix/MeshFix.cpp")
    target_link_libraries(OrthoScanDeform PUBLIC Ortho argparse::argparse)
    if (TARGET CGAL::Ceres_support)
        target_link_libraries(OrthoScanDeform PUBLIC CGAL::Ceres_support)
    endif ()

    add_executable(OrthoScanBase "OrthoScanBase/OrthoScanBase.cpp" "MeshFix/MeshFix.cpp")
    target_link_libraries(OrthoScanBase PUBLIC Ortho argparse::argparse)

    add_executable(FakeToothRoot "FakeToothRoot/FakeToothRoot.cpp" "FakeToothRoot/FakeToothRoot.h" "MeshFix/MeshFix.cpp")
    target_link_libraries(FakeToothRoot PUBLIC Ortho argparse::argparse)

    add_executable(OrthoScanPath "OrthoScanPath/OrthoScanPath.cpp" "OrthoScanPath/OrthoScanPathApp.cpp")
    target_link_libraries(OrthoScanPath PUBLIC Ortho argparse::argparse)

    add_executable(MoveTeeth "MoveTeeth/MoveTeeth.cpp" "OrthoScanPath/OrthoScanPath.cpp")
    target_link_libraries(MoveTeeth PUBLIC Ortho argparse::argparse)
endif ()
