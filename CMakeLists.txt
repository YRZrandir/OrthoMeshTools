cmake_minimum_required(VERSION 3.24)
set(CMAKE_CXX_STANDARD 20)
project("OrthoMeshTools" C CXX)

include(FetchContent)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)
set(FETCHCONTENT_QUIET OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-fpermissive -fPIC -fopenmp")
option(BUILD_EXECUTABLES "" ON)
option(BUILD_PYTHON_LIB "" OFF)
add_compile_definitions(NOMINMAX)
add_compile_options("/EHsc")

FetchContent_Declare(Eigen3 URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
 CMAKE_ARGS -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_TEST_NOQT=ON
)
FetchContent_MakeAvailable(Eigen3)

FetchContent_Declare(assimp URL https://github.com/assimp/assimp/archive/refs/tags/v5.4.0.tar.gz
 CMAKE_ARGS -DASSIMP_BUILD_TESTS=OFF -DASSIMP_INSTALL=OFF -DASSIMP_BUILD_ASSIMP_TOOLS=OFF -DASSIMP_BUILD_ZLIB=OFF)
FetchContent_MakeAvailable(assimp)

FetchContent_Declare(nlohmann_json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(nlohmann_json)

# find_package(Eigen3 REQUIRED)
# ind_package(assimp REQUIRED)
# find_package(nlohmann_json REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Ceres)
include(CGAL_Eigen3_support)
include(CGAL_Ceres_support)

include_directories("./")
add_library(Ortho "Polyhedron.cpp" "Ortho.cpp" "Polyhedron.h" "Ortho.h")
target_link_libraries(Ortho PUBLIC CGAL::CGAL CGAL::Eigen3_support nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX assimp::assimp)

if (BUILD_PYTHON_LIB)
    FetchContent_Declare(pybind11 URL https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.tar.gz)
    FetchContent_MakeAvailable(pybind11)
    # find_package(pybind11 CONFIG)
    add_compile_definitions(FOUND_PYBIND11)
    pybind11_add_module(OrthoMeshTools SHARED "PyBind.cpp" "MeshFix/MeshFix.cpp" "SegClean/SegClean.cpp"
            "ReSegment/ReSegment.cpp" "HoleMerge/HoleMerge.cpp" "GumTrimLine/GumTrimLine.cpp" "ColorMeshByLabel/ColorMeshByLabel.cpp"
            "FakeToothRoot/FakeToothRoot.cpp" "OrthoScanBase/OrthoScanBase.cpp"
    )
    target_link_libraries(OrthoMeshTools PUBLIC Ortho pybind11::module)
    target_compile_definitions(OrthoMeshTools PUBLIC FOUND_PYBIND11)
endif ()

if (BUILD_EXECUTABLES)
    # find_package(argparse REQUIRED)
    FetchContent_Declare(argparse URL https://github.com/p-ranav/argparse/archive/refs/tags/v3.0.tar.gz)
    FetchContent_MakeAvailable(argparse)
    add_executable(MeshFix "MeshFix/MeshFix.cpp" "MeshFix/MeshFix.h" "MeshFix/MeshFixApp.cpp")
    target_compile_definitions(MeshFix PRIVATE _USE_MATH_DEFINES)
    target_link_libraries(MeshFix PUBLIC Ortho argparse::argparse)

    add_executable(SegClean "SegClean/SegClean.cpp" "SegClean/SegClean.h" "SegClean/SegCleanApp.cpp")
    target_link_libraries(SegClean PUBLIC Ortho argparse::argparse)

    add_executable(ReSegment "ReSegment/ReSegment.cpp" "ReSegment/ReSegment.h" "ReSegment/ReSegmentApp.cpp")
    target_link_libraries(ReSegment PUBLIC Ortho argparse::argparse)

    add_executable(HoleMerge "HoleMerge/HoleMerge.cpp" "HoleMerge/HoleMerge.h" "HoleMerge/HoleMergeApp.cpp")
    target_link_libraries(HoleMerge PUBLIC Ortho argparse::argparse)

    add_executable(GumTrimLine "GumTrimLine/GumTrimLine.cpp" "GumTrimLine/GumTrimLineApp.cpp" "GumTrimLine/GumTrimLine.h" "MeshFix/MeshFix.cpp")
    target_link_libraries(GumTrimLine PUBLIC Ortho argparse::argparse OpenMP::OpenMP_CXX)

    add_executable(ColorMeshByLabel "ColorMeshByLabel/ColorMeshByLabel.cpp" "ColorMeshByLabel/ColorMeshByLabel.h" "ColorMeshByLabel/ColorMeshByLabelApp.cpp")
    target_link_libraries(ColorMeshByLabel PUBLIC Ortho argparse::argparse)

    add_executable(OrthoScanDeform "OrthoScanDeform/OrthoScanDeformApp.cpp" "OrthoScanDeform/OrthoScanDeform.h" "MeshFix/MeshFix.cpp")
    target_link_libraries(OrthoScanDeform PUBLIC Ortho argparse::argparse)
    if (TARGET CGAL::Ceres_support)
        target_link_libraries(OrthoScanDeform PUBLIC CGAL::Ceres_support)
    endif ()

    add_executable(OrthoScanBase "OrthoScanBase/OrthoScanBase.cpp" "OrthoScanBase/OrthoScanBaseApp.cpp" "MeshFix/MeshFix.cpp")
    target_link_libraries(OrthoScanBase PUBLIC Ortho argparse::argparse)

    add_executable(FakeToothRoot "FakeToothRoot/FakeToothRoot.cpp" "FakeToothRoot/FakeToothRoot.h" "FakeToothRoot/FakeToothRootApp.cpp" "MeshFix/MeshFix.cpp")
    target_link_libraries(FakeToothRoot PUBLIC Ortho argparse::argparse)

    add_executable(OrthoScanPath "OrthoScanPath/OrthoScanPath.cpp" "OrthoScanPath/OrthoScanPathApp.cpp")
    target_link_libraries(OrthoScanPath PUBLIC Ortho argparse::argparse)

    add_executable(MoveTeeth "MoveTeeth/MoveTeeth.cpp" "OrthoScanPath/OrthoScanPath.cpp")
    target_link_libraries(MoveTeeth PUBLIC Ortho argparse::argparse)

    add_executable(GumLineCompare "GumLineCompare/GumLineCompare.cpp")
    target_link_libraries(GumLineCompare PUBLIC Ortho argparse::argparse)

    add_executable(ColorMeshByCurv "ColorMeshByCurv/ColorMeshByCurv.cpp")
    target_link_libraries(ColorMeshByCurv PUBLIC Ortho argparse::argparse)
endif ()
