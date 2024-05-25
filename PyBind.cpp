#ifdef FOUND_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ColorMeshByLabel/ColorMeshByLabel.h"
#include "GumTrimLine/GumTrimLine.h"
#include "HoleMerge/HoleMerge.h"
#include "MeshFix/MeshFix.h"
#include "ReSegment/ReSegment.h"
#include "SegClean/SegClean.h"
#include "FakeToothRoot/FakeToothRoot.h"
#include "OrthoScanBase/OrthoScanBase.h"
#include "Polyhedron.h"

namespace py = pybind11;

PYBIND11_MODULE(OrthoMeshTools, m)
{
    m.doc() = "A set of tools for ortho scan meshes.";

    m.def("ColorMeshByLabel", &ColorMeshByLabel, "Give each mesh vertex a color according to its label.",
        py::arg("input_mesh"),
        py::arg("input_labels"),
        py::arg("output_mesh"));

    m.def("SegClean", &SegCleanF,
        "Try to make sure each label is applied to single connected component. Areas that are smaller than size_threshold are processsed.",
        py::arg("input_mesh"),
        py::arg("input_labels"),
        py::arg("output_labels"),
        py::arg("size_threshold"));

    m.def("FixMesh", &FixMeshFile, "Fix non-manifold vertices & edges",
        py::arg("path"),
        py::arg("output_path"),
        py::arg("cfg"));

    m.def("FixMeshWithLabel", &FixMeshFileWithLabel, "Fix non-manifold vertices & edges and output manipulated vertex label file.",
        py::arg("path"),
        py::arg("output_path"),
        py::arg("input_label"),
        py::arg("output_label"),
        py::arg("cfg"));

    m.def("ReSegment", &ReSegmentLabels, "ReSegment the mesh using the provided splitlines.",
        py::arg("input_mesh"),
        py::arg("splitlines"),
        py::arg("splitline_labels"),
        py::arg("output_json"),
        py::arg("intersection_width"),
        py::arg("cutface_orit_smooth"),
        py::arg("upper"));

    m.def("HoleMerge", &HoleMerge, "Merge the holes on the mesh",
        py::arg("input_mesh"),
        py::arg("output_mesh"),
        py::arg("threshold"));

    m.def("GumTrimLine", &GumTrimLine, "Get gum trim line",
        py::arg("input_mesh"),
        py::arg("input_labels"),
        py::arg("frame_file"),
        py::arg("output_file"),
        py::arg("smooth"),
        py::arg("fix"),
        py::arg("debug_output"));

    m.def("FakeToothRoot", &FakeToothRoot, "Generate fake tooth root for tooth crowns",
        py::arg("input_mesh"),
        py::arg("output_path"),
        py::arg("frame_path"),
        py::arg("label_path"),
        py::arg("fair"));

    m.def("GenerateGum", &GenerateGum, "Generate fake gum",
        py::arg("input_file"),
        py::arg("input_label"),
        py::arg("crown_frame"),
        py::arg("output_gum"));

    py::register_local_exception<IOError>(m, "IOError");
    py::register_local_exception<MeshError>(m, "MeshError");
    py::register_local_exception<AlgError>(m, "AlgError");
}
#endif
