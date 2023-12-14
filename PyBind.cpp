#ifdef FOUND_PYBIND11
#include <pybind11/pybind11.h>
#include "ColorMeshByLabel/ColorMeshByLabel.h"
#include "GumTrimLine/GumTrimLine.h"
#include "HoleMerge/HoleMerge.h"
#include "MeshFix/MeshFix.h"
#include "ReSegment/ReSegment.h"
#include "SegClean/SegClean.h"
#include "Polyhedron.h"

namespace py = pybind11;

PYBIND11_MODULE(OrthoMeshTools, m)
{
    m.doc() = "A set of tools for ortho scan meshes.";

    m.def("ColorMeshByLabel", &ColorMeshByLabel, "Give each mesh vertex a color according to its label.",
        py::arg("input_mesh"),
        py::arg("input_labels"),
        py::arg("output_mesh"));

    m.def("SegClean", &SegClean,
        "Try to make sure each label is applied to single connected component. Areas that are smaller than size_threshold are processsed.",
        py::arg("input_mesh"),
        py::arg("input_labels"),
        py::arg("output_labels"),
        py::arg("size_threshold"));
    
    m.def("FixMesh", &FixMesh, "Fix non-manifold vertices & edges",
        py::arg("path"),
        py::arg("output_path"),
        py::arg("keep_largest_connected_component"),
        py::arg("large_cc_threshold"),
        py::arg("fix_self_intersection"),
        py::arg("filter_small_holes"),
        py::arg("max_hole_edges"),
        py::arg("max_hole_diam"),
        py::arg("refine"),
        py::arg("max_retry"));

    m.def("FixMeshWithLabel", &FixMeshWithLabel, "Fix non-manifold vertices & edges and output manipulated vertex label file.",
        py::arg("path"),
        py::arg("output_path"),
        py::arg("input_label"),
        py::arg("output_label"),
        py::arg("keep_largest_connected_component"),
        py::arg("large_cc_threshold"),
        py::arg("fix_self_intersection"),
        py::arg("filter_small_holes"),
        py::arg("max_hole_edges"),
        py::arg("max_hole_diam"),
        py::arg("refine"),
        py::arg("max_retry"));  

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
        py::arg("output_file"),
        py::arg("smooth"));

    py::register_local_exception<IOError>(m, "IOError");
    py::register_local_exception<MeshError>(m, "MeshError");
    py::register_local_exception<AlgError>(m, "AlgError");
}
#endif