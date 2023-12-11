#ifndef MESH_FIX_H
#define MESH_FIX_H
#include <string>

void FixMesh(
    std::string path,
    std::string output_path, 
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry);

#endif