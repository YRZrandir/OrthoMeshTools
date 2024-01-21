### OrthoMeshTools
A set of tools for ortho scan mesh processing.

See `PyBind.cpp` for all python interface.

### Dependencies
`CGAL, Eigen, Assimp, nlohmann_json`

`Ceres` Is needed for some programs.

`pybind11` is needed to build python binding.

`argparse` is needed to build executables.

### MeshFix
The MeshFix program performs the following steps:
1. Detecting non-manifold ( and self-intersecting ) parts of the input mesh.
2. Deleting those faces & vertices.
3. Detecting again. If no invalid elements are found, goto 4, otherwise goto 2. 
4. Close the holes.

**Arguments:**

`--input, -i` specify the input mesh file.

`--output, -o` specify the output mesh file.

`--small_component_threshold, -sc` Disconnected components whose face number are smaller than this will be removed.

`--fix_self_intersection` Remove self-intersecting faces.

`--smallhole_edge_num, -sh` Holes whose edge number is smaller than this value will be closed. Just use a very large number if you want to close all holes.

`--smallhole_size, -ss` If the hole's AABB is smaller than this value in all dimension, the hole will be closed. Just use a very large number if you want to close all holes.

`--refine, -r` Refine the patch after closing the hole.

`--max_retry, -m` Sometimes deleting faces can generate new invalid elements, in this case we try to fix the mesh again. This value specifies the max retry time.

**Example:**

1. Fix non-manifold, close all holes, and remove disconnected components that have less than 10 faces:

`MeshFix.exe -i bunny.obj -o bunny_fixed.obj -sh 999999999 -ss 999999999 -sc 10`

1. If the input mesh has some holes itself and they shouldn't be closed, please try adjusting `-sh` and `-ss` to filter them.

`MeshFix.exe -i bunny.obj -o bunny_fixed.obj -sc 10 -sh 256 -ss 34.5`

**Fix mesh with labels** 

We allow each vertex to have a integer label, the labels should be input as a `.json` file looks like: `{"labels": [0, 1, 1, 1, 2, 2]}`. 
   
   The label number must be the same as the vertex number, and the order of labels must match the vertex order in mesh file. Specify the input label file using `-li or --input_label`, and specify the output label file by `-lo or --output_label`.

   When closing hole with `--refine`, new vertices will be added. Their labels are computed according to the nearest 'labeled' vertex. The c++ interface returns all vertices & faces of the hole patch, please use it if more control is needed.

   Note that we need the vertex order to work with labels, so formats like `.stl` cannot be used. For other formats, Assimp can keep the vertex order most of the time, but fails in some situation. (I may add new mesh importer in the furture to solve this).