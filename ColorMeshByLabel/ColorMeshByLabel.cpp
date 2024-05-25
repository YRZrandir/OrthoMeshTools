#include <filesystem>
#include <iostream>
#include <CGAL/boost/graph/IO/OBJ.h>
#include "../Polyhedron.h"
#include "CGAL/Simple_cartesian.h"

bool ColorMeshByLabel( std::string input_file, std::string input_labels, std::string output_file )
{
    std::vector<CGAL::Simple_cartesian<double>::Point_3> vertices;
    std::vector<TTriangle<size_t>> indices;
    if(input_file.ends_with(".obj"))
    {
        LoadVFObj<CGAL::Simple_cartesian<double>, size_t>( input_file, vertices, indices );
    }
    else
    {
        LoadVFAssimp<CGAL::Simple_cartesian<double>, size_t>(input_file, vertices, indices);
    }
    
    auto labels = LoadLabels(input_labels);
    if(vertices.size() != labels.size())
    {
        std::cout << "Error: number of vertices != number of labels." << std::endl;
        return false;
    }
    WriteVFAssimp<CGAL::Simple_cartesian<double>, size_t>(output_file, vertices, indices, labels);
    return true;
}

