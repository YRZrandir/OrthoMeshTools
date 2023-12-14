#include <filesystem>
#include <iostream>
#include <CGAL/boost/graph/IO/OBJ.h>
#include "../Polyhedron.h"
#include "CGAL/Simple_cartesian.h"

std::vector<int> LoadLabels(std::string path)
{
    using namespace nlohmann;
    std::ifstream label_ifs( path );
    json data = json::parse( label_ifs );
    if (data.find( "labels" ) == data.end())
    {
        std::cout << "Invalid Json" << std::endl;
        std::terminate();
    }
    std::vector<int> labels = data["labels"].get<std::vector<int>>();
    return labels;
}

bool ColorMeshByLabel( std::string input_file, std::string input_labels, std::string output_file )
{
    std::vector<CGAL::Simple_cartesian<double>::Point_3> vertices;
    std::vector<TTriangle<size_t>> indices;
    if(!LoadVFAssimp<CGAL::Simple_cartesian<double>, size_t>(input_file, vertices, indices))
    {
        std::cout << "Error: failed to load mesh." << std::endl;
        return false;
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

#ifndef FOUND_PYBIND11
int main(int argc, char* argv[])
{
    std::string input;
    std::string labels;
    std::string output;

    for(int i = 1; i < argc; i++)
    {
        if(std::strcmp(argv[i], "-i") == 0)
        {
            input = std::string(argv[i+1]);
        }
        else if(std::strcmp(argv[i], "-l") == 0)
        {
            labels = std::string(argv[i+1]);
        }
        else if(std::strcmp(argv[i], "-o") == 0)
        {
            output = std::string(argv[i+1]);
        }
    }

    std::cout << "input: " << input << std::endl;
    std::cout << "labels: " << labels << std::endl;
    std::cout << "output: " << output << std::endl;

    if(ColorMeshByLabel(input, labels, output))
    {
        return 0;
    }
    return -1;
}
#endif
