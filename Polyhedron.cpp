//
// Created by yrz on 7/7/22.
//
#include "Polyhedron.h"
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

template <>
class TPolyhedron<CGAL::Polyhedron_items_with_id_3, CGAL::Simple_cartesian<float>>;
template <>
class TPolyhedron<CGAL::Polyhedron_items_with_id_3, CGAL::Simple_cartesian<double>>;
template <>
class TPolyhedron<CGAL::Polyhedron_items_with_id_3, CGAL::Exact_predicates_inexact_constructions_kernel>;

std::vector<int> LoadLabels( std::string path )
{
    using namespace nlohmann;
    std::ifstream label_ifs( path );
    if(label_ifs.fail())
    {
        throw IOError("Cannot open file: " + path);
    }
    json data = json::parse( label_ifs );
    if (data.find( "labels" ) == data.end())
    {
        throw IOError("Cannot find key 'labels' in json file " + path);
    }
    return data["labels"].get<std::vector<int>>();
}