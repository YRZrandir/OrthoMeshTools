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

