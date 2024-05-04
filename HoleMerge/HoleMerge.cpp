#include <iostream>
#include <chrono>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include "../Polyhedron.h"

namespace
{
  using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Polyhedron = TPolyhedron<CGAL::Polyhedron_items_with_id_3, KernelEpick>;
  using hHalfedge = Polyhedron::Halfedge_handle;
  using hVertex = Polyhedron::Vertex_handle;
  using hFacet = Polyhedron::Facet_handle;
  using Halfedge = Polyhedron::Halfedge;
  using CVertex = Polyhedron::Vertex;
  using Facet = Polyhedron::Facet;
  using iHalfedge = Polyhedron::Halfedge_iterator;
  using iVertex = Polyhedron::Vertex_iterator;
  using iFacet = Polyhedron::Facet_iterator;
  using Point_3 = Polyhedron::Point_3;
  using Vector_3 = Polyhedron::Traits::Vector_3;
  using SurfaceShortestPathTraits = CGAL::Surface_mesh_shortest_path_traits<KernelEpick, Polyhedron::Base>;
  using Surface_mesh_shortest_path = CGAL::Surface_mesh_shortest_path<SurfaceShortestPathTraits>;

  struct Sequence_collector
  {
    std::vector<hVertex> vertices;

    void operator()(hHalfedge he, double alpha)
    {
      vertices.push_back(he->vertex());
      vertices.push_back(he->prev()->vertex());
    }

    void operator()(hVertex v)
    {
      vertices.push_back(v);
    }

    void operator()(hFacet f, SurfaceShortestPathTraits::Barycentric_coordinates alpha)
    {
      vertices.push_back(f->halfedge()->vertex());
      vertices.push_back(f->halfedge()->next()->vertex());
      vertices.push_back(f->halfedge()->prev()->vertex());
    }
  };

  void Merge2Holes(const std::vector<hHalfedge> &hole0, const std::vector<hHalfedge> &hole1, Polyhedron &mesh)
  {
    std::pair<hHalfedge, hHalfedge> nearest_pair;
    double nearest_pair_dist = std::numeric_limits<double>::max();
    // Find nearest point pair
#pragma omp parallel for
    for (int i = 0; i < hole0.size(); i++)
    {
      hHalfedge hh0 = hole0[i];
      hVertex hv0 = hh0->vertex();
      Point_3 p0 = hv0->point();
      double min_dist_hv0 = std::numeric_limits<double>::max();
      hHalfedge nearest_to_hv0 = hole1[0];
      for (int j = 0; j < hole1.size(); j++)
      {
        hHalfedge hh1 = hole1[j];
        hVertex hv1 = hh1->vertex();
        Point_3 p1 = hv1->point();
        double dist_hv0_hv1 = CGAL::squared_distance(hv0->point(), hv1->point());

        if (dist_hv0_hv1 < min_dist_hv0)
        {
          nearest_to_hv0 = hh1;
          min_dist_hv0 = dist_hv0_hv1;
        }
      }

#pragma omp critical
      {
        if (min_dist_hv0 < nearest_pair_dist)
        {
          nearest_pair_dist = min_dist_hv0;
          nearest_pair = {hh0, nearest_to_hv0};
        }
      }
    }

    CGAL::set_halfedgeds_items_id(mesh);

    auto start = std::chrono::high_resolution_clock::now();

    CGAL::Surface_mesh_shortest_path<SurfaceShortestPathTraits> shortest_paths(mesh);
    shortest_paths.add_source_point(nearest_pair.first->vertex());

    Sequence_collector sc;
    shortest_paths.shortest_path_sequence_to_source_points(nearest_pair.second->vertex(), sc);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Distance Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start) << std::endl;

    std::unordered_set<hFacet> faces_to_erase;
    for (auto hv : sc.vertices)
    {
      for (auto nei : CGAL::faces_around_target(hv->halfedge(), mesh))
      {
        if (nei != nullptr)
          faces_to_erase.insert(nei);
      }
    }

    for (auto hh : faces_to_erase)
    {
      mesh.erase_facet(hh->halfedge());
    }
  }

  bool MergeLargestHoles(Polyhedron &mesh, int threshold)
  {
    std::vector<hHalfedge> border_starters;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_starters));

    std::vector<std::vector<hHalfedge>> holes;
    for (auto hh : border_starters)
    {
      holes.push_back(std::vector<hHalfedge>());
      for (auto b : CGAL::halfedges_around_face(hh, mesh))
        holes.back().push_back(b);
      if (holes.back().size() < threshold)
        holes.pop_back();
    }

    if (holes.size() < 2)
      return false;

    std::sort(holes.begin(), holes.end(), [](auto &left, auto &right)
              { return left.size() > right.size(); });

    Merge2Holes(holes[0], holes[1], mesh);
    return true;
  }

}

bool HoleMerge(std::string input_mesh, std::string output_mesh, int threshold)
{
  auto start = std::chrono::high_resolution_clock::now();
  Polyhedron mesh;
  if(!CGAL::IO::read_polygon_mesh(input_mesh, mesh))
  {
    std::cout << "Failed to read mesh: " << input_mesh << std::endl;
    return false;
  }
  std::cout << "Read Mesh: " << mesh.size_of_vertices() << ", " << mesh.size_of_facets() << std::endl;

  while (MergeLargestHoles(mesh, threshold)) {}

  if(!CGAL::IO::write_polygon_mesh(output_mesh, mesh))
  {
    std::cout << "Failed to output mesh: " << output_mesh << std::endl;
    return false;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start) << std::endl;
  return true;
}

#ifndef FOUND_PYBIND11

#endif