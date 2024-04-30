#include <iostream>
#include <queue>
#include <vector>
#include <unordered_map>
#include <CGAL/Weights/cotangent_weights.h>
#include <CGAL/Weights/mixed_voronoi_region_weights.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#define TINYCOLORMAP_WITH_EIGEN
#include <tinycolormap.hpp>
#include "Polyhedron.h"

using Kernel = CGAL::Simple_cartesian<double>;
using Polyhedron = TPolyhedron<CGAL::Polyhedron_items_with_id_3, Kernel>;
using hVertex = Polyhedron::Vertex_handle;

void ColorMeshByCurv(const std::string& input, const std::string& output)
{
    std::vector<Kernel::Point_3> vertices;
    std::vector<TTriangle<size_t>> indices;
    if(input.ends_with(".obj"))
    {
        LoadVFObj<CGAL::Simple_cartesian<double>, size_t>( input, vertices, indices );
    }
    else
    {
        LoadVFAssimp<CGAL::Simple_cartesian<double>, size_t>( input, vertices, indices);
    }
    std::vector<Eigen::Vector3f> colors;
    Polyhedron mesh;
    mesh.BuildFromVerticesFaces(vertices, indices);

    
    for(auto hv : CGAL::vertices(mesh))
    {
        auto xi = hv->point();
        Kernel::Vector_3 dx(0.0, 0.0, 0.0);
        double A = 0.0;
        for(auto nei : CGAL::halfedges_around_source(hv, mesh))
        {
            auto xj = nei->vertex()->point();
            auto xb = nei->next()->vertex()->point();
            auto xa = nei->opposite()->next()->vertex()->point();
            auto xij = xj - xi;
            double w = 0.5 * CGAL::Weights::cotangent_weight(xa, xj, xb, xi);
            dx += w * xij;
        }
        for(auto t : CGAL::faces_around_target(hv->halfedge(), mesh))
        {
            auto h = t->halfedge();
            while(h->vertex() != hv)
            {
                h = h->next();
            }
            auto q = h->vertex()->point();
            auto r = h->next()->vertex()->point();
            auto p = h->next()->next()->vertex()->point();
            A += CGAL::Weights::mixed_voronoi_area(p, q, r);
        }
        auto normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(hv, mesh);
        double curv = 0.0;
        if(CGAL::angle(normal, dx) == CGAL::Angle::ACUTE)
        {
            dx /= 2 * A;
            curv = 0.5 * std::sqrt(dx.squared_length());
        }
        if(curv < 0.5)
            curv = 0.0;
        Eigen::Vector3f color = tinycolormap::GetJetColor(curv).ConvertToEigen().cast<float>();
        colors.push_back(color);
    }
    WriteVCFAssimp<Kernel, size_t>(output, vertices, colors, indices);
}

void ColorMeshByCurvForce( const std::string& input, const std::string& output)
{
    std::vector<Kernel::Point_3> vertices;
    std::vector<TTriangle<size_t>> indices;
    if(input.ends_with(".obj"))
    {
        LoadVFObj<CGAL::Simple_cartesian<double>, size_t>( input, vertices, indices );
    }
    else
    {
        LoadVFAssimp<CGAL::Simple_cartesian<double>, size_t>( input, vertices, indices);
    }
    Polyhedron mesh;
    mesh.BuildFromVerticesFaces(vertices, indices);

    std::unordered_map<hVertex, Kernel::Vector_3> dx_map;
    std::unordered_map<hVertex, double> curv_map;
    for(auto hv : CGAL::vertices(mesh))
    {
        auto xi = hv->point();
        Kernel::Vector_3 dx(0.0, 0.0, 0.0);
        double A = 0.0;
        for(auto nei : CGAL::halfedges_around_source(hv, mesh))
        {
            auto xj = nei->vertex()->point();
            auto xb = nei->next()->vertex()->point();
            auto xa = nei->opposite()->next()->vertex()->point();
            auto xij = xj - xi;
            double w = 0.5 * CGAL::Weights::cotangent_weight(xa, xj, xb, xi);
            dx += w * xij;
        }
        for(auto t : CGAL::faces_around_target(hv->halfedge(), mesh))
        {
            auto h = t->halfedge();
            while(h->vertex() != hv)
            {
                h = h->next();
            }
            auto q = h->vertex()->point();
            auto r = h->next()->vertex()->point();
            auto p = h->next()->next()->vertex()->point();
            A += CGAL::Weights::mixed_voronoi_area(p, q, r);
        }
        auto normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(hv, mesh);
        double curv = 0.0;
        if(CGAL::angle(normal, dx) == CGAL::Angle::ACUTE)
        {
            dx /= 2 * A;
            curv = 0.5 * std::sqrt(dx.squared_length());
        }
        dx_map[hv] = dx;
        if(curv >= 1.0)
        {
            curv = 1.0;
        }
        else
        {
            curv = 0.0;
        }
        curv_map[hv] = curv;
    }

    std::unordered_map<hVertex, Kernel::Vector_3> force_map;
    for(auto hv : CGAL::vertices(mesh))
    {
        force_map[hv] = Kernel::Vector_3(0.0, 0.0, 0.0);
    }
    for(auto hv : CGAL::vertices(mesh))
    {
        std::unordered_set<hVertex> neighbors;
        std::queue<hVertex> q;
        q.push(hv);
        neighbors.insert(hv);
        while(!q.empty())
        {
            hVertex curr = q.front();
            q.pop();
            for(auto nei : CGAL::vertices_around_target(curr, mesh))
            {
                if(neighbors.count(nei) == 0 && CGAL::squared_distance(nei->point(), hv->point()) < 3.0 * 3.0)
                {
                    neighbors.insert(nei);
                    q.push(nei);
                }
            }
        }

        for(auto nei : neighbors)
        {
            if(curv_map[nei] > 1.0)
            {
                //continue;
            }
            auto d = hv->point() - nei->point();
            if(curv_map[hv] > curv_map[nei])            
            {
                force_map[nei] += (curv_map[hv] - curv_map[nei]) / d.squared_length() * d;
            }
        }
    }

    double max_f = 0.0;
    double avg_f = 0.0;
    for(auto& [hv, f] : force_map)
    {
        max_f = std::max(max_f, f.squared_length());
        avg_f += std::sqrt(f.squared_length());
    }
    max_f = std::sqrt(max_f);
    avg_f /= force_map.size();

    std::cout << "max force: " << max_f << ". avg force: " << avg_f << std::endl;
    std::vector<Eigen::Vector3f> colors;
    for(auto hv : CGAL::vertices(mesh))
    {
        Eigen::Vector3f color = tinycolormap::GetJetColor(std::sqrt(force_map[hv].squared_length()) / avg_f / 2.0).ConvertToEigen().cast<float>();
        colors.push_back(color);
    }

    //output forces
    std::ofstream ofs(R"(D:\dev\Ortho\OrthoMeshTools\test\ColorCurv\force.obj)");
    int idx = 1;
    for(auto hv : CGAL::vertices(mesh))
    {
        auto f = force_map[hv];
        if(f.squared_length() <= 1e-8)
        {
            continue;
        }
        Eigen::Vector3f color = tinycolormap::GetJetColor(std::sqrt(f.squared_length()) / avg_f / 2.0).ConvertToEigen().cast<float>();
        auto p = hv->point();
        auto p2 = hv->point() + f / std::sqrt(f.squared_length()) * 0.3;
        ofs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << ' ' << color.x() << ' ' << color.y() << ' ' << color.z() << '\n';
        ofs << "v " << p2.x() << ' ' << p2.y() << ' ' << p2.z() << " 1 1 1\n";
        ofs << "l " << idx << ' ' << idx + 1 << '\n';
        idx += 2;
    }

    WriteVCFAssimp<Kernel, size_t>(output, vertices, colors, indices);
}

int main(int argc, char* argv[])
{
    ColorMeshByCurv(R"(D:\dev\Ortho\OrthoMeshTools\test\ColorCurv\0.obj)", R"(D:\dev\Ortho\OrthoMeshTools\test\ColorCurv\0out.obj)");
    return 0;
}