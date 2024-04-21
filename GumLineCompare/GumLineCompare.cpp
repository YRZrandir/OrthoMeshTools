#include <iostream>
#include <fstream>
#include <vector>
#include <argparse/argparse.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_traits.h>

using Kernel = CGAL::Simple_cartesian<double>;

std::vector<Kernel::Point_3> LoadGumLine(const std::string& path)
{
    std::vector<Kernel::Point_3> points;
    std::ifstream ifs(path);
    std::string line;
    if(ifs.is_open())
    {
        while(std::getline(ifs, line))
        {
            if(line.starts_with("v "))
            {
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                std::stringstream ss(line);
                ss.seekg(2);
                ss >> x >> y >> z;
                if(ss.fail())
                {
                    throw std::runtime_error("file format error: " + path);
                }
                points.emplace_back(x, y, z);
            }
        }
    }
    else
    {
        throw std::runtime_error("Cannot open file: " + path);
    }
    return points;
}

double CalcDist(const std::vector<Kernel::Point_3>& source, const std::vector<Kernel::Point_3>& target)
{
    double dist = 0.0;
    std::vector<Kernel::Segment_3> target_segs;
    for(size_t i = 0; i < target.size(); i++)
    {
        target_segs.emplace_back(target[i], target[(i + 1) % target.size()]);
    }
    std::vector<Kernel::Segment_3> source_segs;
    for(size_t i = 0; i < source.size(); i++)
    {
        source_segs.emplace_back(source[i], source[(i + 1) % source.size()]);
    }

    using AABBPrimitive = CGAL::AABB_segment_primitive<Kernel, std::vector<Kernel::Segment_3>::iterator>;
    using AABBTraits = CGAL::AABB_traits<Kernel, AABBPrimitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;
    AABBTree aabb(target_segs.begin(), target_segs.end());
    for(size_t i = 0; i < source_segs.size(); i++)
    {
        auto curr = source_segs[i];
        auto next = source_segs[(i + 1) % source_segs.size()];
        auto proj = aabb.closest_point(next.start());
        double d = std::sqrt(CGAL::squared_distance(proj, next.start()));
        dist += d * (std::sqrt(curr.squared_length()) + std::sqrt(next.squared_length()));
    }
    return dist;
}

int main(int argc, char* argv[])
{
    argparse::ArgumentParser parser("GumLineCompare");
    parser.add_description("Compare 2 gumline. The input files have to be .obj format, and the vertices should be in correct order.");
    parser.add_argument("-s").help("source gumline.").required();
    parser.add_argument("-t").help("target gumline.").required();
    parser.parse_args(argc, argv);
    try
    {
        auto src_line = LoadGumLine(parser.get("-s"));
        auto tgt_line = LoadGumLine(parser.get("-t"));
        std::cout << "Dist = " << CalcDist(src_line, tgt_line) << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}