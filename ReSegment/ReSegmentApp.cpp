#include <iostream>
#include <format>
#include <fstream>
#include <filesystem>
#include <CGAL/version.h>
#include "ReSegment.h"

int main(int argc, char *argv[])
{
    // std::cout << "Not implemented. Please use python interface." << std::endl;
    // return 0;
    std::cout << std::format("CGAL: {}", CGAL_STR(CGAL_VERSION)) << std::endl;
    std::filesystem::current_path("../../test/reseg2/");
    // std::ifstream ifs("../../test/case1/seg.json");
    // nlohmann::json json = nlohmann::json::parse(ifs);
    // ifs.close();
    // std::vector<std::vector<double>> splitlines = json["ctrl"];
    // std::vector<int> labels = json["labels"];

    std::vector<std::vector<std::vector<double>>> splitlines;
    std::vector<int> splitline_labels;
    splitlines.resize(1);
    splitline_labels.push_back(35);

    {
        std::ifstream ifs("points.xyz");
        std::vector<double> xyz;
        while (ifs)
        {
            double v = 0.0;
            ifs >> v;
            xyz.push_back(v);
        }
        splitlines[0].push_back(xyz);
    }

    ReSegmentLabels("oral_scan_U.obj", splitlines, splitline_labels, "out.json", 0.2, 5, false);
    return 0;
}