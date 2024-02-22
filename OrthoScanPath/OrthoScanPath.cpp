#include "OrthoScanPath.h"
#include <argparse/argparse.hpp>

std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> OrthoScanPath(
    std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> input_path,
    const CBCTRegis<double>& cbct_regis,
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> cbct_frames )
{
    auto paths = input_path;
    for(auto& [_, T] : cbct_frames)
        T.translation() *= -1;
    for(auto& step : paths)
    {
        for(auto& [label, T] : step)
        {
            T = Eigen::Affine3d(cbct_regis.CBCT_to_IOS(label).matrix() * T.matrix() * cbct_frames.at(label).matrix() * cbct_regis.IOS_to_CBCT(label).matrix());
        }
    }
    return paths;
}
