#ifndef ORTHO_SCAN_PATH_H
#define ORTHO_SCAN_PATH_H
#include "Ortho.h"
std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> OrthoScanPath(
    std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> input_path,
    const CBCTRegis<double>& cbct_regis,
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> cbct_frames );
#endif