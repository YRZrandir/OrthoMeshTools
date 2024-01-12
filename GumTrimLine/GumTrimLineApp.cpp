#include "GumTrimLine.h"
#include <chrono>
#include <gflags/gflags.h>
                                                                                                                                                                                                                                                       
DEFINE_string(input_file, "", "file path of input mesh");
DEFINE_string(label_file, "", "");
DEFINE_string(frame_file, "", "");
DEFINE_string(output_file, "", "file path of output mesh");
DEFINE_int32(smooth, 10, "a non-nagetive integer that specifies the iteration number of trim line smoothing");
DEFINE_double(fix_factor, 0.0, "");

int main(int argc, char *argv[])
{
    gflags::SetUsageMessage("Extract gum trim line and output it as an .obj file.");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_input_file.empty() || FLAGS_label_file.empty() || FLAGS_output_file.empty() || FLAGS_smooth < 0)
    {
        std::cout << "Invalid paramters. Use --help for help." << std::endl;
        return -1;
    }

    try
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        GumTrimLine(FLAGS_input_file, FLAGS_label_file, FLAGS_frame_file, FLAGS_output_file, FLAGS_smooth, FLAGS_fix_factor);
        std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_time) << std::endl;
        std::cout << "===============================" << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}