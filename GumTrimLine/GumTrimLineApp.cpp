#include "GumTrimLine.h"
#include <chrono>

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();
    std::string input_file = "";
    std::string label_file = "";
    std::string frame_file = "";
    std::string output_file = "";
    int smooth = 20;
    double fix_factor = 0.0;
    for (int i = 1; i < argc; i++)
    {
        if (std::strcmp(argv[i], "-i") == 0)
        {
            input_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-o") == 0)
        {
            output_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-l") == 0)
        {
            label_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-s") == 0)
        {
            smooth = std::atoi(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-a") == 0)
        {
            fix_factor = std::atof(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-f") == 0)
        {
            frame_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-h") == 0)
        {
            std::cout << "Extract gum trim line and output it as an .obj file.\n"
                         "-i : file path of input mesh\n"
                         "-o : file path of output mesh\n"
                         "-l : file path of label file.\n"
                         "-s : a non-nagetive integer that specifies the iteration number of trim line smoothing."
                      << std::endl;
            return 0;
        }
    }
    if (input_file.empty() || output_file.empty() || label_file.empty() || smooth < 0)
    {
        std::cout << "Invalid paramters. Use -h for help." << std::endl;
        return -1;
    }
    try
    {
        GumTrimLine(input_file, label_file, frame_file, output_file, smooth, fix_factor);
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