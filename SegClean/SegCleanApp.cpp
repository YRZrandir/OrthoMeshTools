#include "SegClean.h"

struct Config
    {
        std::string input;
        std::string output;
        std::string labels;
        std::string outmesh;
        int size_threshold;
    };

Config LoadConfig(int argc, char *argv[])
{
    Config cfg;
    for (int i = 1; i < argc; i++)
    {
        if (std::strcmp(argv[i], "-i") == 0)
        {
            cfg.input = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-o") == 0)
        {
            cfg.output = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-l") == 0)
        {
            cfg.labels = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-m") == 0)
        {
            cfg.outmesh = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-t") == 0)
        {
            cfg.size_threshold = std::atoi(argv[i + 1]);
        }
    }
    return cfg;
}

int main(int argc, char *argv[])
{
    Config cfg = LoadConfig(argc, argv);
    SegCleanF(cfg.input, cfg.labels, cfg.output, cfg.size_threshold);
    return 0;
}