#include "FakeToothRoot.h"
int main(int argc, char* argv[])
{
    std::string input_path;
    std::string output_path;
    std::string frame_path;
    std::string label_path;
    int fair = 1; // larger == more fair

    for(int i = 1; i < argc; i++)
    {
        if(std::strcmp(argv[i], "-i") == 0)
        {
            input_path = std::string(argv[i+1]);
        }
        if(std::strcmp(argv[i], "-o") == 0)
        {
            output_path = std::string(argv[i+1]);
        }
        if(std::strcmp(argv[i], "-f") == 0)
        {
            frame_path = std::string(argv[i+1]);
        }
        if(std::strcmp(argv[i], "-l") == 0)
        {
            label_path = std::string(argv[i+1]);
        }
        if(std::strcmp(argv[i], "-s") == 0)
        {
            fair = std::atoi(argv[i + 1]);
        }
    }
    FakeToothRoot(input_path, output_path, frame_path, label_path, fair);
    return 0;
}
