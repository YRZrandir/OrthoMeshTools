#include "ColorMeshByLabel.h"
#include <iostream>
int main(int argc, char* argv[])
{
    std::string input;
    std::string labels;
    std::string output;

    for(int i = 1; i < argc; i++)
    {
        if(std::strcmp(argv[i], "-i") == 0)
        {
            input = std::string(argv[i+1]);
        }
        else if(std::strcmp(argv[i], "-l") == 0)
        {
            labels = std::string(argv[i+1]);
        }
        else if(std::strcmp(argv[i], "-o") == 0)
        {
            output = std::string(argv[i+1]);
        }
    }

    std::cout << "input: " << input << std::endl;
    std::cout << "labels: " << labels << std::endl;
    std::cout << "output: " << output << std::endl;

    if(ColorMeshByLabel(input, labels, output))
    {
        return 0;
    }
    return -1;
}