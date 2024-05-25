#include "HoleMerge.h"
#include <iostream>
#include <CGAL/version.h>

int main(int argc, char *argv[])
{
  std::cout << "CGAL " << CGAL_STR(CGAL_VERSION) << std::endl;
  std::string input_mesh;
  std::string output_mesh;
  int threshold = 0;
  for (int i = 1; i < argc; i++)
  {
    if (std::strcmp(argv[i], "-i") == 0)
    {
      input_mesh = std::string(argv[i + 1]);
    }
    else if (std::strcmp(argv[i], "-o") == 0)
    {
      output_mesh = std::string(argv[i + 1]);
    }
    else if (std::strcmp(argv[i], "-t") == 0)
    {
      threshold = std::atoi(argv[i + 1]);
    }
    else if (std::strcmp(argv[i], "-h") == 0)
    {
      std::cout << "Merge multiple holes on the mesh into one hole.\n"
      "-i : file path of input mesh\n"
      "-o : file path of output mesh\n"
      "-t : an integer threshold value. Holes with edge number smaller than this value will NOT be processed." << std::endl;
      return 0;
    }
  }
  if(input_mesh.empty() || output_mesh.empty() || threshold < 0)
  {
    std::cout << "Invalid paramters. Use -h for help." << std::endl;
    return -1;
  }
  if(!HoleMerge(input_mesh, output_mesh, threshold))
  {
    return -1;
  }
  return 0;
}