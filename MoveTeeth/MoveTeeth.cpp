#include "Polyhedron.h"
#include <argparse/argparse.hpp>

int main(int argc, char* argv[])
{
    argparse::ArgumentParser args("MoveTeeth");
    args.add_argument("-w").help("folder of teeth meshes").required();
    args.add_argument("-f").help("ortho scan frame").required();
    args.add_argument("-p").help("path file").required();
    args.add_argument("-r").help("cbct registration").required();
    args.add_argument("-t").help("cbct teeth").required();
    
    return 0;
}