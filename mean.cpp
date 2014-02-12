#include "mesh_1d.h"

#include <iostream>

int main()
{
    Mesh_1D<double> mesh;

    mesh.load_data("output.dat");

    mesh.dump();

    return 0;
}
