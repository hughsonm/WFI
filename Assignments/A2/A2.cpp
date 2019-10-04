#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/algorithm/string.hpp>

#include "../A1/mom_driver.h"

bool char_is_comma(char c)
{
    return(c==',');
}

bool char_is_plusminus(char c)
{
    return(c=='+' || c=='-');
}

int main(int argc, char** argv)
{
    gmsh::initialize();
    assert(argc == 5);
    std::string mesh_filename(argv[1]);
    std::string antennas_filename(argv[2]);
    double frequency = std::stod(argv[3]);
    std::string probes_filename(argv[4]);

    Chamber img_chamber(mesh_filename);
    img_chamber.setupAntennas(antennas_filename);
    img_chamber.setFrequency(frequency);
    img_chamber.setupProbes(probes_filename);
    gmsh::finalize();
    return(0);
}
