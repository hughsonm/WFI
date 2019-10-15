#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>


#include "../A1/mom_driver.h"

int main(int argc, char** argv)
{
    gmsh::initialize();
    assert(argc == 6);
    std::cerr << "Readimg mesh" << std::endl;
    std::string mesh_filename(argv[1]);
    std::cerr << "Reading Antennas" << std::endl;
    std::string antennas_filename(argv[2]);
    std::cerr << "Setting frequency" << argv[3] << std::endl;
    double frequency = std::stod(argv[3]);
    std::cerr << "Reading Probes" << std::endl;
    std::string probes_filename(argv[4]);
    std::cerr << "Reading data" << std::endl;
    std::string data_filename(argv[5]);

    Chamber img_chamber(mesh_filename);
    std::cerr << "Antennas" << std::endl;
    img_chamber.setupAntennas(antennas_filename);
    std::cerr << "Freq" << std::endl;
    img_chamber.setFrequency(frequency);
    std::cerr << "Probes" << std::endl;
    img_chamber.setupProbes(probes_filename);
    std::cerr << "Measurement" << std::endl;
    img_chamber.readMeasuredData(data_filename);

    WriteMatrixToFile(
        "Points.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        "Tri.txt",
        img_chamber.mesh.tri
    );

    img_chamber.buildDataGreen();
    std::cerr << "Gon try to build Annihilator!" << std::endl;
    img_chamber.buildAnnihilator();

    gmsh::finalize();
    return(0);
}
