#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>

#include "mom_driver.h"

#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675


int main (int argc, char **argv)
{
    gmsh::initialize();
    std::cerr << "Chamber img_chamber(argv[1]);" << std::endl;
    Chamber img_chamber(argv[1]);

    std::cerr << "img_chamber.addTarget(argv[2]);" << std::endl;
    img_chamber.addTarget(argv[2]);

    std::cerr << "img_chamber.setupAntennas(argv[3]);" << std::endl;
    img_chamber.setupAntennas(argv[3]);

    std::cerr << "img_chamber.setFrequency(std::atof(argv[4]));" << std::endl;
    img_chamber.setFrequency(std::atof(argv[4]));

    std::cerr << "Eigen::MatrixXcd Eztot;" << std::endl;
    Eigen::MatrixXcd Eztot;

    std::cerr << "img_chamber.getDomainEzTot(Eztot);" << std::endl;
    img_chamber.getDomainEzTot(Eztot);
	
	WriteMatrixToFile(
		"Eztot.txt",
		Eztot
	);

    std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    std::cerr << " done!" << std::endl;
    return(0);
}
