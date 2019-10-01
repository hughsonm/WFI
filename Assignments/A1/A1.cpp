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

    img_chamber.setupProbes(argv[5]);

    std::cerr << "img_chamber.calcDomainEzTot(Eztot);" << std::endl;
    img_chamber.calcDomainEzTot();

    img_chamber.calcDataEzTot();

    WriteMatrixToFile(
        "Ez_inc.txt",
        img_chamber.Ez_inc
    );
    WriteMatrixToFile(
        "Ez_sct.txt",
        img_chamber.Ez_sct
    );
	WriteMatrixToFile(
		"Ez_tot.txt",
		img_chamber.Ez_tot
	);

    WriteMatrixToFile(
        "Test_pts.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        "Test_tri.txt",
        img_chamber.mesh.tri
    );
    WriteVectorToFile(
        "tri_areas.txt",
        img_chamber.mesh.areas
    );
    WriteMatrixToFile(
        "tri_centroids.txt",
        img_chamber.mesh.centroids
    );

    WriteVectorToFile(
        "k2_fgd.txt",
        img_chamber.k2_f
    );

    WriteMatrixToFile(
        "probe_xyz.txt",
        img_chamber.probe_points
    );

    WriteMatrixToFile(
        "Ez_inc_d.txt",
        img_chamber.Ez_inc_d
    );

    WriteMatrixToFile(
        "Ez_sct_d.txt",
        img_chamber.Ez_sct_d
    );

    WriteMatrixToFile(
        "Ez_tot_d.txt",
        img_chamber.Ez_tot_d
    );

    std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    std::cerr << " done!" << std::endl;
    return(0);
}
