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
    Chamber img_chamber(argv[1]);

    img_chamber.setTarget(argv[2]);

    img_chamber.setupAntennas(argv[3]);

    img_chamber.setFrequency(std::atof(argv[4]));

    img_chamber.setupProbes(argv[5]);

    img_chamber.calcDomainEzTot();

    img_chamber.calcDataEzTot();

    Eigen::MatrixXcd M_Ez_out;

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_inc
    );
    WriteMatrixToFile(
        "Ez_inc.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_sct
    );
    WriteMatrixToFile(
        "Ez_sct.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_tot
    );
    WriteMatrixToFile(
        "Ez_tot.txt",
        M_Ez_out
    );

    WriteMatrixToFile(
        "tri_pts.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        "tri_tri.txt",
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

    img_chamber.target.eps_r.WriteValsToFile(
        "eps_r.txt"
    );
    WriteMatrixToFile(
        "probe_xyz.txt",
        img_chamber.probe_points
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_inc_d
    );
    WriteMatrixToFile(
        "Ez_inc_d.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_sct_d
    );
    WriteMatrixToFile(
        "Ez_sct_d.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_tot_d
    );
    WriteMatrixToFile(
        "Ez_tot_d.txt",
        M_Ez_out
    );


    std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    std::cerr << " done!" << std::endl;
    return(0);
}
