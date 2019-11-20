#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>
//#include <filesystem>

#include "mom_driver.h"

#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675


int main (int argc, char **argv)
{
    gmsh::initialize();
    assert(argc == 7);
    Chamber img_chamber(argv[1]);

    img_chamber.setTarget(argv[2]);

    img_chamber.setupAntennas(argv[3]);

    img_chamber.setFrequency(std::atof(argv[4]));

    img_chamber.setupProbes(argv[5]);

    img_chamber.calcDomainEzTot();

    img_chamber.calcDataEzTot();

    std::string outdir(argv[6]);
    if(not (outdir.back()=='/' or outdir.back()=='\\'))
    {
        outdir += "/";
    }
    if(std::filesystem::exists(outdir))
    {
        std::cerr << "Directory: " << outdir << " already exists" << std::endl;
    }
    else
    {
        std::cerr << "Directory: " << outdir << " does not exist." << std::endl;
        auto dir_success{std::filesystem::create_directory(outdir)};

        if(not dir_success)
        {
            std::cerr << dir_success << std::endl;
            std::cerr << "Failed to make directory" << std::endl;
            assert(false);
        }
        else
        {
            std::cerr << "Now it does." << std::endl;
        }
    }

    Eigen::MatrixXcd M_Ez_out;

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_inc
    );
    WriteMatrixToFile(
        outdir + "Ez_inc.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_sct
    );
    WriteMatrixToFile(
        outdir + "Ez_sct.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_tot
    );
    WriteMatrixToFile(
        outdir + "Ez_tot.txt",
        M_Ez_out
    );

    WriteMatrixToFile(
        outdir + "tri_pts.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        outdir + "tri_tri.txt",
        img_chamber.mesh.tri
    );
    WriteVectorToFile(
        outdir + "tri_areas.txt",
        img_chamber.mesh.areas
    );
    WriteMatrixToFile(
        outdir + "tri_centroids.txt",
        img_chamber.mesh.centroids
    );

    img_chamber.target.eps_r.WriteValsToFile(
        outdir + "eps_r.txt"
    );
    WriteMatrixToFile(
        outdir + "probe_xyz.txt",
        img_chamber.probe_points
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_inc_d
    );
    WriteMatrixToFile(
        outdir + "Ez_inc_d.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_sct_d
    );
    WriteMatrixToFile(
        outdir + "Ez_sct_d.txt",
        M_Ez_out
    );

    FormFieldMatrix(
        M_Ez_out,
        img_chamber.Ez_tot_d
    );
    WriteMatrixToFile(
        outdir + "Ez_tot_d.txt",
        M_Ez_out
    );


    std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    std::cerr << " done!" << std::endl;
    return(0);
}
