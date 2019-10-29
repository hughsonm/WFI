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
    assert(argc == 8);
    std::cerr << "Readimg mesh: " << argv[1] << std::endl;
    std::string mesh_filename(argv[1]);
    std::cerr << "Reading Antennas: " << argv[2] << std::endl;
    std::string antennas_filename(argv[2]);
    std::cerr << "Reading frequency: " << argv[3] << std::endl;
    double frequency = std::stod(argv[3]);
    std::cerr << "Reading Probes: " << argv[4] << std::endl;
    std::string probes_filename(argv[4]);
    std::cerr << "Reading total field data: " << argv[5] << std::endl;
    std::string data_filename(argv[5]);
    std::cerr << "Reading noise: " << argv[6] << " percent" << std::endl;
    double noise_pct{std::stod(argv[6])};


    std::string outdir(argv[7]);
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

    Chamber img_chamber(mesh_filename);
    std::cerr << "Antennas" << std::endl;
    img_chamber.setupAntennas(antennas_filename);
    std::cerr << "Freq" << std::endl;
    img_chamber.setFrequency(frequency);
    std::cerr << "Probes" << std::endl;
    img_chamber.setupProbes(probes_filename);
    std::cerr << "Measurement" << std::endl;
    img_chamber.readMeasuredData(data_filename,noise_pct);



    // Q3
    Eigen::MatrixXcd q3_w;
    Eigen::VectorXcd q3_X;
    img_chamber.A2Q3(
        q3_w,
        q3_X
    );
    WriteMatrixToFile(
        outdir + "q3_w.txt",
        q3_w
    );
    WriteVectorToFile(
        outdir + "q3_X.txt",
        q3_X
    );



    // Q5
    img_chamber.buildDataGreen();
    std::cerr << "Gon try to build Annihilator!" << std::endl;

    Eigen::MatrixXcd q5_w;
    Eigen::MatrixXcd q5_u;
    Eigen::VectorXcd q5_X;
    img_chamber.A2Q5(
        q5_w,
        q5_u,
        q5_X
    );
    WriteMatrixToFile(
        outdir + "q5_w.txt",
        q5_w
    );
    WriteMatrixToFile(
        outdir + "q5_u.txt",
        q5_u
    );
    WriteVectorToFile(
        outdir + "q5_X.txt",
        q5_X
    );


    WriteMatrixToFile(
        outdir + "Points.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        outdir + "Tri.txt",
        img_chamber.mesh.tri
    );
    WriteMatrixToFile(
      outdir + "Gb.txt",
      img_chamber.G_b_domain
    );
    WriteMatrixToFile(
      outdir + "Gbd.txt",
      img_chamber.G_b_data
    );
    gmsh::finalize();
    return(0);
}
