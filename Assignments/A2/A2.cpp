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

    Eigen::MatrixXcd q5_w_tikh;
    Eigen::MatrixXcd q5_u_tikh;
    Eigen::VectorXcd q5_X_tikh;
    std::vector<std::vector<Eigen::Vector2d> > q5_curves;
    img_chamber.A2Q5(
        q5_w_tikh,
        q5_u_tikh,
        q5_X_tikh,
        q5_curves,
        true
    );
    for(auto tt{0}; tt<q5_curves.size(); ++tt)
    {
        Eigen::MatrixXd curve_mat(2,q5_curves[tt].size());
        for(auto pp{0}; pp<q5_curves[tt].size(); ++pp)
        {
            curve_mat.col(pp) = q5_curves[tt][pp];
        }
        WriteMatrixToFile(
            outdir + "q5_curve_" + std::to_string(tt) + ".txt",
            curve_mat
        );
    }
    WriteMatrixToFile(
        outdir + "q5_w_tikh.txt",
        q5_w_tikh
    );
    WriteMatrixToFile(
        outdir + "q5_u_tikh.txt",
        q5_u_tikh
    );
    WriteVectorToFile(
        outdir + "q5_X_tikh.txt",
        q5_X_tikh
    );

    Eigen::MatrixXcd q5_w_unreg;
    Eigen::MatrixXcd q5_u_unreg;
    Eigen::VectorXcd q5_X_unreg;
    img_chamber.A2Q5(
        q5_w_unreg,
        q5_u_unreg,
        q5_X_unreg,
        q5_curves,
        false
    );
    WriteMatrixToFile(
        outdir + "q5_w_unreg.txt",
        q5_w_unreg
    );
    WriteMatrixToFile(
        outdir + "q5_u_unreg.txt",
        q5_u_unreg
    );
    WriteVectorToFile(
        outdir + "q5_X_unreg.txt",
        q5_X_unreg
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
