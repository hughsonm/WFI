#include <iostream>
#include <vector>
#include <filesystem>
#include <gmsh.h>
#include "../A1/mom_driver.h"

double rand_double(double low, double high){
    double rr = ((double)std::rand())/((double)RAND_MAX);
    return(low+(high-low)*rr);
}

int main(int argc, char** argv){
    gmsh::initialize();
    std::cout << "Hello from A3\n";

    // Pass in arguments as a single directory which must contain certain files.
    //  a list of frequencies
    //  a list of transmitter positions
    //  a list of probe positions
    //  a list of scattered-field data at the probe locations
    assert(argc==8);
    std::string meshfile{argv[1]};
    std::string p4_freqfile{argv[2]};
    double p3_freq{std::stod(argv[3])};
    std::string antennafile{argv[4]};
    std::string probefile{argv[5]};
    std::string datatotalfile{argv[6]};
    std::string outdir{argv[7]};

    std::cout << "Running A3...\n";
    std::cout << "Mesh File        : " << meshfile << "\n";
    std::cout << "P4 Freq File     : " << p4_freqfile << "\n";
    std::cout << "Antenna File     : " << antennafile << "\n";
    std::cout << "Probe File       : " << probefile << "\n";
    std::cout << "Total Field File : " << datatotalfile << "\n";
    std::cout << "Output Directory : " << outdir << "\n";

    Chamber img_chamber(meshfile);
    img_chamber.setFrequency(p3_freq);
    img_chamber.setupAntennas(antennafile);
    img_chamber.setupProbes(probefile);
    img_chamber.readMeasuredData(datatotalfile,0.0);
    
    Eigen::MatrixXd kpts;
    Eigen::VectorXcd kvals;
    Eigen::VectorXcd chi;

    img_chamber.A3P3(
        kpts,
        kvals,
        chi
    );



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

    WriteVectorToFile(
        outdir + "chi.txt",
        chi
    );
    WriteMatrixToFile(
        outdir + "pts.txt",
        img_chamber.mesh.points
    );
    WriteMatrixToFile(
        outdir + "tri.txt",
        img_chamber.mesh.tri
    );

    WriteMatrixToFile(
        outdir + "kpts.txt",
        kpts
    );

    WriteVectorToFile(
        outdir + "kvals.txt",
        kvals
    );


    gmsh::finalize();
    return(0);
}
