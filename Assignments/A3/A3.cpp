#include <iostream>
#include <vector>
#include <filesystem>
#include <gmsh.h>
#include "../A1/mom_driver.h"

double rand_double(double low, double high){
    double rr = ((double)std::rand())/((double)RAND_MAX);
    return(low+(high-low)*rr);
}

bool make_outdir(std::string dirname){
    bool success{false};
    if(not (dirname.back()=='/' or dirname.back()=='\\'))
    {
        dirname += "/";
    }
    if(std::filesystem::exists(dirname))
    {
        std::cerr << "Directory: " << dirname << " already exists" << std::endl;
        success = true;
    }
    else
    {
        std::cerr << "Directory: " << dirname << " does not exist." << std::endl;
        success = std::filesystem::create_directory(dirname);
        if(not success)
        {
            std::cerr << success << std::endl;
            std::cerr << "Failed to make directory" << std::endl;
            assert(false);
        }
        else
        {
            success = true;
            std::cerr << "Now it does." << std::endl;
        }
    }
    return(success);
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
    std::string p3_freqfile{argv[2]};
    std::string p4_freqfile{argv[3]};
    std::string antennafile{argv[4]};
    std::string probefile{argv[5]};
    std::string datatotalprefix{argv[6]};
    std::string outdir{argv[7]};

    std::cout << "Running A3...\n";
    std::cout << "Mesh File         : " << meshfile << "\n";
    std::cout << "P4 Freq File      : " << p4_freqfile << "\n";
    std::cout << "Antenna File      : " << antennafile << "\n";
    std::cout << "Probe File        : " << probefile << "\n";
    std::cout << "Total Data Prefix : " << datatotalprefix << "\n";
    std::cout << "Output Directory  : " << outdir << "\n";

    Chamber img_chamber(meshfile);
    img_chamber.setFrequencies(p3_freqfile);
    img_chamber.setupAntennas(antennafile);
    img_chamber.setupProbes(probefile);
    img_chamber.readMeasuredData(datatotalprefix,0.0);

    Eigen::MatrixXd kpts;
    Eigen::VectorXcd kvals;
    Eigen::VectorXcd chi;

    img_chamber.A3P3(
        kpts,
        kvals,
        chi
    );


    make_outdir(outdir);

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

    // img_chamber.A3P4(p4_freqfile);


    gmsh::finalize();
    return(0);
}
