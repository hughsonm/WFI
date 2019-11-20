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
    assert(argc==7);
    std::string meshfile{argv[1]};
    std::string freqfile{argv[2]};
    std::string antennafile{argv[3]};
    std::string probefile{argv[4]};
    std::string datatotalfile{argv[5]};
    std::string outdir{argv[6]};

    std::cout << "Running A3...\n";
    std::cout << "Mesh File        : " << meshfile << "\n";
    std::cout << "Freq File        : " << freqfile << "\n";
    std::cout << "Antenna File     : " << antennafile << "\n";
    std::cout << "Probe File       : " << probefile << "\n";
    std::cout << "Total Field File : " << datatotalfile << "\n";
    std::cout << "Output Directory : " << outdir << "\n";

    Chamber img_chamber(meshfile);
    img_chamber.setupAntennas(antennafile);
    img_chamber.setupProbes(probefile);
    img_chamber.readMeasuredData(datatotalfile,10);
    gmsh::finalize();
    return(0);
}
