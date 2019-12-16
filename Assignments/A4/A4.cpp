#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <gmsh.h>
#include "../A1/mom_driver.h"

double rand_double(double low, double high){
    double rr = ((double)std::rand())/((double)RAND_MAX);
    return(low+(high-low)*rr);
}

bool make_outdir(std::string& dirname){
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
    assert(argc==8);
    std::string meshfile{argv[1]};
    std::string antennafile{argv[2]};
    std::string probefile{argv[3]};
    std::string tx2rxfile{argv[4]};
    std::string bimfreqfile{argv[5]};
    std::string total_data_prefix{argv[6]};
    std::string output_directory_name{argv[7]};


    std::cout << "Mesh File     : " << meshfile << "\n";

    Chamber img_chamber(meshfile);
    img_chamber.setFrequencies(bimfreqfile);
    img_chamber.setupAntennas(antennafile);
    img_chamber.setupProbes(probefile);
    img_chamber.setupTx2RxMap(tx2rxfile);
    img_chamber.readMeasuredData(total_data_prefix,0.0001);

    DBIMInversion inv_dbim = img_chamber.distortedBornIterativeMethod();
    // BIMInversion inv_bim = img_chamber.bornIterativeMethod();


    if(make_outdir(output_directory_name)){
        inv_dbim.WriteResults(output_directory_name);
        // inv_bim.WriteResults(output_directory_name);
    }

    gmsh::finalize();
    return(0);
}
