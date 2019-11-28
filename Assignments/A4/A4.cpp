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
    std::string meshfile{argv[1]};
    std::string tx2rxfile{argv[2]};
    std::cout << "Mesh File     : " << meshfile << "\n";

    Chamber img_chamber(meshfile);
    img_chamber.setupTx2RxMap(tx2rxfile);
    gmsh::finalize();
    return(0);
}