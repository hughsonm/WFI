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
    BIMInversion inv_bim = img_chamber.bornIterativeMethod();


    if(make_outdir(output_directory_name)){
        inv_bim.imaging_mesh.WriteMeshToFile(output_directory_name);
        auto iter_count{0};
        for(auto& step : inv_bim.steps){
            step.chi.WriteValsToFile(
                output_directory_name +
                "chi_iter_" + std::to_string(iter_count) +
                ".txt"
            );
            auto step_freq_count{0};
            for(auto& vec_of_tot_fields : step.Utot){
                auto field_count{0};
                for(auto& tot_field:vec_of_tot_fields){
                    tot_field.WriteValsToFile(
                        output_directory_name +
                        "Ez_tot_" +
                        "iter_" + std::to_string(iter_count) +
                        "freq_" + std::to_string(step_freq_count) +
                        "tx_"   + std::to_string(field_count) +
                        ".txt"
                    );
                    field_count++;
                }
                step_freq_count++;
            }
            iter_count++;
        }
        auto freq_count{0};
        for(auto& vec_of_sct_fields : inv_bim.Ez_sct_meas){
            Eigen::MatrixXcd meas_mat;
            FormFieldMatrix(
                meas_mat,
                vec_of_sct_fields
            );
            WriteMatrixToFile(
                output_directory_name +
                "Ez_sct_meas_" +
                "freq_" +
                std::to_string(freq_count) +
                ".txt",
                meas_mat
            );
            freq_count++;
        }
        Eigen::VectorXcd Fs_vec(inv_bim.steps.size());
        for(auto ii{0}; ii<Fs_vec.size(); ++ii){
            Fs_vec(ii) = inv_bim.steps[ii].Fs;
        }
        WriteVectorToFile(
            output_directory_name + "Fs.txt",
            Fs_vec
        );

    }

    gmsh::finalize();
    return(0);
}
