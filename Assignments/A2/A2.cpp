#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>


#include "../A1/mom_driver.h"

bool char_is_comma(char c)
{
    return(c==',');
}

bool char_is_plusminus(char c)
{
    return(c=='+' || c=='-');
}

int main(int argc, char** argv)
{
    gmsh::initialize();
    assert(argc == 6);
    std::string mesh_filename(argv[1]);
    std::string antennas_filename(argv[2]);
    double frequency = std::stod(argv[3]);
    std::string probes_filename(argv[4]);
    std::string data_filename(argv[5]);

    Chamber img_chamber(mesh_filename);
    img_chamber.setupAntennas(antennas_filename);
    img_chamber.setFrequency(frequency);
    img_chamber.setupProbes(probes_filename);
    img_chamber.readMeasuredData(data_filename);

    img_chamber.buildDataGreen();
    img_chamber.buildAnnihilator();


    Eigen::MatrixXcd Ez_tot_d_meas;
    ReadMatrixFromFile(
        "Ez_tot_d.txt",
        Ez_tot_d_meas
    );
    std::cerr << "Here is our data!" << std::endl;
    std::cerr << Ez_tot_d_meas << std::endl;

    std::cerr << "Measured data is " << Ez_tot_d_meas.rows() << " by " << Ez_tot_d_meas.cols() << std::endl;
    std::cerr << "We have " << img_chamber.probe_points.rows() << " probes" << std::endl;
    std::cerr << "We have " << img_chamber.antennas.size() << " antennas" << std::endl;

    assert(Ez_tot_d_meas.rows()==img_chamber.probe_points.rows());
    assert(Ez_tot_d_meas.cols()==img_chamber.antennas.size());

    gmsh::finalize();
    return(0);
}
