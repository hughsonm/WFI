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
    assert(argc == 6);
    std::string mesh_filename(argv[1]);
    std::string antennas_filename(argv[2]);
    double frequency = std::stod(argv[3]);
    std::string probes_filename(argv[4]);
    std::string data_filename(argv[5]);

    Chamber img_chamber(mesh_filename);
    std::cerr << "Antennas" << std::endl;
    img_chamber.setupAntennas(antennas_filename);
    std::cerr << "Freq" << std::endl;
    img_chamber.setFrequency(frequency);
    std::cerr << "Probes" << std::endl;
    img_chamber.setupProbes(probes_filename);
    std::cerr << "Measurement" << std::endl;
    img_chamber.readMeasuredData(data_filename);
    std::cerr << "CalcEzTot" << std::endl;
    img_chamber.calcDataEzTot();


    img_chamber.buildDataGreen();
    std::cerr << "Gon try to build Annihilator!" << std::endl;
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
