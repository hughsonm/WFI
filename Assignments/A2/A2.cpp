#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/algorithm/string.hpp>

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
    //gmsh::initialize();

    std::ofstream writer;
    writer.open("myfile.txt", std::ofstream::out);
    for(int ll = 0; ll < 100; ll++)
    {
        for(int nn = 0; nn < 13; nn++)
        {
            std::complex<double> dd(
                ((double)std::rand())/RAND_MAX-0.5,
                ((double)std::rand())/RAND_MAX-0.5
            );
            writer << dd << "\t";
        }
        std::complex<double> dd(
            ((double)std::rand())/RAND_MAX-0.5,
            ((double)std::rand())/RAND_MAX-0.5
        );
        writer << dd << std::endl;
    }
    writer.close();

    std::ifstream reader;
    reader.open("myfile.txt", std::ifstream::in);
    while (!reader.eof())
    {
        std::complex<double> id;
        reader >> id;
        std::cout << id << std::endl;
    }


    /*Chamber img_chamber(argv[1]);
    img_chamber.setupAntennas(argv[3]);
    img_chamber.setFrequency(std::atof(argv[4]));
    img_chamber.setupProbes(argv[5]);*/
    reader.close();
    //gmsh::finalize();
    return(0);
}
