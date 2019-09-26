#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>

#include "mom_driver.h"

#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675


int main (int argc, char **argv)
{
    gmsh::initialize();
    std::cerr << "Chamber img_chamber(argv[1]);" << std::endl;
    Chamber img_chamber(argv[1]);

    std::cerr << "img_chamber.addTarget(argv[2]);" << std::endl;
    img_chamber.addTarget(argv[2]);

    std::cerr << "img_chamber.setupAntennas(argv[3]);" << std::endl;
    img_chamber.setupAntennas(argv[3]);

    std::cerr << "img_chamber.setFrequency(std::atof(argv[4]));" << std::endl;
    img_chamber.setFrequency(std::atof(argv[4]));

    std::cerr << "Eigen::MatrixXcd Eztot;" << std::endl;
    Eigen::MatrixXcd Eztot;

    std::cerr << "img_chamber.getDomainEzTot(Eztot);" << std::endl;
    img_chamber.getDomainEzTot(Eztot);


    // assert(argc==4);
    // std::cerr << "\nInitializing gmsh...";
    // gmsh::initialize();
    // std::cerr << " done!" << std::endl;
    //
    // double frequency = std::atof(argv[3]);
    // Eigen::MatrixXd tri;
    // Eigen::MatrixXd points;
    // std::cerr << "Building triangulation...";
    // BuildTriangulation(
    //     argv[1],
    //     tri,
    //     points,
    //     true
    // );
    // std::cerr << " done!" << std::endl;
    //
    // Eigen::VectorXcd eps_r;
    //
    // AssignRelativeConstitutives(
    //     eps_r,
    //     tri,
    //     "Tag2EpsRel.txt"
    // );
    //
    // std::cerr << "\nCalculating triangle areas...";
    // Eigen::VectorXd areas;
    // CalculateTriAreas(
    //     areas,
    //     tri,
    //     points
    // );
    //
    // std::cerr << " done!" << std::endl;
    //
    // std::cerr << "\nCalculating triangle centroids...";
    // Eigen::MatrixXd centroids;
    // CalculateTriCentroids(
    //     centroids,
    //     tri,
    //     points
    // );
    //
    // std::cerr << " done!" << std::endl;
    //
    // std::cerr << "\nReading antenna file...";
    // Eigen::MatrixXd tx_locations;
    // Eigen::VectorXcd tx_coefficients;
    // ReadAntennaFile(
    //     tx_locations,
    //     tx_coefficients,
    //     argv[2]
    // );
    // std::cerr << " done!" << std::endl;
    //
    // std::cerr << "\nCalculating incident fields...";
    // Eigen::MatrixXcd Ez_inc_all(centroids.rows(),tx_locations.rows());
    // Eigen::VectorXcd Ez_inc_one;
    // for(int itx = 0; itx < tx_locations.rows(); itx++)
    // {
    //     EvaluateIncidentField(
    //         Ez_inc_one,
    //         tx_locations.row(itx),
    //         tx_coefficients(itx),
    //         frequency,
    //         centroids
    //     );
    //     std::cerr << itx << " ";
    //     Ez_inc_all.col(itx) = Ez_inc_one;
    // }
    // std::cerr << " done!" << std::endl;
    //
    // double omega = 2*M_PI*frequency;
    // double k2_b = omega*omega/CNAUGHT/CNAUGHT;
    //
    // Eigen::VectorXcd k2_f(eps_r.size());
    // k2_f = eps_r.array()*k2_b;
    //
    // std::cerr << "\nCalculating background Green function...";
    // Eigen::MatrixXcd G_b;
    // BuildDomainGreen(
    //     G_b,
    //     centroids,
    //     areas,
    //     k2_b
    // );
    // std::cerr << " done!" << std::endl;
    // std::cerr << "\nFilling Chi's diagonal with k2...";
    // Eigen::MatrixXcd Chi(areas.size(),areas.size());
    // Chi.setZero();
    // for(int dd = 0; dd < Chi.rows(); dd++)
    // {
    //     Chi(dd,dd) = k2_f(dd) - k2_b;
    // }
    // std::cerr << " done!" << std::endl;
    //
    // std::cerr << "\nAllocating space for L...";
    // Eigen::MatrixXcd L(G_b.rows(),G_b.cols());
    // std::cerr << "\nMultiplying G_b by Chi...";
    // // L.setIdentity();
    // L = G_b;
    // for(int cc = 0; cc<Chi.cols();cc++)
    // {
    //     L.col(cc) = L.col(cc)*Chi(cc,cc);
    //     L(cc,cc) += 1.0;
    // }
    // std::cerr << " done!" << std::endl;
    //
    //
    // std::cerr << "\nMaking Ez_sct_all...";
    // Eigen::MatrixXcd Ez_sct_all(
    //     Ez_inc_all.rows(),
    //     Ez_inc_all.cols()
    // );
    // std::cerr << " done!" << std::endl;
    //
    //
    // std::cerr << "\nFactoring L(" << areas.size() << " by " << areas.size() << ")...";
    // double tick_start = std::clock();
    // Eigen::PartialPivLU<Eigen::MatrixXcd> LU_L;
    // LU_L.compute(L);
    // double tick_stop = std::clock();
    // std::cerr << " done in " << (tick_stop-tick_start)/CLOCKS_PER_SEC << " sec" << std::endl;
    //
    // std::cerr << "\nSolving scattered fields for each source...";
    // for(int jj = 0; jj < Ez_inc_all.cols(); jj++)
    // {
    //     Eigen::VectorXcd b = -G_b*(Chi*Ez_inc_all.col(jj));
    //     Eigen::VectorXcd x = LU_L.solve(b);
    //     Ez_sct_all.col(jj) = x;
    //     std::cerr << jj << " ";
    // }
    // std::cerr << "done!" << std::endl;
    //
    // Eigen::MatrixXd rxlocations;
    // Eigen::VectorXcd rxcoefficients;
    //
    // ReadAntennaFile(
    //     rxlocations,
    //     rxcoefficients,
    //     argv[2]
    // );
    //
    // Eigen::MatrixXcd Gd_b;
    // BuildDataGreen(
    //     Gd_b,
    //     centroids,
    //     areas,
    //     rxlocations,
    //     k2_b
    // );
    //
    //
    // std::cerr << "\nWriting results to file...";
    // WriteMatrixToFile(
    //     "Test_tri.txt",
    //     tri
    // );
    //
    // WriteMatrixToFile(
    //     "Test_pts.txt",
    //     points
    // );
    //
    // WriteVectorToFile(
    //     "k2_fgd.txt",
    //     k2_f
    // );
    //
    // WriteMatrixToFile(
    //     "Ez_inc.txt",
    //     Ez_inc_all
    // );
    // WriteMatrixToFile(
    //     "Ez_sct.txt",
    //     Ez_sct_all
    // );
    // std::cerr <<" done!" << std::endl;
    //
    // std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    // std::cerr << " done!" << std::endl;
    return(0);
}
