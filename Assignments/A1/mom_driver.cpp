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


// Antenna::Antenna(void)
// {
//     location << 0,0,0;
//     direction << 0,0,0;
//     coefficient = 1;
//     frequency = 1;
// }

void Antenna::getField(
    Eigen::VectorXcd & Ez,
    const Eigen::MatrixXd & points
)
{
    int n_points = points.rows();
    Eigen::VectorXd distances(n_points);
    for(int dd = 0; dd < distances.size(); dd++)
    {
        distances(dd) = std::sqrt((location - points.row(dd).transpose()).array().pow(2).sum());
    }
    double k = 2*M_PI*frequency/CNAUGHT;
    Ez.resize(n_points);
    for(int ipt = 0; ipt < Ez.size(); ipt++)
    {
        Ez(ipt) = boost::math::cyl_hankel_2(
            0,
            k*distances(ipt)
        );
    }
}

void Mesh::buildTriangulation(
    std::string filename,
    bool verbose
)
{
    gmsh::open(filename);
    tri.resize(0,4);
    std::vector<Eigen::VectorXd> v_tri;
    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords;
    std::vector<double> node_para_coord;
    gmsh::model::mesh::getNodes(
        node_tags,
        node_coords,
        node_para_coord
    );

    //std::cerr << "n_nodes = " << node_tags.size() << std::endl;
    int max_tag = *(std::max_element(node_tags.begin(),node_tags.end()));
    //std::cerr << "Max tag: " << max_tag << std::endl;
    points.resize(max_tag+1,3);
    points.setZero();
    //std::cerr << points << std::endl;
    for(int nn = 0; nn < node_tags.size(); nn++)
    {
        int n_tag = node_tags[nn];
        //std::cerr << "Node with tag = " << n_tag << std::endl;
        Eigen::Vector3d point;
        point << node_coords[3*nn+0],node_coords[3*nn+1],node_coords[3*nn+2];
        //std::cerr << point << std::endl;
        points.row(n_tag) = point.transpose();
        //std::cerr << "assigned a row in points" << std::endl;

    }


    // Now try going physical group based
    gmsh::vectorpair physdimtags;
    gmsh::model::getPhysicalGroups(physdimtags);
    //std::cerr << "Phys Group Based..." << std::endl;
    for(int ii = 0; ii < physdimtags.size(); ii++)
    {
        int physdim = physdimtags[ii].first;
        int phystag = physdimtags[ii].second;
        std::vector<int> ent_tags;
        gmsh::model::getEntitiesForPhysicalGroup(
            physdim,
            phystag,
            ent_tags
        );
        //std::cerr << "\t" <<physdim << "," << phystag << std::endl;
        for(int jj = 0; jj < ent_tags.size(); jj++)
        {
            //std::cerr << "\t\t" << ent_tags[jj] << std::endl;
            std::vector<int> types;
            std::vector<std::vector<std::size_t> > eletags;
            std::vector<std::vector<std::size_t> > nodetags;
            gmsh::model::mesh::getElements(
                types,
                eletags,
                nodetags,
                physdim,
                ent_tags[jj]
            );
            for(int tt = 0; tt < eletags.size(); tt++)
            {
                std::string type_name;
                int type_dim;
                int type_order;
                int type_num_nodes;
                std::vector<double> type_node_coords;
                gmsh::model::mesh::getElementProperties(
                    types[tt],
                    type_name,
                    type_dim,
                    type_order,
                    type_num_nodes,
                    type_node_coords
                );
                //std::cerr << "\t\t\tT=" << types[tt] << ":";
                for(int ee = 0; ee < eletags[tt].size();ee++)
                {
                    //std::cerr << eletags[tt][ee] << "(";
                    Eigen::VectorXd element(type_num_nodes+1);
                    for(int nn = 0; nn < type_num_nodes; nn++)
                    {
                        element(nn) = nodetags[tt][ee*type_num_nodes+nn];
                        //std::cerr << nodetags[tt][ee*type_num_nodes+nn] << ",";
                    }
                    //std::cerr << "),";
                    element(type_num_nodes) = phystag;
                    v_tri.push_back(element);

                }
                //std::cerr << std::endl;
            }
        }
    }
    int n_2d_eles = 0;
    for(int ii = 0; ii < v_tri.size(); ii++)
    {
        //std::cerr << v_tri[ii] << std::endl;
        if(v_tri[ii].size() == 4)
        {
            n_2d_eles++;
        }
    }
    tri.resize(n_2d_eles,4);
    int tri_ptr = 0;
    for(int ii = 0; ii < v_tri.size(); ii++)
    {
        //std::cerr << v_tri[ii] << std::endl;
        if(v_tri[ii].size() == 4)
        {
            tri.row(tri_ptr++) = v_tri[ii].transpose();
            //std::cerr << "\t" << tri.row(tri_ptr-1);
        }
    }
}

Chamber::Chamber(std::string meshfile)
{
    Ez_inc_ready = false;
    Ez_tot_ready = false;
    Ez_sct_ready = false;
    mesh.buildTriangulation(meshfile);
}

void Chamber::addTarget(std::string targetfile)
{
    std::ifstream reader;
    reader.open(targetfile,std::ifstream::in);
    std::string ins;
    reader >> ins;
    assert(ins.compare("tag") == 0);
    reader >> ins;
    assert(ins.compare("eps_rel_real") == 0);
    reader >> ins;
    assert(ins.compare("eps_rel_imag") == 0);

    int tag;
    double eps_rel_real,eps_rel_imag;
    std::complex<double> eps_rel_complex;
    eps_r.resize(mesh.tri.rows());
    for(int rr = 0; rr < eps_r.size(); rr++)
    {
        eps_r(rr) = 1.0;
    }
    while(reader)
    {
        reader >>tag;
        if(reader.eof()) break;
        reader >> eps_rel_real;
        reader >> eps_rel_imag;
        eps_rel_complex.real(eps_rel_real);
        eps_rel_complex.imag(eps_rel_imag);

        for(int rr = 0; rr < mesh.tri.rows(); rr++)
        {
            if(mesh.tri(rr,3) == tag)
            {
                eps_r(rr) = eps_rel_complex;
            }
        }
    }
    reader.close();
}

void Chamber::setupAntennas(std::string antennafile)
{
    std::ifstream reader;
    reader.open(antennafile,std::ifstream::in);
    std::string ins;
    reader >> ins;
    assert(ins.compare("x") == 0);
    reader >> ins;
    assert(ins.compare("y") == 0);
    reader >> ins;
    assert(ins.compare("z") == 0);
    reader >> ins;
    assert(ins.compare("magnitude") == 0);
    reader >> ins;
    assert(ins.compare("phase") == 0);
    reader >> ins;
    assert(ins.compare("style") == 0);
    double x,y,z,mag,phs;
    AntennaStyle_t sty;
    std::string sty_string;
    antennas.resize(0);
    Antenna iant;
    while(reader)
    {
        reader >> x;
        if(reader.eof()) break;
        reader >> y;
        reader >> z;
        reader >> mag;
        reader >> phs;
        reader >> sty_string;
        iant.coefficient = std::polar(mag,phs);
        iant.frequency = frequency;
        if(sty_string.compare("planewave") == 0)
        {
            iant.style = PlaneWave;
            iant.direction << x,y,z;
        }
        else if(sty_string.compare("linesource") == 0)
        {
            iant.style = LineSource;
            iant.location << x,y,z;
        }
        else if(sty_string.compare("Patch") == 0)
        {
            iant.style = Patch;
            iant.location << x,y,z;
        }
        antennas.push_back(iant);
    }
    reader.close();
}

void Chamber::setFrequency(double freq)
{
    frequency = freq;
    for(int aa = 0; aa < antennas.size(); aa++)
    {
        antennas[aa].frequency = freq;
    }
}
void Chamber::getEzTot(Eigen::MatrixXcd & Ezdest)
{

}


void WriteMatrixToFile(
    std::string filename,
    Eigen::MatrixXcd matrix,
    bool append
)
{
    std::ofstream writer;
    if(append)
    {
        writer.open(filename,std::ofstream::out|std::ofstream::app);
    }
    else
    {
        writer.open(filename,std::ofstream::out);
    }
    for(int ii = 0; ii < matrix.rows(); ii++)
    {
        int n_cols= matrix.cols();
        for(int jj = 0; jj <n_cols; jj++)
        {
            std::string sep = ",";
            if (jj == (n_cols-1)) sep = "";
            std::string imag_prefix = "+";
            if(matrix(ii,jj).imag()<0) imag_prefix = "";

            writer << matrix(ii,jj).real() << imag_prefix << matrix(ii,jj).imag() << "i" << sep;
        }
        writer << std::endl;
    }

    writer.close();
}

void WriteMatrixToFile(
    std::string filename,
    Eigen::MatrixXd matrix,
    bool append
)
{
    std::ofstream writer;
    if(append)
    {
        writer.open(filename,std::ofstream::out|std::ofstream::app);
    }
    else
    {
        writer.open(filename,std::ofstream::out);
    }
    for(int ii = 0; ii < matrix.rows(); ii++)
    {
        int n_cols= matrix.cols();
        for(int jj = 0; jj <n_cols; jj++)
        {
            std::string sep = ",";
            if (jj == (n_cols-1)) sep = "";

            writer << matrix(ii,jj) << sep;
        }
        writer << std::endl;
    }
    writer.close();
}

void WriteVectorToFile(
    std::string filename,
    Eigen::VectorXd vec,
    bool append
)
{
    std::ofstream writer;
    if(append)
    {
        writer.open(filename,std::ofstream::out|std::ofstream::app);
    }
    else
    {
        writer.open(filename,std::ofstream::out);
    }
    for(int ii = 0; ii < vec.size(); ii++)
    {
        writer << vec(ii) << std::endl;
    }
    writer.close();
}

void WriteVectorToFile(
    std::string filename,
    Eigen::VectorXcd vec,
    bool append
)
{
    std::ofstream writer;
    if(append)
    {
        writer.open(filename,std::ofstream::out|std::ofstream::app);
    }
    else
    {
        writer.open(filename,std::ofstream::out);
    }
    for(int ii = 0; ii < vec.size(); ii++)
    {
        std::string imag_prefix = "+";
        if(vec(ii).imag()<0) imag_prefix = "";

        writer << vec(ii).real() << imag_prefix << vec(ii).imag() << "i" << std::endl;
    }
    writer.close();
}

void ReadAntennaFile(
    Eigen::MatrixXd & locations,
    Eigen::VectorXcd &  coefficients,
    const std::string filename
)
{
    std::ifstream reader;
    reader.open(filename,std::ifstream::in);
    std::vector<Eigen::VectorXd> position_vec;
    std::vector<std::complex<double> > coeff_vec;
    std::string dummy_line;
    // Demand that we have X Y Z Magnitude Phase
    reader >> dummy_line;
    assert(dummy_line.compare("x") == 0);
    reader >> dummy_line;
    assert(dummy_line.compare("y") == 0);
    reader >> dummy_line;
    assert(dummy_line.compare("z") == 0);
    reader >> dummy_line;
    assert(dummy_line.compare("magnitude") == 0);
    reader >> dummy_line;
    assert(dummy_line.compare("phase") == 0);
    int itx = 0;
    while(reader)
    {
        Eigen::VectorXd position(3);
        double tx,ty,tz;

        reader >> tx;

        if(reader.eof()) break;

        reader >> ty;
        reader >> tz;

        position << tx,ty,tz;
        position_vec.push_back(position);
        double mag,phase;
        reader >> mag;
        reader >> phase;
        coeff_vec.push_back(std::polar(mag,phase));
        itx++;
    }

    int ntx = itx;
    locations.resize(ntx,3);
    coefficients.resize(ntx);

    for (int irow = 0; irow < ntx; irow++)
    {
        Eigen::MatrixXd loc_row(1,3);
        loc_row << position_vec[irow][0],position_vec[irow][1],position_vec[irow][2];
        locations.block(irow,0,1,3) = loc_row;
        coefficients(irow) = coeff_vec[irow];
    }

    reader.close();
}

// void BuildTriangulation(
//     std::string filename,
//     Eigen::MatrixXd & tri,
//     Eigen::MatrixXd & points,
//     bool verbose
// )
// {
//     // tri is a triangulation with N rows, 4 columns
//     // Column 1: index of point 1 in points
//     // Column 2: index of point 2 in points
//     // Column 3: index of point 3 in points
//     // Column 4: group tag
//     // points is the xyz coordinates of each point referenced by the
//     // triangulation.
//
//     // points will have N rows and 3 columns
//     gmsh::open(filename);
//     tri.resize(0,4);
//     std::vector<Eigen::VectorXd> v_tri;
//     std::vector<std::size_t> node_tags;
//     std::vector<double> node_coords;
//     std::vector<double> node_para_coord;
//     gmsh::model::mesh::getNodes(
//         node_tags,
//         node_coords,
//         node_para_coord
//     );
//
//     //std::cerr << "n_nodes = " << node_tags.size() << std::endl;
//     int max_tag = *(std::max_element(node_tags.begin(),node_tags.end()));
//     //std::cerr << "Max tag: " << max_tag << std::endl;
//     points.resize(max_tag+1,3);
//     points.setZero();
//     //std::cerr << points << std::endl;
//     for(int nn = 0; nn < node_tags.size(); nn++)
//     {
//         int n_tag = node_tags[nn];
//         //std::cerr << "Node with tag = " << n_tag << std::endl;
//         Eigen::Vector3d point;
//         point << node_coords[3*nn+0],node_coords[3*nn+1],node_coords[3*nn+2];
//         //std::cerr << point << std::endl;
//         points.row(n_tag) = point.transpose();
//         //std::cerr << "assigned a row in points" << std::endl;
//
//     }
//
//
//     // Now try going physical group based
//     gmsh::vectorpair physdimtags;
//     gmsh::model::getPhysicalGroups(physdimtags);
//     //std::cerr << "Phys Group Based..." << std::endl;
//     for(int ii = 0; ii < physdimtags.size(); ii++)
//     {
//         int physdim = physdimtags[ii].first;
//         int phystag = physdimtags[ii].second;
//         std::vector<int> ent_tags;
//         gmsh::model::getEntitiesForPhysicalGroup(
//             physdim,
//             phystag,
//             ent_tags
//         );
//         //std::cerr << "\t" <<physdim << "," << phystag << std::endl;
//         for(int jj = 0; jj < ent_tags.size(); jj++)
//         {
//             //std::cerr << "\t\t" << ent_tags[jj] << std::endl;
//             std::vector<int> types;
//             std::vector<std::vector<std::size_t> > eletags;
//             std::vector<std::vector<std::size_t> > nodetags;
//             gmsh::model::mesh::getElements(
//                 types,
//                 eletags,
//                 nodetags,
//                 physdim,
//                 ent_tags[jj]
//             );
//             for(int tt = 0; tt < eletags.size(); tt++)
//             {
//                 std::string type_name;
//                 int type_dim;
//                 int type_order;
//                 int type_num_nodes;
//                 std::vector<double> type_node_coords;
//                 gmsh::model::mesh::getElementProperties(
//                     types[tt],
//                     type_name,
//                     type_dim,
//                     type_order,
//                     type_num_nodes,
//                     type_node_coords
//                 );
//                 //std::cerr << "\t\t\tT=" << types[tt] << ":";
//                 for(int ee = 0; ee < eletags[tt].size();ee++)
//                 {
//                     //std::cerr << eletags[tt][ee] << "(";
//                     Eigen::VectorXd element(type_num_nodes+1);
//                     for(int nn = 0; nn < type_num_nodes; nn++)
//                     {
//                         element(nn) = nodetags[tt][ee*type_num_nodes+nn];
//                         //std::cerr << nodetags[tt][ee*type_num_nodes+nn] << ",";
//                     }
//                     //std::cerr << "),";
//                     element(type_num_nodes) = phystag;
//                     v_tri.push_back(element);
//
//                 }
//                 //std::cerr << std::endl;
//             }
//         }
//     }
//     int n_2d_eles = 0;
//     for(int ii = 0; ii < v_tri.size(); ii++)
//     {
//         //std::cerr << v_tri[ii] << std::endl;
//         if(v_tri[ii].size() == 4)
//         {
//             n_2d_eles++;
//         }
//     }
//     tri.resize(n_2d_eles,4);
//     int tri_ptr = 0;
//     for(int ii = 0; ii < v_tri.size(); ii++)
//     {
//         //std::cerr << v_tri[ii] << std::endl;
//         if(v_tri[ii].size() == 4)
//         {
//             tri.row(tri_ptr++) = v_tri[ii].transpose();
//             //std::cerr << "\t" << tri.row(tri_ptr-1);
//         }
//     }
// }

void CalculateTriAreas(
    Eigen::VectorXd & areas,
    const Eigen::MatrixXd & tri,
    const Eigen::MatrixXd & points,
    bool verbose
)
{
    assert(tri.cols() == 4);
    assert(points.cols() == 3);
    areas.resize(tri.rows());
    for(int tri_idx = 0; tri_idx < tri.rows(); tri_idx++)
    {
        Eigen::MatrixXd itri = tri.block(tri_idx,0,1,3);
        int p_idx = itri(0);
        int q_idx = itri(1);
        int r_idx = itri(2);
        double px = points(p_idx,0);
        double py = points(p_idx,1);
        double qx = points(q_idx,0);
        double qy = points(q_idx,1);
        double rx = points(r_idx,0);
        double ry = points(r_idx,1);
        double v1x = qx-px;
        double v1y = qy-py;
        double v2x = rx-px;
        double v2y = ry-py;
        areas(tri_idx) = std::abs((v1x*v2y-v1y*v2x)/2);
    }
}

void CalculateTriCentroids(
    Eigen::MatrixXd & centroids,
    const Eigen::MatrixXd & tri,
    const Eigen::MatrixXd & points
)
{
    assert(tri.cols() == 4);
    assert(points.cols() == 3);
    centroids.resize(tri.rows(),3);
    centroids.setZero();
    for(int tri_idx = 0; tri_idx < tri.rows(); tri_idx++)
    {
        Eigen::VectorXd itri = tri.block(tri_idx,0,1,3).transpose();
        Eigen::VectorXd tx(3),ty(3),tz(3);
        tx << points(itri(0),0),points(itri(1),0),points(itri(2),0);
        ty << points(itri(0),1),points(itri(1),1),points(itri(2),1);
        tz << points(itri(0),2),points(itri(1),2),points(itri(2),2);
        Eigen::VectorXd centroid(3);
        centroid << tx.array().mean(),ty.array().mean(),tz.array().mean();
        centroids.block(tri_idx,0,1,3) = centroid.transpose();
    }
}





void BuildDomainGreen(
    Eigen::MatrixXcd & G,
    const Eigen::MatrixXd & centroids,
    const Eigen::VectorXd & areas,
    double k2_b
)
{
    int n_tri = areas.size();

    G.resize(n_tri,n_tri);

    std::complex<double> j(0,1);

    double k_b = std::sqrt(k2_b);

    for(int mm = 0; mm < n_tri; mm++)
    {
        for(int nn = mm; nn < n_tri; nn++)
        {
            //std::cerr << mm << " " << nn;
            Eigen::VectorXd dxyz = (centroids.row(nn)-centroids.row(mm)).transpose();
            dxyz = dxyz.array().pow(2);
            double dmn = std::sqrt(dxyz.array().sum());

            std::complex<double> Gmn;
            double a_n = std::sqrt(areas[nn]/M_PI);
            //std::cerr << "dmn = " << dmn << ", a_n = " << a_n << ",";
            if(mm==nn)
            {
                //std::cerr << " " << k_b*a_n;
                std::complex<double> H12 = boost::math::cyl_hankel_2(
                    1,
                    k_b*a_n
                );
                Gmn = -j/(2.0*k_b*k_b)*(M_PI*k_b*a_n*H12 - 2.0*j);
            }
            else
            {
                //std::cerr << " " << k_b*a_n;
                std::complex<double> J1 = boost::math::cyl_bessel_j(
                    1,
                    k_b*a_n
                );
                //std::cerr << " " << k_b*dmn;
                std::complex<double> H02 = boost::math::cyl_hankel_2(
                    0,
                    k_b*dmn
                );
                // Gmn = j*M_PI*k_b*a_n*J1*H02/2.0;
                Gmn = -j*M_PI*a_n*J1*H02/2.0/k_b;
            }
            //std::cerr << std::endl;
            // Gmn *= -j/4.0;
            G(mm,nn) = Gmn;
            G(nn,mm) = Gmn;
        }
    }
}

void BuildDataGreen(
    Eigen::MatrixXcd & G,
    const Eigen::MatrixXd & centroids,
    const Eigen::VectorXd & areas,
    const Eigen::MatrixXd & rxlocations,
    double k2_b
)
{
    int n_rx = rxlocations.rows();
    int n_ele = areas.size();
    G.resize(n_rx,n_ele);
    std::complex<double> j(0,1);
    double k_b = std::sqrt(k2_b);
    for(int rr = 0; rr < n_rx; rr++)
    {

        for(int ee = 0; ee < n_ele; ee++)
        {
            Eigen::VectorXd dxyz = (rxlocations.row(rr)-centroids.row(ee)).transpose();
            dxyz = dxyz.array().pow(2);
            double d_re = std::sqrt(dxyz.array().sum());
            double a_e = std::sqrt(areas[ee]/M_PI);

            std::complex<double> J1 = boost::math::cyl_bessel_j(
                1,
                k_b*a_e
            );
            //std::cerr << " " << k_b*dmn;
            std::complex<double> H02 = boost::math::cyl_hankel_2(
                0,
                k_b*d_re
            );
            // Gmn = j*M_PI*k_b*a_n*J1*H02/2.0;
            std::complex<double> G_re = -j*M_PI*a_e*J1*H02/2.0/k_b;
            G(rr,ee) = G_re;
        }
    }
}
