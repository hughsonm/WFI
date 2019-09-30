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
    double k = 2*M_PI*frequency/CNAUGHT;
    Ez.resize(n_points);
    std::complex<double> j_imag(0,1);

    switch(style)
    {
        case LineSource:
        {
            Eigen::VectorXd distances(n_points);
            for(int dd = 0; dd < distances.size(); dd++)
            {
                distances(dd) = std::sqrt((location - points.row(dd).transpose()).array().pow(2).sum());
            }
            for(int ipt = 0; ipt < Ez.size(); ipt++)
            {
                Ez(ipt) = boost::math::cyl_hankel_2(
                    0,
                    k*distances(ipt)
                );
            }
            break;
        }
        case PlaneWave:
        {
            for(int ipt = 0; ipt < Ez.size(); ipt++)
            {
                Eigen::Vector3cd k_vector = direction*k;
                std::complex<double> phase = points.row(ipt)*k_vector;
                Ez(ipt) = std::exp(j_imag*phase);
            }
            break;
        }
        case Patch:
        {
            assert(0==1);
        }
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
    std::cerr << "Built triangulation. Found ";
    std::cerr << tri.rows();
    std::cerr << " triangles, and ";
    std::cerr << points.rows();
    std::cerr << " points." << std::endl;

    std::cerr << "Calculating tri areas...";
    CalculateTriAreas();
    std::cerr << "done!" << std::endl;
    std::cerr << "Calculating tri centroids...";
    CalculateTriCentroids();
    std::cerr << "done!" << std::endl;
}

void Mesh::CalculateTriAreas()
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

void Mesh::CalculateTriCentroids()
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

void Mesh::buildDataGreen(
    Eigen::MatrixXcd & G,
    std::complex<double> k2_b,
    const Eigen::MatrixXd & locations
)
{
    // For now, assume we go from all transmitters to all receivers
    std::cerr << "Resizing data green to " << locations.rows() << " by " << areas.size() << std::endl;
    G.resize(locations.rows(),areas.size());
    std::cerr << "success!" << std::endl;

    double k_b = std::sqrt(k2_b).real();
    std::complex<double> j_imag(0,1);

    for(int ll = 0; ll < locations.rows(); ll++)
    {
        std::cerr << "Calculating data green row for receiver at position:" << std::endl;
        std::cerr << locations.row(ll) << std::endl;
        for(int aa = 0; aa < areas.size(); aa++)
        {
            std::cerr << "element positon: " << std::endl;
            std::cerr << centroids.row(aa) << std::endl;
            Eigen::VectorXd diff = locations.row(ll) - centroids.row(aa);
            double distance = std::sqrt((diff.transpose())*diff);
            double radius = std::sqrt(areas[aa]/M_PI);
            std::cerr << "distance = " << distance << ", radius = " << radius << std::endl;
            std::complex<double> J1 = boost::math::cyl_bessel_j(
                1,
                k_b*radius
            );
            //std::cerr << " " << k_b*dmn;
            std::complex<double> H02 = boost::math::cyl_hankel_2(
                0,
                k_b*distance
            );
            // Gmn = j*M_PI*k_b*a_n*J1*H02/2.0;
            G(ll,aa) = -j_imag*M_PI*radius*J1*H02/2.0/k_b;
        }
    }
}

void Mesh::buildDomainGreen(
    Eigen::MatrixXcd & G,
    std::complex<double> k2_b
)
{
    int n_tri = areas.size();
    G.resize(n_tri,n_tri);

    std::complex<double> j(0,1);

    double k_b = (std::sqrt(k2_b)).real();

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


Chamber::Chamber(std::string meshfile)
{
    Ez_inc_ready = false;
    Ez_tot_ready = false;
    Ez_sct_ready = false;
    G_b_domain_ready = false;
    G_b_data_ready = false;
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
    k2_f.resize(mesh.tri.rows());
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
                k2_f(rr) = k2_b*eps_r(rr);
            }
        }
    }
    reader.close();
}

void Chamber::setupProbes(std::string probefile)
{
    std::ifstream reader;
    reader.open(probefile,std::ifstream::in);
    std::string ins;
    reader >> ins;
    assert(ins.compare("x") == 0);
    reader >> ins;
    assert(ins.compare("y") == 0);
    reader >> ins;
    assert(ins.compare("z") == 0);
    double x,y,z;
    std::vector<Eigen::Vector3d> vec_of_probes;
    while(reader)
    {
        Eigen::Vector3d probe_xyz;
        reader >> x;
        if(reader.eof()) break;
        reader >> y;
        reader >> z;
        probe_xyz << x,y,z;
        vec_of_probes.push_back(probe_xyz);
    }
    probe_points.resize(vec_of_probes.size(),3);
    for(int pp = 0; pp < vec_of_probes.size(); pp++)
    {
        probe_points.row(pp) = vec_of_probes[pp].transpose();
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
            iant.direction /= std::sqrt((iant.direction.transpose()*iant.direction));
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
    Ez_inc_ready = false;
    Ez_sct_ready = false;
    Ez_tot_ready = false;
    G_b_domain_ready = false;
}

void Chamber::setFrequency(double freq)
{
    frequency = freq;
    for(int aa = 0; aa < antennas.size(); aa++)
    {
        antennas[aa].frequency = freq;
    }
    double omega = 2*M_PI*freq;
    k2_b = omega*omega/CNAUGHT/CNAUGHT;
    k2_f.resize(eps_r.size());
    for(int ee = 0; ee<eps_r.size(); ee++)
    {
        k2_f(ee) = k2_b*eps_r(ee);
    }
    Ez_inc_ready = false;
    Ez_sct_ready = false;
    Ez_tot_ready = false;
    G_b_domain_ready = false;
}
void Chamber::calcDomainEzTot(void)
{
    int entry_point = 0;
    if(Ez_tot_ready)
    {
        entry_point = 3;
    }
    else if(Ez_inc_ready && Ez_sct_ready)
    {
        entry_point = 2;
    }
    else if(G_b_domain_ready)
    {
        entry_point = 1;
    }
    else
    {
        entry_point = 0;
    }
    switch(entry_point)
    {
        case 0:
        {
            // Calculate incident fields.
            Ez_inc.resize(mesh.areas.size(),antennas.size());
            for(int aa = 0; aa < antennas.size(); aa++)
            {
                Eigen::VectorXcd Ez_inc_a;
                antennas[aa].getField(
                    Ez_inc_a,
                    mesh.centroids
                );
                Ez_inc.col(aa) = Ez_inc_a;
            }
            Ez_inc_ready = true;

            // Build domain green matrix
            std::cerr << "Building domain green..." << std::endl;
            mesh.buildDomainGreen(
                G_b_domain,
                k2_b
            );
            std::cerr << "G:(" << G_b_domain.rows() << "," << G_b_domain.cols() << ")" << std::endl;
            std::cerr << "Building Chi..." << std::endl;

            // Build contrast matrix
            Chi.resize(k2_f.size(),k2_f.size());
            Chi.setZero();
            for(int dd = 0; dd < Chi.cols();dd++)
            {
                Chi(dd,dd) = k2_f(dd)-k2_b;
            }

            // Build domain L operator
            std::cerr << "Building L_domain" << std::endl;
            L_domain.resize(G_b_domain.rows(),G_b_domain.cols());
            L_domain = G_b_domain;
            for(int cc = 0; cc< Chi.cols(); cc++)
            {
                L_domain.col(cc) *= Chi(cc,cc);
                L_domain(cc,cc) += 1.0;
            }
            std::cerr << "filled,";

            // Perform LU factorization of domain L operator
            LU_L.compute(L_domain);
            std::cerr << "factored" << std::endl;
        }
        case 1:
        {
            // Allocate space for incident and total fields.
            Ez_tot.resize(mesh.areas.size(),antennas.size());
            Ez_tot.setZero();
            Ez_sct.resize(mesh.areas.size(),antennas.size());
            Ez_sct.setZero();
            // Calculate all the scattered fields
            for(int jj = 0; jj < Ez_tot.cols(); jj++)
            {
                Ez_sct.col(jj) = LU_L.solve(
                    -G_b_domain*(
                        Chi*Ez_inc.col(jj)
                    )
                );
            }
            Ez_sct_ready = true;
        }
        case 2:
        {
            // Add scattered to indident to get total fields
            Ez_tot = Ez_inc + Ez_sct;
            Ez_tot_ready = true;
        }
        case 3:
        {
            // Write out the total fields
            // Ezdest.resize(Ez_tot.rows(),Ez_tot.cols());
            // Ezdest = Ez_tot;
        }
    }
}

void Chamber::calcDataEzTot(void)
{
    std::cerr << "Probe points are..." << std::endl;
    std::cerr << probe_points << std::endl;
    if(!Ez_tot_ready)
    {
        std::cerr << "Ez_tot was not built. Getting it now." << std::endl;
        calcDomainEzTot();
    }
    if(!G_b_data_ready)
    {
        std::cerr << "Data green was not built. Getting it now" << std::endl;
        mesh.buildDataGreen(
            G_b_data,
            k2_b,
            probe_points
        );
        G_b_data_ready = true;
    }
    std::cerr << "Resizing d-field holders to " << probe_points.rows() << " by " << antennas.size() << std::endl;
    Ez_inc_d.resize(probe_points.rows(),antennas.size());
    Ez_sct_d.resize(probe_points.rows(),antennas.size());
    Ez_tot_d.resize(probe_points.rows(),antennas.size());

    for(int aa = 0; aa < antennas.size(); aa++)
    {
        std::cerr << "Calculating d-inc for antenna " << aa;
        Eigen::VectorXcd Ez_inc_a;
        antennas[aa].getField(
            Ez_inc_a,
            probe_points
        );
        Ez_inc_d.col(aa) = Ez_inc_a;
        std::cerr << " done " << std::endl;
    }
    for(int tt = 0; tt < Ez_tot.cols(); tt++)
    {
        std::cerr << "Calculating d-sct for antenna " << tt;
        Ez_sct_d.col(tt) = G_b_data*(Chi*(Ez_tot.col(tt)));
        std::cerr << "done." << std::endl;

    }
    Ez_tot_d = Ez_inc_d + Ez_sct_d;

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
