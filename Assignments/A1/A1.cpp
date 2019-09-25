#include <eigen3/Eigen/Eigen>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>

#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675

void WriteMatrixToFile(
    std::string filename,
    Eigen::MatrixXcd matrix,
    bool append=false
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
    bool append=false
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
    bool append=false
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
    bool append=false
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

void BuildTriangulation(
    std::string filename,
    Eigen::MatrixXd & tri,
    Eigen::MatrixXd & points,
    bool verbose = false
)
{
    // tri is a triangulation with N rows, 4 columns
    // Column 1: index of point 1 in points
    // Column 2: index of point 2 in points
    // Column 3: index of point 3 in points
    // Column 4: group tag
    // points is the xyz coordinates of each point referenced by the
    // triangulation.

    // points will have N rows and 3 columns
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

void CalculateTriAreas(
    Eigen::VectorXd & areas,
    const Eigen::MatrixXd & tri,
    const Eigen::MatrixXd & points,
    bool verbose = false
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

void EvaluateIncidentField(
    Eigen::VectorXcd & obs_Ez,
    const Eigen::VectorXd source_loc,
    const std::complex<double> source_coeff,
    const double frequency,
    const Eigen::MatrixXd & obs_pts
)
{
    assert(source_loc.size() == obs_pts.cols());

    Eigen::MatrixXd d2 = obs_pts;
    for(int obs_idx = 0; obs_idx < d2.rows(); obs_idx++)
    {
        d2.block(obs_idx,0,1,3) -= source_loc.transpose();
    }

    d2 = d2.array().pow(2);
    Eigen::VectorXd r2 = d2.array().rowwise().sum();

    double k = 2*M_PI*frequency/CNAUGHT;

    Eigen::VectorXd hankelarg = r2.array().sqrt();

    hankelarg *= k;
    obs_Ez.resize(hankelarg.size());
    obs_Ez.setZero();
    for (int ihankel = 0; ihankel < hankelarg.size(); ihankel++)
    {
        obs_Ez[ihankel] = boost::math::cyl_hankel_2(0,hankelarg[ihankel]);
    }

}

void AssignRelativeConstitutives(
    Eigen::VectorXcd & eps_r,
    const Eigen::MatrixXd & tri,
    const std::string constfilename
)
{
    std::ifstream reader;
    reader.open(constfilename,std::ifstream::in);
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
    eps_r.resize(tri.rows());
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

        for(int rr = 0; rr < tri.rows(); rr++)
        {
            if(tri(rr,3) == tag)
            {
                eps_r(rr) = eps_rel_complex;
            }
        }
    }
    reader.close();
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

int main (int argc, char **argv)
{
    assert(argc==4);
    std::cerr << "\nInitializing gmsh...";
    gmsh::initialize();
    std::cerr << " done!" << std::endl;

    double frequency = std::atof(argv[3]);
    Eigen::MatrixXd tri;
    Eigen::MatrixXd points;
    std::cerr << "Building triangulation...";
    BuildTriangulation(
        argv[1],
        tri,
        points,
        true
    );
    std::cerr << " done!" << std::endl;

    Eigen::VectorXcd eps_r;

    AssignRelativeConstitutives(
        eps_r,
        tri,
        "Tag2EpsRel.txt"
    );

    std::cerr << "\nCalculating triangle areas...";
    Eigen::VectorXd areas;
    CalculateTriAreas(
        areas,
        tri,
        points
    );

    std::cerr << " done!" << std::endl;

    std::cerr << "\nCalculating triangle centroids...";
    Eigen::MatrixXd centroids;
    CalculateTriCentroids(
        centroids,
        tri,
        points
    );

    std::cerr << " done!" << std::endl;

    std::cerr << "\nReading antenna file...";
    Eigen::MatrixXd tx_locations;
    Eigen::VectorXcd tx_coefficients;
    ReadAntennaFile(
        tx_locations,
        tx_coefficients,
        argv[2]
    );
    std::cerr << " done!" << std::endl;

    std::cerr << "\nCalculating incident fields...";
    Eigen::MatrixXcd Ez_inc_all(centroids.rows(),tx_locations.rows());
    Eigen::VectorXcd Ez_inc_one;
    for(int itx = 0; itx < tx_locations.rows(); itx++)
    {
        EvaluateIncidentField(
            Ez_inc_one,
            tx_locations.row(itx),
            tx_coefficients(itx),
            frequency,
            centroids
        );
        std::cerr << itx << " ";
        Ez_inc_all.col(itx) = Ez_inc_one;
    }
    std::cerr << " done!" << std::endl;

    double omega = 2*M_PI*frequency;
    double k2_b = omega*omega/CNAUGHT/CNAUGHT;

    Eigen::VectorXcd k2_f(eps_r.size());
    k2_f = eps_r.array()*k2_b;

    std::cerr << "\nCalculating background Green function...";
    Eigen::MatrixXcd G_b;
    BuildDomainGreen(
        G_b,
        centroids,
        areas,
        k2_b
    );
    std::cerr << " done!" << std::endl;
    std::cerr << "\nFilling Chi's diagonal with k2...";
    Eigen::MatrixXcd Chi(areas.size(),areas.size());
    Chi.setZero();
    for(int dd = 0; dd < Chi.rows(); dd++)
    {
        Chi(dd,dd) = k2_f(dd) - k2_b;
    }
    std::cerr << " done!" << std::endl;

    std::cerr << "\nAllocating space for L...";
    Eigen::MatrixXcd L(G_b.rows(),G_b.cols());
    std::cerr << "\nMultiplying G_b by Chi...";
    // L.setIdentity();
    L = G_b;
    for(int cc = 0; cc<Chi.cols();cc++)
    {
        L.col(cc) = L.col(cc)*Chi(cc,cc);
        L(cc,cc) += 1.0;
    }
    std::cerr << " done!" << std::endl;


    std::cerr << "\nMaking Ez_sct_all...";
    Eigen::MatrixXcd Ez_sct_all(
        Ez_inc_all.rows(),
        Ez_inc_all.cols()
    );
    std::cerr << " done!" << std::endl;


    std::cerr << "\nFactoring L(" << areas.size() << " by " << areas.size() << ")...";
    double tick_start = std::clock();
    Eigen::PartialPivLU<Eigen::MatrixXcd> LU_L;
    LU_L.compute(L);
    double tick_stop = std::clock();
    std::cerr << " done in " << (tick_stop-tick_start)/CLOCKS_PER_SEC << " sec" << std::endl;

    std::cerr << "\nSolving scattered fields for each source...";
    for(int jj = 0; jj < Ez_inc_all.cols(); jj++)
    {
        Eigen::VectorXcd b = -G_b*(Chi*Ez_inc_all.col(jj));
        Eigen::VectorXcd x = LU_L.solve(b);
        Ez_sct_all.col(jj) = x;
        std::cerr << jj << " ";
    }
    std::cerr << "done!" << std::endl;

    Eigen::MatrixXd rxlocations;
    Eigen::VectorXcd rxcoefficients;

    ReadAntennaFile(
        rxlocations,
        rxcoefficients,
        argv[2]
    );

    Eigen::MatrixXcd Gd_b;
    BuildDataGreen(
        Gd_b,
        centroids,
        areas,
        rxlocations,
        k2_b
    );


    std::cerr << "\nWriting results to file...";
    WriteMatrixToFile(
        "Test_tri.txt",
        tri
    );

    WriteMatrixToFile(
        "Test_pts.txt",
        points
    );

    WriteVectorToFile(
        "k2_fgd.txt",
        k2_f
    );

    WriteMatrixToFile(
        "Ez_inc.txt",
        Ez_inc_all
    );
    WriteMatrixToFile(
        "Ez_sct.txt",
        Ez_sct_all
    );
    std::cerr <<" done!" << std::endl;

    std::cerr << "\nFinalizing gmsh...";
    gmsh::finalize();
    std::cerr << " done!" << std::endl;
    return(0);
}
