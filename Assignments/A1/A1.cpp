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
    // use MatrixXd member function conservativeResize(m,n).
    tri.resize(0,4);
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);

    // Read in the 2d physical groups
    for(int phys_idx = 0; phys_idx < dimTags.size(); phys_idx++)
    {
        if(verbose)
        {
            std::cerr << "Physical Group Number: " << phys_idx;
            std::cerr << std::endl;
        }
        int phys_dim = dimTags[phys_idx].first;
        int phys_tag = dimTags[phys_idx].second;
        if(phys_dim != 2)
        {
            continue;
        }
        if(verbose)
        {
            std::cerr << "\tDimension: " << phys_dim << std::endl;
            std::cerr << "\tTag: " << phys_tag << std::endl;
        }
        std::string phys_name;
        gmsh::model::getPhysicalName(
            phys_dim,
            phys_tag,
            phys_name
        );

        if(verbose)
        {
            std::cerr << "\tGetting entities for this physical group's";
            std::cerr << " dimension" << std::endl;
        }
        gmsh::vectorpair entities_dimtags_for_phys;
        gmsh::model::getEntities(
            entities_dimtags_for_phys,
            phys_dim
        );
        if(verbose)
        {
            std::cerr << "\tPhysical Group Name: ";
            std::cerr << phys_name << std::endl;
        }
        std::vector<int> ent_tags;
        gmsh::model::getEntitiesForPhysicalGroup(
            phys_dim,
            phys_tag,
            ent_tags
        );
        if(verbose)
        {
            std::cerr << "\tEntity Count: " << ent_tags.size();
            std::cerr << std::endl;
        }

        std::vector<std::vector<int> > node_tags_for_ele_for_phys;

        for(int ent_idx = 0; ent_idx < ent_tags.size(); ent_idx++)
        {
            int ent_for_phys_dim = entities_dimtags_for_phys[ent_idx].first;
            int ent_for_phys_tag = entities_dimtags_for_phys[ent_idx].second;

            if(verbose)
            {
                std::cerr << "\t\tThis physical group has the entity";
                std::cerr << " with dim = ";
                std::cerr << ent_for_phys_dim;
                std::cerr << " and tag = ";
                std::cerr << ent_for_phys_tag;
                std::cerr << std::endl;
            }

            std::string ent_for_phys_type;
            gmsh::model::getType(
                ent_for_phys_dim,
                ent_for_phys_tag,
                ent_for_phys_type
            );

            if(verbose)
            {
                std::cerr << "\t\tGetting elements and nodes for this";
                std::cerr << " entity" << std::endl;
            }

            std::vector<int> ele_for_ent_types;
            std::vector<std::vector<std::size_t> > ele_for_ent_tags;
            std::vector<std::vector<std::size_t> > node_for_ent_tags;

            gmsh::model::mesh::getElements(
                ele_for_ent_types,
                ele_for_ent_tags,
                node_for_ent_tags,
                ent_for_phys_dim,
                ent_for_phys_tag
            );
            int node_tag_ptr = 0;
            if(verbose)
            {
                std::cerr << "\t\tGot " << ele_for_ent_types.size();
                std::cerr << " different types of elements" << std::endl;
            }
            for (int type_idx = 0; type_idx < ele_for_ent_types.size(); type_idx++)
            {
                if(verbose)
                {
                    std::cerr << "\t\tType no. " << type_idx;
                    std::cerr << " is " << ele_for_ent_types[type_idx] << std::endl;
                }

                std::string type_name;
                int type_dim;
                int type_order;
                int type_num_nodes;
                std::vector<double> type_node_coords;
                gmsh::model::mesh::getElementProperties(
                    ele_for_ent_types[type_idx],
                    type_name,
                    type_dim,
                    type_order,
                    type_num_nodes,
                    type_node_coords
                );

                if(verbose)
                {
                    std::cerr << "\t\tType Name: " << type_name << std::endl;
                    std::cerr << "\t\tDimension: " << type_dim << std::endl;
                    std::cerr << "\t\tOrder    : " << type_order << std::endl;
                    std::cerr << "\t\tNum Nodes: " << type_num_nodes << std::endl;
                    std::cerr << "\t\tCount    : " << ele_for_ent_tags[type_idx].size() << std::endl;
                    std::cerr << "\t\tTags     : " << std::endl;
                }
                // How many elements we got? A lot.
                int n_new_elements = ele_for_ent_tags[type_idx].size();
                Eigen::Index cur_n_rows = tri.rows();
                if(verbose)
                {
                    std::cerr << "\t\tTri has " << cur_n_rows << " rows";
                    std::cerr << " and we add " << n_new_elements;
                    std::cerr << " for a total of ";
                    std::cerr << (cur_n_rows + n_new_elements);
                    std::cerr << std::endl;
                }


                tri.conservativeResize(
                    cur_n_rows + n_new_elements,
                    Eigen::NoChange
                );
                for(int ele_idx = 0; ele_idx < n_new_elements; ele_idx++)
                {
                    //std::cerr << "\t\t" << ele_for_ent_tags[type_idx][ele_idx] << std::endl;
                    //ele_tags_for_phys.push_back(ele_for_ent_tags[type_idx][ele_idx]);
                    int add_row = cur_n_rows + ele_idx;

                    for(int cc = 0; cc < type_num_nodes; cc++)
                    {
                        tri(add_row,cc) = node_for_ent_tags[type_idx][node_tag_ptr++];
                    }
                    tri(add_row,3) = phys_tag;
                }
            }
        }

        if(verbose)
        {
            std::cerr << "\tThis physical's node tags for each ele:";
            std::cerr << std::endl;
            for(int ele_idx = 0; ele_idx < node_tags_for_ele_for_phys.size(); ele_idx++)
            {
                std::cerr << "\t";
                for(int node_idx = 0; node_idx < node_tags_for_ele_for_phys[ele_idx].size(); node_idx++)
                {
                    std::cerr << node_tags_for_ele_for_phys[ele_idx][node_idx] << "\t";
                }
                std::cerr << std::endl;
            }
        }
    }

    // Read in all of the points (gmsh calls them nodes)
    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords;
    std::vector<double> node_para_coords;
    gmsh::model::mesh::getNodes(
        node_tags,
        node_coords,
        node_para_coords
    );

    int max_tag = *(max_element(node_tags.begin(),node_tags.end()));
    points.resize(max_tag+1,3);
    points.setZero();
    int node_coord_ptr = 0;
    for(int node_idx = 0; node_idx < node_tags.size(); node_idx++)
    {
        int node_row = node_tags[node_idx];
        for(int cc = 0; cc < 3; cc++)
        {
            points(node_row,cc) = node_coords[node_coord_ptr++];
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
        areas(tri_idx) = (v1x*v2y-v1y*v2x)/2;
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

void BuildBackgroundGreen(
    Eigen::MatrixXcd & G,
    Eigen::MatrixXd & centroids,
    Eigen::VectorXd & areas,
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
            Eigen::VectorXd dxyz = (centroids.row(nn)-centroids.row(mm)).transpose();
            dxyz = dxyz.array().pow(2);
            double dmn = std::sqrt(dxyz.array().sum());
            std::complex<double> Gmn;
            double a_n = std::sqrt(areas[nn]/M_PI);
            if(mm==nn)
            {
                std::complex<double> H12 = boost::math::cyl_hankel_2(
                    1,
                    k_b*a_n
                );
                Gmn = j*(M_PI*k_b*a_n*H12 - 2.0*j)/2.0;
            }
            else
            {
                std::complex<double> J1 = boost::math::cyl_bessel_j(
                    1,
                    k_b*a_n
                );
                std::complex<double> H02 = boost::math::cyl_hankel_2(
                    0,
                    k_b*dmn
                );
                Gmn = j*M_PI*k_b*a_n*J1*H02/2.0;
            }
            Gmn *= -j/4.0;
            G(mm,nn) = Gmn;
            G(nn,mm) = Gmn;
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
        false
    );
    std::cerr << " done!" << std::endl;

    std::cerr << "\nCalculating triangle areas...";
    Eigen::VectorXd areas;
    CalculateTriAreas(
        areas,
        tri,
        points
    );
    std::cerr << " done!" << std::endl;

    std::cerr << "\nCalculating triangle areas...";
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

    Eigen::VectorXcd k2_f(areas.size());
    k2_f.setZero();
    k2_f = k2_f.array() + (1.00*k2_b);
    for(int ik2 = 0; ik2 < k2_f.size(); ik2++)
    {
        if(
            (-0.05 < centroids(ik2,0)) &&
            (centroids(ik2,0) < 0.05)  &&
            (-0.05 < centroids(ik2,1)) &&
            (centroids(ik2,1) < 0.05)
        )
        {
            k2_f(ik2) = k2_b*1.1;
        }
    }

    std::cerr << "\nCalculating background Green function...";
    Eigen::MatrixXcd G_b;
    BuildBackgroundGreen(
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
