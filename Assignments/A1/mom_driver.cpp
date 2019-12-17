#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675

#include "mom_driver.h"

#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <boost/math/special_functions/hankel.hpp>
#include <delaunator.hpp>

static void delaunay(
    Eigen::MatrixXi & tri,
    const Eigen::MatrixXd & pts
){
    assert(pts.cols()==2);
    std::vector<double> coordinates;
    for(auto rr{0}; rr<pts.rows(); ++rr){
        coordinates.push_back(pts(rr,0));
        coordinates.push_back(pts(rr,1));
    }
    delaunator::Delaunator del(coordinates);
    const auto ntri {del.triangles.size()/3};
    assert(ntri*3 == del.triangles.size());
    tri.resize(ntri,3);
    for(auto tt{0};tt<ntri;++tt){
        for(auto pp{0}; pp<3; ++pp){
            tri(tt,pp) = del.triangles[3*tt+pp];
        }
    }
}

static void tift(
    Eigen::VectorXcd& xx,
    const Eigen::MatrixXi& tri,
    const Eigen::MatrixXd& xhatpts,
    const Eigen::VectorXcd& xhat,
    const Eigen::MatrixXd& xpts
)
{
    xx.resize(xpts.rows());
    std::complex<double> j_imag(0,1);
    xx.setZero();
    for(auto tt{0}; tt<tri.rows();++tt){
        const Eigen::VectorXd pp{xhatpts.row(tri(tt,0))};
        const Eigen::VectorXd vv{xhatpts.row(tri(tt,1))-pp.transpose()};
        const Eigen::VectorXd uu{xhatpts.row(tri(tt,2))-pp.transpose()};
        const auto J{std::abs(vv(0)*uu(1)-vv(1)*uu(0))};
        std::complex<double> chi_hat_on_tri(0.0,0.0);
        for(auto index{0};index < 3; ++index){
            auto chi_hat_index{tri(tt,index)};
            chi_hat_on_tri += xhat(chi_hat_index);
        }
        chi_hat_on_tri /= 3.0;
        for(auto ii{0}; ii<xpts.rows(); ++ii){
            const Eigen::VectorXd rr{xpts.row(ii)};
            const auto pdr{pp.dot(rr)};
            const auto udr{uu.dot(rr)};
            const auto vdr{vv.dot(rr)};
            const auto umvdr{udr-vdr};
            auto gg{
                (std::exp(j_imag*udr)-std::exp(j_imag*vdr))/
                (j_imag*umvdr)
            };
            gg -= (std::exp(j_imag*udr)-1.0)/(j_imag*udr);
            gg *= J*chi_hat_on_tri*std::exp(j_imag*pdr);
            gg /= (j_imag*vdr);
            xx(ii) += gg;
        }
    }
    xx /= (2*M_PI)*(2*M_PI);
}

void Antenna::getEz(
    Eigen::VectorXcd & Ez,
    const Eigen::MatrixXd & points,
    double frequency
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
                if(1e-10 < std::abs(distances(ipt))){
                    Ez(ipt) = -j_imag/4.0*boost::math::cyl_hankel_2(
                        0,
                        k*distances(ipt)
                    );
                } else{
                    Ez(ipt) = 0;
                }
            }
            break;
        }
        case PlaneWave:
        {
            for(int ipt = 0; ipt < Ez.size(); ipt++)
            {
                Eigen::Vector3cd k_vector = direction*k;
                std::complex<double> phase = points.row(ipt)*k_vector;
                Ez(ipt) = std::exp(-j_imag*phase);
            }
            break;
        }
        case Patch:
        {
            std::cout << "Not yet implemented\n";
            assert(false);
        }
    }
}

void Mesh::buildTriangulation(
    std::string filename,
    bool verbose
)
{
    meshfilename = filename;
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

    int max_tag = *(std::max_element(
        node_tags.begin(),node_tags.end()
    ));
    points.resize(max_tag+1,3);
    points.setZero();
    for(int nn = 0; nn < node_tags.size(); nn++)
    {
        int n_tag = node_tags[nn];
        Eigen::Vector3d point;
        point << node_coords[3*nn+0],node_coords[3*nn+1],node_coords[3*nn+2];
        points.row(n_tag) = point.transpose();
    }

    gmsh::vectorpair physdimtags;
    gmsh::model::getPhysicalGroups(physdimtags);
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
        for(int jj = 0; jj < ent_tags.size(); jj++)
        {
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
                for(int ee = 0; ee < eletags[tt].size();ee++)
                {
                    Eigen::VectorXd element(type_num_nodes+1);
                    for(int nn = 0; nn < type_num_nodes; nn++)
                    {
                        element(nn) = nodetags[tt][ee*type_num_nodes+nn];
                    }
                    element(type_num_nodes) = phystag;
                    v_tri.push_back(element);
                }
            }
        }
    }
    int n_2d_eles = 0;
    for(int ii = 0; ii < v_tri.size(); ii++)
    {
        if(v_tri[ii].size() == 4)
        {
            n_2d_eles++;
        }
    }
    tri.resize(n_2d_eles,4);
    int tri_ptr = 0;
    for(int ii = 0; ii < v_tri.size(); ii++)
    {
        if(v_tri[ii].size() == 4)
        {
            tri.row(tri_ptr++) = v_tri[ii].transpose();
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
    ready = true;
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
    G.resize(locations.rows(),areas.size());

    double k_b = std::sqrt(k2_b).real();
    std::complex<double> j_imag(0,1);

    for(int ll = 0; ll < locations.rows(); ll++)
    {
        for(int aa = 0; aa < areas.size(); aa++)
        {
            Eigen::VectorXd diff = locations.row(ll) - centroids.row(aa);
            double distance = std::sqrt((diff.transpose())*diff);
            double radius = std::sqrt(areas(aa)/M_PI);
            std::complex<double> J1 = boost::math::cyl_bessel_j(
                1,
                k_b*radius
            );
            std::complex<double> H02 = boost::math::cyl_hankel_2(
                0,
                k_b*distance
            );
            G(ll,aa) = -j_imag*M_PI*radius*J1*H02/2.0/k_b;
        }
    }
}

void Mesh::buildDomainGreen(
    Eigen::MatrixXcd & G,
    std::complex<double> k2_b
)
{
    auto n_tri{areas.size()};
    G.resize(n_tri,n_tri);

    std::complex<double> j(0,1);
    double k_b = (std::sqrt(k2_b)).real();
    for(int mm = 0; mm < n_tri; mm++)
    {
        for(int nn = 0; nn < n_tri; nn++)
        {
            Eigen::VectorXd dxyz = (centroids.row(nn)-centroids.row(mm)).transpose();
            double dmn = std::sqrt(dxyz.transpose() * dxyz);

            std::complex<double> Gmn;
            double a_n = std::sqrt(areas(nn)/M_PI);
            if(mm==nn)
            {
                std::complex<double> H12 = boost::math::cyl_hankel_2(
                    1,
                    k_b*a_n
                );
                Gmn = -j/(2.0*k_b*k_b)*(M_PI*k_b*a_n*H12 - 2.0*j);
            }
            else
            {
                assert(dmn > a_n);
                // In order to apply integral rule, the distance must
                // be greater than the radius.
                std::complex<double> J1 = boost::math::cyl_bessel_j(
                    1,
                    k_b*a_n
                );
                std::complex<double> H02 = boost::math::cyl_hankel_2(
                    0,
                    k_b*dmn
                );
                Gmn = -j*M_PI*a_n*J1*H02/2.0/k_b;
            }
            // Background green's function is not symmetric for tris
            G(mm,nn) = Gmn;
        }
    }
}


Chamber::Chamber(const std::string& meshfile)
{
    mesh.buildTriangulation(meshfile);
    // Set the target_tot to be zero contrast.
    target_tot.eps_r.setLocations(mesh.centroids);
    target_tot.eps_r.fillOnes();
    target_tot.contrast.setLocations(mesh.centroids);
    target_tot.contrast.fillZeros();
    target_tot.ready = true;
}

void Chamber::A3P3(
    Eigen::MatrixXd& kpts,
    Eigen::VectorXcd& kvals,
    Eigen::VectorXcd& chi
){
    // Check for a mesh
    if(not(mesh.ready)){
        std::cerr << "Build a mesh before calling A3P3\n";
        assert(false);
    }

    // Check for antennas.
    if(not antennas.size()){
        std::cerr << "Set up antennas before calling A3P3\n";
        assert(false);
    }
    // Check for probes.
    if(not probes.size()){
        std::cerr << "Set up probes before calling A3P3\n";
        assert(false);
    }

    // Make scattered data from total data.
    if(not(Ez_tot_meas.size()==frequencies.size())){
        std::cerr << "Missing total field data\n";
        assert(false);
    } else{
        calcDataEzInc();
        Ez_sct_meas.resize(Ez_tot_meas.size());
        auto vf_idx{0};
        for(auto& vf:Ez_tot_meas){
            auto ff_idx{0};
            Ez_sct_meas[vf_idx].resize(vf.size());
            for(auto& ff:vf){
                Ez_sct_meas[vf_idx][ff_idx].setVals(
                    ff.getValRef()-Ez_inc_d[vf_idx][ff_idx].getValRef()
                );
                ff_idx++;
            }
            vf_idx++;
        }
    }

    auto n_data{0};
    for(auto& vec : Ez_sct_meas){
        assert(vec.size() == antennas.size());
        for(auto& field : vec){
            auto field_size{field.getValRef().size()};
            assert(field_size == probes.size());
            n_data += field_size;
        }
    }

    fillKSpace(
        kpts,
        kvals
    );


    Eigen::MatrixXi tri;
    delaunay(
        tri,
        kpts
    );
    tift(
        chi,
        tri,
        kpts,
        kvals,
        mesh.centroids.block(
            0,0,
            mesh.centroids.rows(),2
        )
    );
}

void Chamber::A3P4(
    std::string freqfile
){
    std::cout << "A3P4\n";
}

void Chamber::setTarget(const std::string& targetfile)
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

    target_tot.eps_r.fillOnes();
    Eigen::VectorXcd eps_to_set;
    target_tot.eps_r.getVals(eps_to_set);
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
                eps_to_set(rr) = eps_rel_complex;
            }
        }
    }
    target_tot.eps_r.setVals(eps_to_set);
    reader.close();
}

void Chamber::setupProbes(const std::string& probefile)
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
    probes.resize(0);
    while(reader)
    {
        Probe pi;
        reader >> x;
        if(reader.eof()) break;
        reader >> y;
        reader >> z;
        pi.location << x,y,z;
        pi.ready = true;
        probes.push_back(pi);
    }
    probe_points.resize(probes.size(),3);
    for(int pp = 0; pp < probes.size(); pp++)
    {
        probe_points.row(pp) = probes[pp].location.transpose();
    }
    reader.close();
}

void Chamber::readMeasuredData(
    std::string dataprefix,
    double noise_pct
)
{
    for(auto ifreq{0}; ifreq<frequencies.size();++ifreq){
        Eigen::MatrixXcd tot_matrix;
        const std::string datafilename{dataprefix+std::to_string(ifreq)+".txt"};
        ReadMatrixFromFile(
            datafilename,
            tot_matrix
        );
        const auto data_abs_avg{
            tot_matrix.array().abs().mean()
        };
        const Eigen::MatrixXd noise_angles{
            M_PI*
            Eigen::MatrixXd::Random(
                tot_matrix.rows(),
                tot_matrix.cols()
            )
        };
        const Eigen::MatrixXd noise_magnitudes{
            data_abs_avg*noise_pct/100.0*
            Eigen::MatrixXd::Random(
                tot_matrix.rows(),
                tot_matrix.cols()
            )
        };        
        Eigen::MatrixXcd noise(
            tot_matrix.rows(),
            tot_matrix.cols()
        );
        std::complex<double> imag_j(0,1);
        for(auto rr{0}; rr< tot_matrix.rows(); ++rr)
        {
            for(auto cc{0}; cc<tot_matrix.cols();++cc)
            {
                noise(rr,cc) = noise_magnitudes(rr,cc)*std::exp(imag_j*noise_angles(rr,cc));
            }
        }

        tot_matrix += noise;

        Ez_tot_meas[ifreq].resize(tot_matrix.cols());
        for(auto cc{0}; cc<tot_matrix.cols(); ++cc)
        {
            Eigen::VectorXcd column{tot_matrix.col(cc)};
            Ez_tot_meas[ifreq][cc].setLocations(probe_points);
            Ez_tot_meas[ifreq][cc].setVals(column);
        }
    }

}


void Chamber::A2Q3(
    Eigen::MatrixXcd & w_calc,
    Eigen::VectorXcd & X_calc
)
{
	// Annihilating the Annihilator

	// Basis is green's functions, centered at receivers

	// Build a matrix H s.t. Hij is the inner product of :
	// (g_i)* and g_j

	// For this, we need:
    // mesh
    // antennas
    // incident fields on domain
    // scattered fields at receivers
    // domain green's function

    if(!mesh.ready)
    {
        mesh.buildTriangulation(mesh.meshfilename);
    }

    if(!antennas.size())
    {
        std::cerr << "Hey, you didn't set up any antennas!" << std::endl;
    }
    assert(antennas.size());

    if(!Ez_inc[0].size())
    {
        calcDomainEzInc();
    }
    assert(Ez_inc[0].size());

    if(!Ez_sct_meas[0].size())
    {
        std::cerr << "We have no scattered data. Attempting to calculate..." << std::endl;
        if(!Ez_tot_meas[0].size())
        {
            std::cerr << "We have to total data either." << std::endl;
        }
        else
        {
            std::cerr << "Calculating Ez inc at probes" << std::endl;
            calcDataEzInc();
            std::cerr << "Got " << Ez_inc_d[0].size() << " fields" << std::endl;
            Ez_sct_meas[0].resize(Ez_inc_d[0].size());
            std::cerr << "Resized Ez_sct_meas[0] to " << Ez_sct_meas[0].size() << std::endl;
            std::cerr << "Gonna pull data from Ez_tot_meas[0](" << Ez_tot_meas[0].size() << "), and Ez_inc_d[0](" << Ez_inc_d[0].size() << ")" << std::endl;
            for(auto pp{0}; pp < Ez_inc_d[0].size();++pp)
            {
                const Eigen::VectorXcd & tot_p{Ez_tot_meas[0][pp].getValRef()};
                const Eigen::VectorXcd & inc_p{Ez_inc_d[0][pp].getValRef()};
                std::cout << "tot has " << tot_p.size() << " values " << std::endl;
                std::cout << "inc has " << inc_p.size() << " values " << std::endl;
                Eigen::VectorXcd sct_p = tot_p-inc_p;
                std::cout << "sct has " << sct_p.size() << " values " << std::endl;
                Ez_sct_meas[0][pp].setVals(sct_p);
            }
            assert(Ez_sct_meas[0].size());
        }
    }

    if(not (G_inc_domain_by_freq[0].rows() and G_inc_domain_by_freq[0].cols()))
    {
        std::cerr << "We have no green function operator. Attempting to build" << std::endl;
        assert(k2_bs[0] != 0.0);
        std::cerr << "frequencies[0] is " << frequencies[0] << std::endl;
        std::cerr << "k2_bs[0] is " << k2_bs[0] << std::endl;
        mesh.buildDomainGreen(
            G_inc_domain_by_freq[0],
            k2_bs[0]
        );
        std::cerr << "Built green function operator!" << std::endl;
    }
    // gg(:,k) is k2_bs[0]*G(r,r_k), where r_k are the probe locations.
	Eigen::MatrixXcd gg(mesh.centroids.rows(),probes.size());
    std::cerr << "Building " << gg.rows() << " by " << gg.cols() << " matrix of little g fields" << std::endl;
	for(auto pp{0};pp<probes.size();++pp){
		Eigen::VectorXcd gp;
        std::cerr << "Getting green function for probe " << pp << std::endl;
        // Set up a dummy antenna so I can getEz
        Antenna dummy_for_probe;
        dummy_for_probe.style = LineSource;
        dummy_for_probe.location = probes[pp].location;
        dummy_for_probe.coefficient = 1.0;
		dummy_for_probe.getEz(
            gp,
            mesh.centroids,
            frequencies[0]
        );
        std::cerr << "Got the hankel function field, multiplying by k2b = " << k2_bs[0] << std::endl;
		gp *= k2_bs[0];
        std::cerr << "Assigning to gg" << std::endl;
		gg.col(pp) = gp;
        std::cerr << "Assigned!" << std::endl;
	}

    if(!G_inc_data_by_freq[0].rows() || !G_inc_data_by_freq[0].cols())
    {
        std::cerr << "Annihilator needs G_inc_data_by_freq[0]. Building it now..." << std::endl;
        mesh.buildDataGreen(
            G_inc_data_by_freq[0],
            k2_bs[0],
            probe_points
        );
    }

    Eigen::MatrixXcd H{
        (k2_bs[0])*G_inc_data_by_freq[0]*(gg.conjugate())
    };

    std::cerr << "LU-decomposing H:" << std::endl;

	Eigen::PartialPivLU<Eigen::MatrixXcd> H_LU;

    H_LU.compute(H);

    std::cerr << "Decomposition complete! " << std::endl;
    std::cerr << "Decomp gave us a matrix of size " << H_LU.matrixLU().rows() << " by " << H_LU.matrixLU().cols() << std::endl;

    // WriteMatrixToFile(
    //     "H_LU.txt",
    //     H_LU.matrixLU()
    // );

    std::cerr << "Now let us make contrast sources..." << std::endl;
	Eigen::MatrixXcd w_basis_coeffs(
        probes.size(),
        antennas.size()
    );

    std::cerr << "Gonna get contrast sources for " << probes.size() << " probes " << std::endl;
    std::cerr << "How much scattered data do I have? "<< Ez_sct_meas[0].size() << std::endl;

	for(auto tt{0};tt<antennas.size();++tt){
        std::cerr << "Solving for contrast sources for tx " << tt << std::endl;
        Eigen::VectorXcd Ezs_t;
        Ez_sct_meas[0][tt].getVals(Ezs_t);
        std::cerr << "Asking H_LU to solve against " << Ezs_t.size() << " values " << std::endl;

        Eigen::VectorXcd w{H_LU.solve(Ezs_t)};
        std::cerr << "LU solve worked!!" << std::endl;

		w_basis_coeffs.col(tt) = w;
	}

    Eigen::MatrixXcd w_on_eles(
        mesh.areas.size(),
        antennas.size()
    );

    for(auto tt{0}; tt<antennas.size();++tt)
    {
        Eigen::VectorXcd w_for_t(mesh.areas.size());
        for(auto bb{0}; bb<w_basis_coeffs.rows();++bb)
        {
            w_for_t += w_basis_coeffs(bb,tt)*(gg.col(bb).conjugate());
        }
        w_on_eles.col(tt) = w_for_t;
    }

    Eigen::MatrixXcd domain_total_fields{
        mesh.areas.size(),
        antennas.size()
    };

    for(auto tt{0}; tt<w_on_eles.cols(); ++tt)
    {
        // My G operators map from k2_bs[0]*w to u^s
        domain_total_fields.col(tt) = Ez_inc[0][tt].getValRef() + k2_bs[0]*G_inc_domain_by_freq[0]*w_on_eles.col(tt);
    }
    Eigen::VectorXcd contrast(mesh.areas.size());
    for(auto kk{0}; kk<contrast.size(); ++kk)
    {
        std::complex<double> numerator{0},denominator{0};
        for(auto tt{0}; tt<antennas.size(); ++tt)
        {
            std::complex<double> ukt{domain_total_fields(kk,tt)};
            std::complex<double> uktstar{std::conj(ukt)};
            numerator += uktstar*w_on_eles(kk,tt);
            denominator += uktstar*ukt;
        }
        contrast(kk) = numerator/denominator;
    }
    w_calc = w_on_eles;
    X_calc = contrast;
}

static int findTikhLeftestTurn(
    const std::vector<Eigen::Vector2d> & curve
)
{
    std::size_t curve_len{curve.size()};
    int turn_index = -1;
    // Need at least five data points to find a maximally left turn.
    if(5 <= curve_len)
    {
        std::cerr << "Curve has " << curve_len << " points, so start searching..." << std::endl;
        std::vector<Eigen::Vector2d> arrows(curve_len-1);
        for(auto cc{0};cc<(curve_len-1);++cc)
        {
            // std::cerr << "Point: " << cc << std::endl;
            // std::cerr << curve[cc] << std::endl;
            arrows[cc] = curve[cc+1] - curve[cc];
            // std::cerr << "Arrow: " << cc << std::endl;
            // std::cerr << arrows[cc] << std::endl;
        }
        std::vector<double> thetas(arrows.size()-1);
        for(auto aa{0};aa<(arrows.size()-1);++aa)
        {
            Eigen::Vector2d u{arrows[aa]},v{arrows[aa+1]};
            double z_of_cross{
                u[0]*v[1]-u[1]*v[0]
            };
            // std::cerr << "z_of_cross = " << z_of_cross << std::endl;
            double sintheta{
                z_of_cross/u.norm()/v.norm()
            };
            // std::cerr << "sintheta = " << sintheta << std::endl;
            double theta{
                std::asin(sintheta)
            };
            // std::cerr << "theta = " << theta << std::endl;
            thetas[aa] = theta;
        }
        for(auto tt{0}; tt<(thetas.size()-2); ++tt)
        {
            double tp0{thetas[tt]},tp1{thetas[tt+1]},tp2{thetas[tt+2]};
            std::cerr << tp0 << ",";
            if((tp0 < tp1) and (tp2 < tp1) and (0 < tp1))
            {
                // Found a peak!
                std::cerr << tp1 << "," << tp2 << ",found peak";
                //std::cerr << "FOUND PEAK AT ANGLE INDEX " << tt+1 << std::endl;
                turn_index = tt+2;
                break;
            }
        }
        std::cerr << std::endl;
    }
    else
    {
        std::cerr << "Curve has " << curve_len << " points, so no bother " << std::endl;
    }

    return(turn_index);
}

void Chamber::A2Q5(
    Eigen::MatrixXcd & w_calc,
    Eigen::MatrixXcd & u_calc,
    Eigen::VectorXcd & X_calc,
    std::vector<std::vector<Eigen::Vector2d> > & curves,
    bool tikhonov
)
{
    auto ntx{antennas.size()};
    auto ntri{mesh.areas.size()};

    Ez_sct_meas[0].resize(ntx);

    if(!Ez_inc_d[0].size())
    {
        calcDataEzInc();
    }

    Eigen::MatrixXcd Ez_sct_meas_mat(G_inc_data_by_freq[0].rows(),ntx);
    for(auto tt{0}; tt<ntx; ++tt)
    {
        Eigen::VectorXcd tot_t;
        Ez_tot_meas[0][tt].getVals(tot_t);
        Eigen::VectorXcd inc_t;
        Ez_inc_d[0][tt].getVals(inc_t);
        Eigen::VectorXcd sct_t{tot_t-inc_t};
        Ez_sct_meas[0][tt].setVals(sct_t);
        Ez_sct_meas_mat.col(tt) = sct_t;
    }


    if(!G_inc_data_by_freq[0].rows() || !G_inc_data_by_freq[0].cols())
    {
        std::cerr << "Annihilator needs G_inc_data_by_freq[0]. Building it now..." << std::endl;
        mesh.buildDataGreen(
            G_inc_data_by_freq[0],
            k2_bs[0],
            probe_points
        );
    }
    Eigen::MatrixXcd P{k2_bs[0]*G_inc_data_by_freq[0]};
    Eigen::MatrixXcd optimal_alphas(P.cols(),ntx);

    if(tikhonov)
    {
        Eigen::MatrixXcd PhP{
            P.adjoint()*P
        };
        Eigen::MatrixXcd I{
            Eigen::MatrixXcd::Identity(
                PhP.rows(),PhP.cols()
            )
        };

        std::cerr << "PhP is " << PhP.rows() << " by " << PhP.cols() << std::endl;

        std::vector<bool> optimum_found_for_tx(ntx,false);
        // std::vector<std::vector<double> > tikh_distances_for_tx(ntx);
        curves.resize(antennas.size());

        Eigen::MatrixXcd Phd_all{
            P.adjoint()*Ez_sct_meas_mat
        };
        std::vector<double> optimal_lambdas(ntx);

        double lambda_exp{0.0};
        double lambda0{std::pow(10,+lambda_exp)};

        Eigen::MatrixXcd A0{PhP+lambda0*I};
        Eigen::PartialPivLU<Eigen::MatrixXcd> LU_0;
        LU_0.compute(A0);
        for(auto tt{0}; tt<antennas.size(); ++tt)
        {
            Eigen::VectorXcd solution{LU_0.solve(Phd_all.col(tt))};
            Eigen::VectorXcd residual{P*solution - Ez_sct_meas_mat.col(tt)};

            std::complex<double> complex_sol_norm_sq{solution.adjoint()*solution};
            std::complex<double> complex_res_norm_sq{residual.adjoint()*residual};

            double sol_norm_sq{complex_sol_norm_sq.real()};
            double res_norm_sq{complex_res_norm_sq.real()};

            double log_sol_norm{std::log(std::sqrt(sol_norm_sq))};
            double log_res_norm{std::log(std::sqrt(res_norm_sq))};

            Eigen::Vector2d tikh_point;
            tikh_point(0) = log_res_norm;
            tikh_point(1) = log_sol_norm;

            curves[tt].resize(1);
            curves[tt][0] = tikh_point;
            std::cerr << "Initial tikh point for tx " << tt << "is: " << std::endl;
            std::cerr << curves[tt][0] << std::endl;
        }
        while(
            !std::all_of(
                optimum_found_for_tx.begin(),
                optimum_found_for_tx.end(),
                [](bool x){return x;}
            )
        )
        {
            lambda_exp++;
            std::cerr << "lambda_exp = " << lambda_exp << std::endl;

            double lambda_big{std::pow(10,+lambda_exp)};
            double lambda_sml{std::pow(10,-lambda_exp)};

            std::cerr << "big = " << lambda_big << " and small = " << lambda_sml << std::endl;

            Eigen::MatrixXcd Abig=PhP+lambda_big*I;
            Eigen::MatrixXcd Asml=PhP+lambda_sml*I;

            std::cerr << "First elements of Abig, Asml:" << Abig(0,0) << "," << Asml(0,0) << std::endl;

            Eigen::PartialPivLU<Eigen::MatrixXcd> LU_big;
            Eigen::PartialPivLU<Eigen::MatrixXcd> LU_sml;

            std::cerr << "Computing LU factorization of A-matrices...(" << Abig.rows() << "," << Abig.cols() << ")";

            LU_big.compute(Abig);
            std::cerr << " done 1...";
            LU_sml.compute(Asml);
            std::cerr << " done 2!" << std::endl;

            for(int itx = 0; itx < ntx; itx++)
            {

                if(!optimum_found_for_tx[itx])
                {
                    std::cerr << "Optimum not yet found for tx " << itx << std::endl;
                    // Solve system
                    Eigen::VectorXcd alpha_big{LU_big.solve(Phd_all.col(itx))};
                    Eigen::VectorXcd alpha_sml{LU_sml.solve(Phd_all.col(itx))};

                    double res_norm_big{(Ez_sct_meas_mat.col(itx)-P*alpha_big).norm()};
                    double alpha_norm_big{alpha_big.norm()};

                    double res_norm_sml{(Ez_sct_meas_mat.col(itx)-P*alpha_sml).norm()};
                    double alpha_norm_sml{alpha_sml.norm()};

                    Eigen::Vector2d tikh_point_big,tikh_point_sml;

                    tikh_point_big[0] = std::log(res_norm_big);
                    tikh_point_big[1] = std::log(alpha_norm_big);

                    tikh_point_sml[0] = std::log(res_norm_sml);
                    tikh_point_sml[1] = std::log(alpha_norm_sml);

                    std::cerr << "Prepending:" << std::endl;
                    std::cerr << tikh_point_sml << std::endl;

                    curves[itx].insert(
                        curves[itx].begin(),
                        tikh_point_sml
                    );

                    std::cerr << "Appending" << std::endl;
                    std::cerr << tikh_point_big << std::endl;

                    curves[itx].push_back(
                        tikh_point_big
                    );

                    // Scan this curve for an optimum(scan for a peak in left-curviness)
                    int turn_index = findTikhLeftestTurn(
                        curves[itx]
                    );
                    if(0<turn_index)
                    {

                        double lambda_opt{
                            lambda_sml*std::pow(10,turn_index)
                        };
                        std::cerr << "Found optimal lambda for tx " << itx << ": " << lambda_opt << std::endl;
                        optimal_lambdas[itx] = lambda_opt;
                        optimum_found_for_tx[itx] = true;
                    }
                }
            }
        }

        std::cerr << "Moving on from lambda calc..." << std::endl;
        // Cool, we got all of our optimal lambdas!
        // Calculate the alphas which correspond to all these optimal
        // values of lambda. This could possibly be done in the loops up
        // above, but I don't want to spend the energy on that.
        // Recalculating them here is simpler and cleaner.
        std::vector<bool> alpha_is_calculated(optimal_lambdas.size(),false);
        for(auto ll{0}; ll < optimal_lambdas.size(); ++ll)
        {
            if(alpha_is_calculated[ll]) continue;
            std::cerr << "Calculating optimal alphas for tx " << ll << std::endl;
            double lambda_opt{optimal_lambdas[ll]};
            std::cerr << "\tOptimal lambda here is " << lambda_opt << std::endl;
            std::cerr << "\tFactor(PhP+lambdaI)..." << std::endl;
            Eigen::PartialPivLU<Eigen::MatrixXcd> LU_opt;
            std::cerr << "\tInitialized" << std::endl;
            LU_opt.compute(PhP+lambda_opt*I);
            std::cerr << "\tFactored" << std::endl;

            for(auto ii{ll};ii<optimal_lambdas.size();++ii)
            {
    	    // Check which alpha-systems use this lambda for regularization.
    	    bool system_uses_this_lambda{
    	        std::abs((optimal_lambdas[ii]-lambda_opt)/lambda_opt) < 0.01
    	    };
                if(system_uses_this_lambda && !alpha_is_calculated[ii])
                {
                    std::cerr << "\t\tSystem " <<ii << " also uses this lambda" << std::endl;
                    optimal_alphas.col(ii) = LU_opt.solve(Phd_all.col(ii));
                    alpha_is_calculated[ii] = true;
                }
            }


        }
    }
    else
    {
        if(P.rows()<P.cols())
        {
            Eigen::MatrixXcd Ph{P.adjoint()};
            Eigen::MatrixXcd PPh{P*Ph};
            Eigen::PartialPivLU<Eigen::MatrixXcd> PPh_LU;
            PPh_LU.compute(PPh);
            for(auto tt{0}; tt<ntx; ++tt)
            {
                optimal_alphas.col(tt) = Ph*(PPh_LU.solve(Ez_sct_meas_mat.col(tt)));
            }
        }
        else
        {
            Eigen::MatrixXcd Phd_all{P.adjoint()*Ez_sct_meas_mat};
            Eigen::MatrixXcd PhP{P.adjoint()*P};
            Eigen::PartialPivLU<Eigen::MatrixXcd> PhP_LU;
            PhP_LU.compute(PhP);

            for(auto tt{0}; tt<ntx; ++tt)
            {
                optimal_alphas.col(tt) = PhP_LU.solve(Phd_all.col(tt));
            }

        }
    }
    if(!Ez_inc[0].size())
    {
        calcDomainEzInc();
    }
    if(!G_inc_domain_by_freq[0].rows() || !G_inc_domain_by_freq[0].cols())
    {
        std::cerr << "Annihilator need G_inc_domain_by_freq[0]. Building it now..." << std::endl;
        mesh.buildDomainGreen(
            G_inc_domain_by_freq[0],
            k2_bs[0]
        );
    }
    std::cerr << "Making the matrix of total fields in imaging domain..." << std::endl;
    Eigen::MatrixXcd Ez_tot_opt(ntri,ntx);

    for(auto cc{0}; cc<Ez_tot_opt.cols(); ++cc)
    {
        Eigen::VectorXcd inc_c;
        Ez_inc[0][cc].getVals(inc_c);
        Ez_tot_opt.col(cc) = inc_c + k2_bs[0]*G_inc_domain_by_freq[0]*optimal_alphas.col(cc);
    }
    std::cerr << "success!" << std::endl;

    Eigen::VectorXcd chi_opt(ntri);
    chi_opt.setZero();

    for(auto ii{0}; ii < ntri; ++ii)
    {
        std::complex<double> numerator{0.0},denominator{0.0};
        for(auto tt{0}; tt<ntx; ++tt)
        {
            std::complex<double> uit{Ez_tot_opt(ii,tt)};
            std::complex<double> uitstar{std::conj(uit)};
            std::complex<double> ait{optimal_alphas(ii,tt)};
            numerator += uitstar*ait;
            denominator += uitstar*uit;
        }
        chi_opt(ii) = numerator/denominator;

    }

    w_calc = optimal_alphas;
    u_calc = Ez_tot_opt;
    X_calc = chi_opt;

}

void Chamber::setupAntennas(const std::string& antennafile)
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
        if(sty_string.compare("planewave") == 0)
        {
            iant.style = PlaneWave;
            iant.direction << x,y,z;
            iant.direction /= std::sqrt(
                (iant.direction.transpose()*iant.direction)
            );
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
        iant.ready = true;
        antennas.push_back(iant);
    }
    reader.close();
}

void Chamber::setFrequencies(
    const std::string& freqfile
)
{
    std::ifstream reader;
    reader.open(freqfile, std::ifstream::in);
    frequencies.resize(0);

    while(true){
        double freq;
        reader >> freq;
        if(reader.eof()){
            break;
        }
        frequencies.push_back(freq);
        double omega = 2*M_PI*freq;
        k2_bs.push_back(omega*omega/CNAUGHT/CNAUGHT);
    }

    G_inc_domain_by_freq.resize(frequencies.size());
    G_inc_data_by_freq.resize(frequencies.size());
    L_domain_by_freq.resize(frequencies.size());

    Ez_inc.resize(frequencies.size());
    Ez_sct.resize(frequencies.size());
    Ez_tot.resize(frequencies.size());

    Ez_inc_d.resize(frequencies.size());
    Ez_sct_d.resize(frequencies.size());
    Ez_tot_d.resize(frequencies.size());

    Ez_sct_meas.resize(frequencies.size());
    Ez_tot_meas.resize(frequencies.size());

    reader.close();
}

void Chamber::calcDomainEzInc(void)
{
    Ez_inc.resize(0);
    for(auto& freq : frequencies){
        std::vector<Field> vec_of_field_for_freq;
        for(auto& ant : antennas){
            Field field_for_ant;
            Eigen::VectorXcd ez_for_ant;
            ant.getEz(
                ez_for_ant,
                mesh.centroids,
                freq
            );
            field_for_ant.setLocations(mesh.centroids);
            field_for_ant.setVals(ez_for_ant);
            vec_of_field_for_freq.push_back(field_for_ant);
        }
        Ez_inc.push_back(vec_of_field_for_freq);
    }
}

void Chamber::calcDataEzInc(void)
{
    Ez_inc_d.resize(0);

    for(auto& freq : frequencies){
        std::vector<Field> vec_of_field_for_freq;
        for(auto& ant : antennas){
            Field field_for_ant;
            Eigen::VectorXcd ez_for_ant;
            ant.getEz(
                ez_for_ant,
                probe_points,
                freq
            );
            field_for_ant.setLocations(probe_points);
            field_for_ant.setVals(ez_for_ant);
            vec_of_field_for_freq.push_back(field_for_ant);
        }
        Ez_inc_d.push_back(vec_of_field_for_freq);
    }
}

void Chamber::calcDomainEzTot(void)
{
    // Calculate incident fields.
    //Ez_inc.resize(mesh.areas.size(),antennas.size());
    calcDomainEzInc();

    for(auto ifreq{0};ifreq<frequencies.size();++ifreq){
        // Build domain green matrix
        std::cerr << "Building domain green in calcDomainEzTot..." << std::endl;
        mesh.buildDomainGreen(
            G_inc_domain_by_freq[ifreq],
            k2_bs[ifreq]
        );
        std::cerr << "G:(" << G_inc_domain_by_freq[ifreq].rows() << ",";
        std::cerr << G_inc_domain_by_freq[ifreq].cols() << ")" << std::endl;
        std::cerr << "Building DeltaK2..." << std::endl;

        // Build contrast matrix
        Eigen::VectorXcd eps_r_from_target;
        target_tot.eps_r.getVals(eps_r_from_target);
        Eigen::MatrixXcd DeltaK2;
        DeltaK2.resize(eps_r_from_target.size(),eps_r_from_target.size());
        DeltaK2.setZero();
        for(int dd = 0; dd < DeltaK2.cols();dd++)
        {
            DeltaK2(dd,dd) = k2_bs[ifreq]*eps_r_from_target(dd)-k2_bs[ifreq];
        }

        std::cerr << "Building L_domain" << std::endl;
        L_domain_by_freq[ifreq] = -G_inc_domain_by_freq[ifreq];

        for(int cc = 0; cc< DeltaK2.cols(); cc++)
        {
            L_domain_by_freq[ifreq].col(cc) *= DeltaK2(cc,cc);
            L_domain_by_freq[ifreq](cc,cc) += 1.0;
        }

        std::cerr << "filled,";
        // Perform LU factorization of domain L operator
        Eigen::FullPivLU<Eigen::MatrixXcd> LU_L;
        LU_L.compute(L_domain_by_freq[ifreq]);
        std::cerr << "factored" << std::endl;

        Ez_sct[ifreq].resize(antennas.size());
        for(auto ss{0}; ss < Ez_sct[ifreq].size(); ss++)
        {
            Ez_sct[ifreq][ss].setLocations(mesh.centroids);
        }
        // Calculate all the scattered fields
        for(int jj = 0; jj < Ez_inc[ifreq].size(); jj++)
        {
            Eigen::VectorXcd inc_vec,sct_vec;
            Ez_inc[ifreq][jj].getVals(inc_vec);
            sct_vec = LU_L.solve(
                G_inc_domain_by_freq[ifreq]*(
                    DeltaK2*inc_vec
                )
            );
            Ez_sct[ifreq][jj].setVals(sct_vec);
        }
        // Add scattered to indident to get total fields
        Ez_tot[ifreq].resize(antennas.size());
        for(auto tt{0}; tt<antennas.size(); ++tt)
        {
            Eigen::VectorXcd inc_vec,sct_vec,tot_vec;
            Ez_inc[ifreq][tt].getVals(inc_vec);
            Ez_sct[ifreq][tt].getVals(sct_vec);
            tot_vec = inc_vec+sct_vec;
            Ez_tot[ifreq][tt].setVals(tot_vec);
        }
    }
}

void Chamber::calcDataEzTot(void)
{
    for(auto ifreq{0};ifreq<frequencies.size();++ifreq){
        if(!(Ez_tot[ifreq].size() == antennas.size()))
        {
            std::cerr << "Ez_tot was not built. Getting it now.";
            std::cerr << std::endl;
            calcDomainEzTot();
        }
        if(!G_inc_data_by_freq[ifreq].rows() || !G_inc_data_by_freq[ifreq].cols())
        {
            std::cerr << "Data green was not built. Getting it now";
            std::cerr << std::endl;
            mesh.buildDataGreen(
                G_inc_data_by_freq[ifreq],
                k2_bs[ifreq],
                probe_points
            );
        }
        std::cerr << "Resizing d-field holders to ";
        std::cerr << probe_points.rows() << " by ";
        std::cerr << antennas.size() << std::endl;

        calcDataEzInc();
        Ez_sct_d[ifreq].resize(antennas.size());
        Ez_tot_d[ifreq].resize(antennas.size());

        Eigen::MatrixXcd DeltaK2;
        const Eigen::VectorXcd& target_epsr{target_tot.eps_r.getValRef()};
        auto n_eps{target_epsr.size()};
        DeltaK2.resize(n_eps,n_eps);

        for(auto dd{0}; dd<n_eps; ++dd){
            DeltaK2(dd,dd) = k2_bs[ifreq]*(target_epsr(dd)-1.0);
        }


        for(int tt = 0; tt < Ez_tot[ifreq].size(); tt++)
        {
            std::cerr << "Calculating d-sct for antenna ";
            std::cerr << tt << "\n";
            Eigen::VectorXcd Ez_tot_t;
            Ez_tot[ifreq][tt].getVals(Ez_tot_t);
            Eigen::VectorXcd Ez_sct_d_t{G_inc_data_by_freq[ifreq]*(DeltaK2*(Ez_tot_t))};
            Ez_sct_d[ifreq][tt].setVals(Ez_sct_d_t);
            Eigen::VectorXcd Ez_inc_d_t;
            Ez_inc_d[ifreq][tt].getVals(Ez_inc_d_t);
            Eigen::VectorXcd Ez_tot_d_t{Ez_sct_d_t+Ez_inc_d_t};
            Ez_tot_d[ifreq][tt].setVals(Ez_tot_d_t);
        }

    }

}

void Chamber::fillKSpace(
    Eigen::MatrixXd& kpts,
    Eigen::VectorXcd& kvals
){

    // Every measured field datum produces a point in k-space.
    assert(frequencies.size());
    assert(frequencies.size() == Ez_sct_meas.size());

    auto n_data{0};
    for(auto& vec : Ez_sct_meas){
        assert(vec.size() == antennas.size());
        for(auto& field : vec){
            auto field_size{field.getValRef().size()};
            assert(field_size == probes.size());
            n_data += field_size;
        }
    }

    kpts.resize(n_data,2);
    kvals.resize(n_data);
    std::complex<double> j_imag(0.0,1.0);
    auto k_ptr{0};
    for(auto ifreq{0}; ifreq<frequencies.size();++ifreq){
        const double kf{std::sqrt(k2_bs[ifreq]).real()};
        for(auto itx{0}; itx<Ez_sct_meas[ifreq].size(); ++itx){
            auto& field_for_tx{Ez_sct_meas[ifreq][itx].getValRef()};
            auto field_size{field_for_tx.size()};
            Eigen::Vector3d rt;
            if(antennas[itx].style == LineSource){
                rt = antennas[itx].location;
            } else if(antennas[itx].style==PlaneWave){
                rt = antennas[itx].direction;
            } else {
                assert(false);
            }

            const auto rt_norm{rt.norm()};
            Eigen::Vector3d rt_hat{rt/rt_norm};
            for(auto irx{0}; irx<field_size;++irx){
                Eigen::Vector3d rk{probes[irx].location};
                const auto rk_norm{rk.norm()};
                Eigen::Vector3d rk_hat{rk/rk_norm};
                Eigen::Vector3d qq;
                std::complex<double> cc;
                if(antennas[itx].style == LineSource){
                    cc = antennas[itx].coefficient*
                        -j_imag*kf*std::exp(
                            -j_imag*kf*(
                                rk_norm+rt_norm
                            )
                        )/(
                            8.0*M_PI*std::sqrt(
                                rk_norm*rt_norm
                            )
                        );
                    qq = -rk_hat-rt_hat;
                } else if(antennas[itx].style == PlaneWave){
                    cc = antennas[itx].coefficient*
                        kf*kf*(1.0-j_imag)*std::exp(
                            -j_imag*kf*rk_norm
                        )/(
                            4*std::sqrt(
                                M_PI*kf*rk_norm
                            )
                        );
                    qq = rt_hat-rk_hat;
                } else {
                    assert(false);
                }
                auto kx{kf*qq(0)};
                auto ky{kf*qq(1)};
                auto kz{field_for_tx(irx)/cc};
                kpts(k_ptr,0) = kx;
                kpts(k_ptr,1) = ky;
                kvals(k_ptr)  = kz;
                k_ptr++;
            }
        }
    }
}

void Chamber::setupTx2RxMap(
    const std::string& mapfile
)
{
    std::ifstream reader;
    reader.open(mapfile, std::ifstream::in);

    std::string ins;
    reader >> ins;

    while(!reader.eof()){
        assert(ins == "FrequencyIndex");
        int ff_index;
        reader >> ff_index;
        ff_index--;
        assert(0<=ff_index);
        tx2rx.resize(ff_index+1);

        reader >> ins;
        assert(ins == "TransmitterIndex");

        int tx_index;
        reader >> tx_index;
        tx_index--;
        assert(0<=tx_index);
        tx2rx[ff_index].resize(tx_index+1);

        reader >> ins;
        assert(ins == "ReceiverIndices");

        while(true){
            reader >> ins;
            if(ins == "FrequencyIndex"){
                break;
            } else if(reader.eof()){
                break;
            } else{
                int rx_index{std::stoi(ins)};
                rx_index--;
                assert(0<=rx_index);
                tx2rx[ff_index][tx_index].push_back(rx_index);
            }
        }
    }
    reader.close();
}

void DomainTotalFieldSolve(
    const std::vector<std::vector<Field> >& U_inc,
    const Target& target,
    const std::vector<Eigen::MatrixXcd>& DomainGreens,
    std::vector<Eigen::PartialPivLU<Eigen::MatrixXcd> >& LU_Ls,
    bool build_Ls,
    const std::vector<std::complex<double> >& k2s,
    std::vector<std::vector<Field> >& U_tot
){
    if(build_Ls){
        const Eigen::VectorXcd& contrast_vec = target.contrast.getValRef();
        LU_Ls.resize(DomainGreens.size());
        for(auto ifreq{0};ifreq<DomainGreens.size();++ifreq){
            const Eigen::MatrixXcd L{
                Eigen::MatrixXcd::Identity(
                    DomainGreens[ifreq].rows(),
                    DomainGreens[ifreq].cols()
                )
                -k2s[ifreq]*
                DomainGreens[ifreq]*
                contrast_vec.asDiagonal()
            };
            // L += Eigen::MatrixXcd::Identity(L.rows(),L.cols());
            LU_Ls[ifreq].compute(L);
        }
    }
    U_tot.resize(U_inc.size());
    for(auto ifreq{0}; ifreq<U_inc.size();++ifreq){
        U_tot[ifreq].resize(U_inc[ifreq].size());
        for(auto itx{0}; itx<U_inc[ifreq].size();++itx){
            U_tot[ifreq][itx].setVals(
                LU_Ls[ifreq].solve(
                    U_inc[ifreq][itx].getValRef()
                )
            );
        }
    }
}

void DataScatteredFieldSolve(
    const std::vector<std::vector<Field> >& U_tot_dom,
    const Target& target,
    const std::vector<Eigen::MatrixXcd>& DataGreens,
    const std::vector<std::complex<double> >&k2s,
    std::vector<std::vector<Field> >& U_sct_dat
){
    U_sct_dat.resize(U_tot_dom.size());
    for(auto ifreq{0};ifreq<U_tot_dom.size();++ifreq){
        U_sct_dat[ifreq].resize(U_tot_dom[ifreq].size());
        for(auto itx{0}; itx<U_tot_dom[ifreq].size();++itx){
            U_sct_dat[ifreq][itx].setVals(
                k2s[ifreq]*(
                    DataGreens[ifreq]*(
                        target.contrast.getValRef().asDiagonal()*U_tot_dom[ifreq][itx].getValRef()
                    )
                )
            );
        }
    }
}

DBIMInversion Chamber::distortedBornIterativeMethod(){

    assert(frequencies.size());
    assert(probes.size());
    assert(antennas.size());
    assert(tx2rx.size() == frequencies.size());
    assert(Ez_tot_meas.size());

    G_inc_domain_by_freq.resize(frequencies.size());

    for(auto ifreq{0};ifreq<frequencies.size();ifreq++){
        mesh.buildDataGreen(
            G_inc_data_by_freq[ifreq],
            k2_bs[ifreq],
            probe_points
        );
        mesh.buildDomainGreen(
            G_inc_domain_by_freq[ifreq],
            k2_bs[ifreq]
        );
    }

    // Build all my Ms matrices
    M_s_data.resize(frequencies.size());
    for(auto& M_s_vec : M_s_data){
        M_s_vec.resize(antennas.size());
        for(auto& M_s : M_s_vec){
            M_s.resize(
                probes.size(),
                probes.size()
            );
        }
    }

    calcDataEzInc();
    calcDomainEzInc();

    for(auto ifreq{0}; ifreq<tx2rx.size(); ++ifreq){
        for(auto itx{0}; itx < tx2rx[ifreq].size(); ++itx){
            auto row_ptr{0};
            M_s_data[ifreq][itx].resize(
                tx2rx[ifreq][itx].size(),
                probes.size()
            );
            M_s_data[ifreq][itx].setZero();
            for(auto& rx_idx : tx2rx[ifreq][itx]){
                M_s_data[ifreq][itx](row_ptr++,rx_idx)=1.0;
            }
        }
    }

    // double sing_val_thresh{-1};
    // while(not(0.0<sing_val_thresh && sing_val_thresh<1.0)){
    //     std::cout << "What threshold shall I use for singular values?[0-1]\n";
    //     std::cin >> sing_val_thresh;
    // }
    int n_dbim_iterations{-1};
    while(n_dbim_iterations<0){
        std::cout << "How many DBIM iterations shall we run?\n";
        std::cin >> n_dbim_iterations;
    }
    int n_lscg_iterations{-1};
    while(n_lscg_iterations<0){
        std::cout << "How many least-squares conjugate-gradient steps?\n";
        std::cin >> n_lscg_iterations;
    }

    double lscg_tolerance{-1.0};
    while(not(0.0<lscg_tolerance and lscg_tolerance < 1.0)){
        std::cout << "Least-squares conjugate-gradient tolerance?\n";
        std::cin >> lscg_tolerance;
    }

    double contrast_update_scale;
    std::cout << "Contrast update scale?\n";
    std::cin >> contrast_update_scale;
    DBIMInversion inv_log;
    inv_log.Ez_tot_meas = Ez_tot_meas;
    inv_log.imaging_mesh = mesh;

    for(
        auto dbim_iter{0};
        dbim_iter<n_dbim_iterations;
        ++dbim_iter
    ){
        // Calculate background fields on the domain.
        //  Use the background contrast as the target.
        //  Solve using the incident green's function and the incident fields.
        //  Ask DomainTotalFieldSolve to please give us the factored L matrices
        std::cout << "DBIM Iteration " << dbim_iter << "\n";
        std::vector<Eigen::PartialPivLU<Eigen::MatrixXcd> > LU_L_inc_domain_by_freq;
        std::vector<std::vector<Field> > Ez_bkg_dom;
        std::cout << "\tCalculating background fields in domain for "
         << frequencies.size() << " freqs, "
         << antennas.size() << " transmitters\n";
        DomainTotalFieldSolve(
            Ez_inc,
            target_tot,
            G_inc_domain_by_freq,
            LU_L_inc_domain_by_freq,
            true,
            k2_bs,
            Ez_bkg_dom
        );

        // DataScatteredFieldSolve to propagate those total fields to scattered
        // fields at receiver points.
        // Add the analytic incident field to get background field at receiver
        // points.
        std::cout << "\tCalculating scattered fields at probe locations, due to background field in domain.\n";
        std::vector<std::vector<Field> > Ez_bkg_sct_dat;
        DataScatteredFieldSolve(
            Ez_bkg_dom,
            target_tot,
            G_inc_data_by_freq,
            k2_bs,
            Ez_bkg_sct_dat
        );

        std::cout << "\tAdding incident field to scattered field to get background field at probe locations\n";
        std::vector<std::vector<Field> > Ez_bkg_dat(Ez_bkg_sct_dat.size());
        for(auto ifreq{0};ifreq<Ez_bkg_sct_dat.size();++ifreq){
            auto new_size{Ez_bkg_sct_dat[ifreq].size()};
            Ez_bkg_dat[ifreq].resize(new_size);
            for(auto itx{0};itx < Ez_bkg_sct_dat[ifreq].size();++itx){
                Ez_bkg_dat[ifreq][itx].setVals(
                    Ez_inc_d[ifreq][itx].getValRef()+
                    Ez_bkg_sct_dat[ifreq][itx].getValRef()
                );
            }
        }

        // Form G^dat_bkg
        //  Impose incident fields located at the probe points

        //  Set up some dummy antennas located at the probe points, and ask them
        //  for their Ez responses.
        //  Set up a target which is the background contrast.
        //  DomainTotalFieldSolve with those incident fields and that target.
        //  Pass in those handy-dandy factored L-matrices.
        std::cout << "\tCreating incident linesources at probe locations.\n";
        std::vector<std::vector<Field> > probe_unit_ez_inc(frequencies.size());
        for(auto ifreq{0}; ifreq<frequencies.size();++ifreq){
            probe_unit_ez_inc[ifreq].resize(probes.size());
            for(auto pp{0}; pp<probes.size();++pp){
                Antenna dummy_antenna;
                dummy_antenna.style = LineSource;
                dummy_antenna.location = probes[pp].location;
                dummy_antenna.coefficient = 1.0;
                Eigen::VectorXcd dummy_ez_inc;
                dummy_antenna.getEz(
                    dummy_ez_inc,
                    mesh.centroids,
                    frequencies[ifreq]
                );
                probe_unit_ez_inc[ifreq][pp].setVals(dummy_ez_inc);
            }
        }

        std::vector<std::vector<Field> > Ez_tot_due_to_probe_units;

        std::cout << "\tCalculating total field response to linesources at probe locations.\n";
        DomainTotalFieldSolve(
            probe_unit_ez_inc,
            target_tot,
            G_inc_domain_by_freq,
            LU_L_inc_domain_by_freq,
            false,
            k2_bs,
            Ez_tot_due_to_probe_units
        );

        std::cout << "\tBuilding data Green's function from total field responses and element areas.\n";
        std::vector<Eigen::MatrixXcd> G_bkg_dat_by_freq(frequencies.size());
        for(auto ifreq{0}; ifreq< G_bkg_dat_by_freq.size();++ifreq){
            G_bkg_dat_by_freq[ifreq].resize(
                probes.size(),
                mesh.centroids.rows()
            );
            for(auto iprobe{0}; iprobe<probes.size();++iprobe){
                G_bkg_dat_by_freq[ifreq].row(iprobe) =(
                    Ez_tot_due_to_probe_units[ifreq][iprobe].getValRef().array()*
                    mesh.areas.array()
                ).matrix().transpose();
            }
            // std::cout << G_bkg_dat_by_freq[ifreq].row(0) << "\n\n\n";
            // std::cout << G_inc_data_by_freq[ifreq].row(0) << "\n\n\n";

        }

        if(dbim_iter==0){
            std::cout << "\tIt is the first time through the loop, so I can use the background analytic data Green's function\n";
            G_bkg_dat_by_freq = G_inc_data_by_freq;
        }


        // Build big vector of Utotmeas, and Ubkg
        std::cout << "\tCounting the number of data. ";
        auto n_data{0};
        for(auto ifreq{0}; ifreq<Ez_tot_meas.size();++ifreq){
            for(auto itx{0};itx<Ez_tot_meas[ifreq].size();++itx){
                n_data += M_s_data[ifreq][itx].rows();
            }
        }
        auto max_possible_n_data{
            frequencies.size()*
            antennas.size()*
            antennas.size()
        };
        std::cout << "Using " << n_data << " of " << max_possible_n_data << " possible data.\n";
        Eigen::VectorXcd Ez_tot_meas_concat(n_data);
        Eigen::VectorXcd Ez_bkg_calc_concat(n_data);

        std::cout << "\tConcatenating all data and all background fields into big vectors.\n";
        auto concat_ptr{0};
        for(auto ifreq{0}; ifreq<Ez_tot_meas.size();++ifreq){
            for(auto itx{0};itx<Ez_tot_meas[ifreq].size();++itx){
                const auto n_data_for_tx{M_s_data[ifreq][itx].rows()};
                Ez_tot_meas_concat.block(
                    concat_ptr,0,
                    n_data_for_tx,1
                ) = M_s_data[ifreq][itx]*Ez_tot_meas[ifreq][itx].getValRef();
                Ez_bkg_calc_concat.block(
                    concat_ptr,0,
                    n_data_for_tx,1
                ) = M_s_data[ifreq][itx]*Ez_bkg_dat[ifreq][itx].getValRef();
                concat_ptr+=n_data_for_tx;
            }
        }
        std::cout << "\tFilled in " << concat_ptr << " rows, and I should have filled " << n_data << " rows.\n";
        const Eigen::VectorXcd MsDeltaUTot{Ez_tot_meas_concat-Ez_bkg_calc_concat};

        std::cout << "\tBuilding D, from distorted Green's functions and current background field.\n ";
        Eigen::MatrixXcd D_bkg(n_data,mesh.centroids.rows());
        auto D_ptr{0};
        for(auto ifreq{0};ifreq<frequencies.size();++ifreq){
            Eigen::MatrixXcd kG{
                k2_bs[ifreq]*
                G_bkg_dat_by_freq[ifreq]
            };
            for(auto itx{0};itx<antennas.size();++itx){
                Eigen::MatrixXcd D_block{
                        M_s_data[ifreq][itx]*
                        kG*
                        Ez_bkg_dom[ifreq][itx].getValRef().asDiagonal()
                };
                D_bkg.block(
                    D_ptr,0,
                    D_block.rows(),D_block.cols()
                ) = D_block;
                D_ptr += D_block.rows();
            }
        }



        // Build a big vertcat of kb^2*G_bkg_dat*[Ubkgdom]

        // Solve for chi^tot_bkg via truncated svd.
        // u^{tot,meas} - u^{bkg}_{dat} =
        // k_b^2*G^dat_bkg*U^bkg_dom*(eps_rel_fgd-eps_rel_bkg)
        // std::cout << "\tDecomposing D.\n";
        // Eigen::BDCSVD<Eigen::MatrixXcd> SVD_D_bkg;
        // SVD_D_bkg.compute(
        //     D_bkg,
        //     Eigen::DecompositionOptions::ComputeThinU|
        //     Eigen::DecompositionOptions::ComputeThinV
        // );
        //
        // Eigen::VectorXd D_singular_values{SVD_D_bkg.singularValues()};
        // //std::cout << D_singular_values.transpose() << "\n";
        // auto n_sing_vals_to_keep{0};
        // for(
        //     n_sing_vals_to_keep=0;
        //     n_sing_vals_to_keep<D_singular_values.size();
        //     ++n_sing_vals_to_keep
        // ){
        //     if(D_singular_values(n_sing_vals_to_keep)<(D_singular_values(0)*sing_val_thresh)){
        //         break;
        //     }
        // }
        // std::cout << "\tTruncating U,S,V. Use " << n_sing_vals_to_keep << " of " << D_singular_values.size() << " singular vectors.\n";
        // Eigen::MatrixXcd U_trunc{SVD_D_bkg.matrixU()};
        // U_trunc = U_trunc.block(
        //     0,0,
        //     U_trunc.rows(),n_sing_vals_to_keep
        // );
        //
        // Eigen::MatrixXcd V_trunc{SVD_D_bkg.matrixV()};
        // V_trunc = V_trunc.block(
        //     0,0,
        //     V_trunc.rows(),n_sing_vals_to_keep
        // );
        //
        // const Eigen::VectorXcd S_trunc_vec{
        //     D_singular_values.block(0,0,n_sing_vals_to_keep,1)
        // };
        // const Eigen::MatrixXcd S_trunc_mat_inv{
        //     S_trunc_vec.array().inverse().matrix().asDiagonal()
        // };
        // std::cout << "\tSolving for delta_eps_rel_tot_bkg\n";
        // Eigen::VectorXcd delta_eps_rel_tot_bkg{
        //     contrast_update_scale*
        //     V_trunc*(
        //         S_trunc_mat_inv*(
        //             U_trunc.adjoint()*MsDeltaUTot
        //         )
        //     )
        // };

        Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXcd> lscg;
        lscg.compute(D_bkg);
        lscg.setMaxIterations(n_lscg_iterations);
        lscg.setTolerance(lscg_tolerance);

        const Eigen::VectorXcd delta_eps_rel_tot_bkg{lscg.solve(MsDeltaUTot)};

        // Assign T/B contrast to B/I contrast
        std::cout << "\tUpdating contrast.\n";
        const Eigen::VectorXcd new_contrast{
            target_tot.contrast.getValRef()+
            delta_eps_rel_tot_bkg
        };
        target_tot.contrast.setVals(new_contrast);
        DBIMStep step;
        step.chi = target_tot.contrast;
        step.Utot = Ez_bkg_dom;
        step.Fs = MsDeltaUTot.norm()/Ez_tot_meas_concat.norm();
        inv_log.steps.push_back(step);

        std::cout << "\tStoring DBIM step.\n";
        std::cout << "\tavg. abs. eps. step : " << delta_eps_rel_tot_bkg.array().abs().mean() <<"\n";
        std::cout << "\tCost functional    : " << step.Fs << "\n";
        std::cout << "\tData norm          : " << Ez_tot_meas_concat.norm() << "\n";
        std::cout << "\tBkg field norm     : " << Ez_bkg_calc_concat.norm() << "\n";
    }
    return(inv_log);
}

BIMInversion Chamber::bornIterativeMethod(){

    //For each frequency, make sure I have my G_inc_dom and G_inc_data.
    assert(frequencies.size());
    assert(probes.size());
    assert(antennas.size());
    assert(tx2rx.size() == frequencies.size());

    assert(Ez_tot_meas.size());

    calcDataEzInc();
    Ez_sct_meas.resize(Ez_tot_meas.size());
    {
        auto vf_idx{0};
        for(auto& vf:Ez_tot_meas){
            auto ff_idx{0};
            Ez_sct_meas[vf_idx].resize(vf.size());
            for(auto& ff:vf){
                // std::cout << ff.getValRef().transpose() << "\n";
                // std::cout << Ez_inc_d[vf_idx][ff_idx].getValRef().transpose() << "\n";
                Ez_sct_meas[vf_idx][ff_idx].setVals(
                    ff.getValRef()-Ez_inc_d[vf_idx][ff_idx].getValRef()
                );
                ff_idx++;
            }
            vf_idx++;
        }
    }


    G_inc_domain_by_freq.resize(frequencies.size());

    for(auto ifreq{0};ifreq<frequencies.size();ifreq++){
        mesh.buildDataGreen(
            G_inc_data_by_freq[ifreq],
            k2_bs[ifreq],
            probe_points
        );
        mesh.buildDomainGreen(
            G_inc_domain_by_freq[ifreq],
            k2_bs[ifreq]
        );
    }

    // Build all my Ms matrices
    M_s_data.resize(frequencies.size());
    for(auto& M_s_vec : M_s_data){
        M_s_vec.resize(antennas.size());
        for(auto& M_s : M_s_vec){
            M_s.resize(
                probes.size(),
                probes.size()
            );
        }
    }

    for(auto ifreq{0}; ifreq<tx2rx.size(); ++ifreq){
        for(auto itx{0}; itx < tx2rx[ifreq].size(); ++itx){
            auto row_ptr{0};
            M_s_data[ifreq][itx].resize(
                tx2rx[ifreq][itx].size(),
                probes.size()
            );
            M_s_data[ifreq][itx].setZero();
            for(auto& rx_idx : tx2rx[ifreq][itx]){
                M_s_data[ifreq][itx](row_ptr++,rx_idx)=1.0;
            }
        }
    }

    double sing_val_thresh{-1};
    while(not(0.0<sing_val_thresh && sing_val_thresh<1.0)){
        std::cout << "What threshold shall I use for singular values?[0-1]\n";
        std::cin >> sing_val_thresh;
    }
    int n_bim_iterations{-1};
    while(n_bim_iterations<0){
        std::cout << "How many BIM iterations shall we run?\n";
        std::cin >> n_bim_iterations;
    }

    BIMInversion inv_log;
    inv_log.imaging_mesh = mesh;
    inv_log.Ez_sct_meas = Ez_sct_meas;
    for(auto bim_iter{0};bim_iter<n_bim_iterations;++bim_iter){
        calcDomainEzTot();
        const auto n_cells{mesh.areas.size()};
        auto total_data_count{0};
        for(auto& tx_rx_map_at_freq : tx2rx){
            for(auto& rx_list : tx_rx_map_at_freq){
                total_data_count += rx_list.size();
            }
        }
        std::cout << "Starting DU product\n";
        auto du_row_ptr{0};
        Eigen::MatrixXcd DU(total_data_count,n_cells);
        for(auto ifreq{0};ifreq<Ez_tot.size();++ifreq){
            for(auto itx{0};itx<Ez_tot[ifreq].size();++itx){
                const auto n_data_for_tx{M_s_data[ifreq][itx].rows()};
                DU.block(
                    du_row_ptr,0,
                    n_data_for_tx,n_cells
                ) =
                M_s_data[ifreq][itx]*
                k2_bs[ifreq]*
                G_inc_data_by_freq[ifreq]*
                Ez_tot[ifreq][itx].getValRef().asDiagonal();
                du_row_ptr += n_data_for_tx;
            }
        }
        std::cout << "Done DU product\n";

        // Now get an svd from DU

        std::cout << "Beginning SV decomposition...\n";

        Eigen::BDCSVD<Eigen::MatrixXcd> SVD_DU;
        SVD_DU.compute(
            DU,
            Eigen::DecompositionOptions::ComputeThinU|
            Eigen::DecompositionOptions::ComputeThinV
        );

        std::cout << "Decomposition complete!\n";

        std::cout << "Truncating SVD matrices...\n";
        Eigen::VectorXd DU_singular_values{SVD_DU.singularValues()};

        auto n_sing_vals_to_keep{0};
        for(
            n_sing_vals_to_keep=0;
            n_sing_vals_to_keep<DU_singular_values.size();
            ++n_sing_vals_to_keep
        ){
            if(DU_singular_values(n_sing_vals_to_keep)<(DU_singular_values(0)*sing_val_thresh)){
                break;
            }
        }
        std::cout << "Gonna use " << n_sing_vals_to_keep << " singular values, out of " << DU_singular_values.size() << "\n";
        Eigen::MatrixXcd U_trunc{SVD_DU.matrixU()};
        U_trunc = U_trunc.block(
            0,0,
            U_trunc.rows(),n_sing_vals_to_keep
        );

        Eigen::MatrixXcd V_trunc{SVD_DU.matrixV()};
        V_trunc = V_trunc.block(
            0,0,
            V_trunc.rows(),n_sing_vals_to_keep
        );

        Eigen::VectorXcd S_trunc_vec{DU_singular_values.block(0,0,n_sing_vals_to_keep,1)};
        std::cout << "Done truncation\n";
        std::cout << U_trunc.rows() << "," << U_trunc.cols() << "," << S_trunc_vec.rows() << "," << V_trunc.cols() << "," << V_trunc.rows() << "\n";

        // Construct the big concatenation of measured data, masked with M_s_data.
        std::cout << "Constructing concatenated data vector\n";
        Eigen::VectorXcd Ez_sct_meas_masked(total_data_count);
        auto Ez_meas_mask_ptr{0};
        for(auto ifreq{0};ifreq<frequencies.size();++ifreq){
            for(auto itx{0};itx<antennas.size();++itx){
                const auto n_data_for_tx{M_s_data[ifreq][itx].rows()};
                Ez_sct_meas_masked.block(
                    Ez_meas_mask_ptr,0,
                    n_data_for_tx,1
                ) = M_s_data[ifreq][itx]*Ez_sct_meas[ifreq][itx].getValRef();
                Ez_meas_mask_ptr += n_data_for_tx;
            }
        }

        std::cout << "Done constructing concatenated data vector\n";

        // Use truncated SVD to solve for chi
        std::cout << "Solving for chi, with truncated svd...\n";
        const Eigen::MatrixXcd S_trunc_mat_inv{
            S_trunc_vec.array().inverse().matrix().asDiagonal()
        };
        Eigen::VectorXcd chi_eps{
            V_trunc*(
                S_trunc_mat_inv*(
                    U_trunc.adjoint()*Ez_sct_meas_masked
                )
            )
        };
        std::cout << "Done chi calculation.\n";
        // Now set the target, given chi:
        Eigen::VectorXcd eps_rel{(chi_eps.array()+1).matrix()};

        target_tot.eps_r.setVals(eps_rel);


        Field chi_for_step;
        chi_for_step.setLocations(inv_log.imaging_mesh.centroids);
        chi_for_step.setVals(chi_eps);
        BIMStep step;
        step.chi = chi_for_step;
        step.Utot = Ez_tot;
        step.Fs = (Ez_sct_meas_masked - DU*chi_eps).norm()/
            (Ez_sct_meas_masked.norm());
        inv_log.steps.push_back(step);
        std::cout << "Error at iteration " << bim_iter << " is "<<  step.Fs << "\n";
    }
    return(inv_log);
}

void DBIMInversion::WriteResults(std::string outdir){
    if(not(outdir.back()=='\\' or outdir.back() == '/')) outdir += "/";
    imaging_mesh.WriteMeshToFile(outdir);
    auto iter_count{0};
    for(auto& step : steps){
        step.chi.WriteValsToFile(
            outdir +
            "chi_iter_" + std::to_string(iter_count) +
            ".txt"
        );
        auto step_freq_count{0};
        for(auto& vec_of_tot_fields : step.Utot){
            auto field_count{0};
            for(auto& tot_field:vec_of_tot_fields){
                tot_field.WriteValsToFile(
                    outdir +
                    "Ez_tot_" +
                    "iter_" + std::to_string(iter_count) +
                    "freq_" + std::to_string(step_freq_count) +
                    "tx_" + std::to_string(field_count) +
                    ".txt"
                );
                field_count++;
            }
            step_freq_count++;
        }
        iter_count++;
    }
    auto freq_count{0};
    for(auto& vec_of_tot_data:Ez_tot_meas){
        Eigen::MatrixXcd meas_mat;
        FormFieldMatrix(
            meas_mat,
            vec_of_tot_data
        );
        WriteMatrixToFile(
            outdir +
            "Ez_tot_meas_" +
            "freq_" + std::to_string(freq_count) +
            ".txt",
            meas_mat
        );
        freq_count++;
    }
    Eigen::VectorXcd Fs_vec(steps.size());
    for(auto ii{0}; ii<steps.size(); ++ii){
        Fs_vec(ii) = steps[ii].Fs;
    }
    WriteVectorToFile(
        outdir + "Fs.txt",
        Fs_vec
    );
}

void BIMInversion::WriteResults(std::string outdir){
    if(not(outdir.back()=='\\' or outdir.back() == '/')) outdir += "/";
    imaging_mesh.WriteMeshToFile(outdir);
    auto iter_count{0};
    for(auto& step : steps){
        step.chi.WriteValsToFile(
            outdir +
            "chi_iter_" + std::to_string(iter_count) +
            ".txt"
        );
        auto step_freq_count{0};
        for(auto& vec_of_tot_fields : step.Utot){
            auto field_count{0};
            for(auto& tot_field:vec_of_tot_fields){
                tot_field.WriteValsToFile(
                    outdir +
                    "Ez_tot_" +
                    "iter_" + std::to_string(iter_count) +
                    "freq_" + std::to_string(step_freq_count) +
                    "tx_" + std::to_string(field_count) +
                    ".txt"
                );
                field_count++;
            }
            step_freq_count++;
        }
        iter_count++;
    }
    auto freq_count{0};
    for(auto& vec_of_tot_data:Ez_sct_meas){
        Eigen::MatrixXcd meas_mat;
        FormFieldMatrix(
            meas_mat,
            vec_of_tot_data
        );
        WriteMatrixToFile(
            outdir +
            "Ez_sct_meas_" +
            "freq_" + std::to_string(freq_count) +
            ".txt",
            meas_mat
        );
        freq_count++;
    }
    Eigen::VectorXcd Fs_vec(steps.size());
    for(auto ii{0}; ii<steps.size(); ++ii){
        Fs_vec(ii) = steps[ii].Fs;
    }
    WriteVectorToFile(
        outdir + "Fs.txt",
        Fs_vec
    );
}

void WriteMatrixToFile(
    std::string filename,
    const Eigen::MatrixXcd & matrix,
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
    writer << matrix.rows() << "\t" << matrix.cols() << std::endl;
    for(int ii = 0; ii < matrix.rows(); ii++)
    {
        int n_cols= matrix.cols();
        for(int jj = 0; jj <n_cols; jj++)
        {
            std::string sep = "\t";
            if (jj == (n_cols-1)) sep = "";
            writer << matrix(ii,jj) << sep;
        }
        writer << std::endl;
    }
    writer.close();
}

void WriteMatrixToFile(
    std::string filename,
    const Eigen::MatrixXd & matrix,
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
    writer << matrix.rows() << "\t" << matrix.cols() << std::endl;
    for(int ii = 0; ii < matrix.rows(); ii++)
    {
        int n_cols= matrix.cols();
        for(int jj = 0; jj <n_cols; jj++)
        {
            std::string sep = "\t";
            if (jj == (n_cols-1)) sep = "";
            writer << matrix(ii,jj) << sep;
        }
        writer << std::endl;
    }
    writer.close();
}

void WriteVectorToFile(
    std::string filename,
    const Eigen::VectorXd & vec,
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
    writer << vec.size() << std::endl;
    for(int ii = 0; ii < vec.size(); ii++)
    {
        writer << vec(ii) << std::endl;
    }
    writer.close();
}

void WriteVectorToFile(
    std::string filename,
    const Eigen::VectorXcd & vec,
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
    writer << vec.size() << std::endl;
    for(int ii = 0; ii < vec.size(); ii++)
    {
        writer << vec(ii) << std::endl;
    }
    writer.close();
}

void ReadVectorFromFile(
    std::string filename,
    Eigen::VectorXd & vec
)
{
    std::ifstream reader;
    reader.open(filename, std::ifstream::in);
    int nrows;
    reader >> nrows;
    vec.resize(nrows);
    for(int ii = 0; ii < nrows; ii++)
    {
        reader >> vec(ii);
    }
    reader.close();
}


void ReadVectorFromFile(
    std::string filename,
    Eigen::VectorXcd & vec
)
{
    std::ifstream reader;
    reader.open(filename, std::ifstream::in);
    int nrows;
    reader >> nrows;
    vec.resize(nrows);
    for(int ii = 0; ii < nrows; ii++)
    {
        reader >> vec(ii);
    }
    reader.close();
}

void ReadMatrixFromFile(
    std::string filename,
    Eigen::MatrixXd & matrix
)
{
    std::ifstream reader;
    reader.open(filename,std::ifstream::in);
    int nrows,ncols;
    reader >> nrows;
    reader >> ncols;
    matrix.resize(nrows,ncols);
    for(int mm = 0; mm < nrows;++mm)
    {
        for(int nn = 0; nn < ncols; ++nn)
        {
            reader >> matrix(mm,nn);
        }
    }
    reader.close();
}


void ReadMatrixFromFile(
    std::string filename,
    Eigen::MatrixXcd & matrix
)
{
    std::ifstream reader;
    reader.open(filename,std::ifstream::in);
    int nrows,ncols;
    reader >> nrows;
    reader >> ncols;

    matrix.resize(nrows,ncols);
    for(int mm = 0; mm < nrows;++mm)
    {
        for(int nn = 0; nn < ncols; ++nn)
        {
            reader >> matrix(mm,nn);
        }
    }
    reader.close();
}

void FormFieldMatrix(
    Eigen::MatrixXcd & M,
    std::vector<Field> & V
)
{
    assert(V.size());
    Eigen::VectorXcd vtest;
    V[0].getVals(vtest);
    auto nvals = vtest.size();
    M.resize(nvals,V.size());
    for(auto ii{0}; ii<V.size();++ii)
    {
        Eigen::VectorXcd vi;
        V[ii].getVals(vi);
        assert(vi.size()==nvals);
        M.col(ii) = vi;
    }
}
