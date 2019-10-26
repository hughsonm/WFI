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
                Ez(ipt) = -j_imag/4.0*boost::math::cyl_hankel_2(
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
    int n_tri = areas.size();
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
            // Background green's function is symmetric
            G(mm,nn) = Gmn;
        }
    }
}


Chamber::Chamber(std::string meshfile)
{
    mesh.buildTriangulation(meshfile);

    // Set the target to be zero contrast.

    target.eps_r.setLocations(mesh.centroids);

    target.eps_r.erase();
    target.ready = true;

    // target.eps_r.resize(mesh.areas.size());
    // for(int rr = 0; rr<eps_r.size();++rr) eps_r(rr) = 1.0;
    // k2_b = 0;
    // k2_f = eps_r*k2_b;
}

void Chamber::setTarget(std::string targetfile)
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
    // Eigen::VectorXcd eps_to_set;
    // //eps_r.resize(mesh.tri.rows());
    // //k2_f.resize(mesh.tri.rows());
    // for(int rr = 0; rr < eps_r.size(); rr++)
    // {
    //     eps_r(rr) = 1.0;
    // }
    target.eps_r.erase();
    Eigen::VectorXcd eps_to_set;
    target.eps_r.getVals(eps_to_set);
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
                //k2_f(rr) = k2_b*eps_r(rr);
            }
        }
    }
    target.eps_r.setVals(eps_to_set);
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
    // std::vector<Eigen::Vector3d> vec_of_probes;
    probes.resize(0);
    while(reader)
    {
        Probe pi;
        // Eigen::Vector3d probe_xyz;
        reader >> x;
        if(reader.eof()) break;
        reader >> y;
        reader >> z;
        // probe_xyz << x,y,z;
        pi.location << x,y,z;
        pi.ready = true;
        // vec_of_probes.push_back(probe_xyz);
        probes.push_back(pi);
    }
    probe_points.resize(probes.size(),3);
    for(int pp = 0; pp < probes.size(); pp++)
    {
        probe_points.row(pp) = probes[pp].location.transpose();
    }

    reader.close();
}

void Chamber::readMeasuredData(std::string datafile)
{
    Eigen::MatrixXcd tot_matrix;
    ReadMatrixFromFile(
        datafile,
        tot_matrix
    );
    Ez_tot_meas.resize(tot_matrix.cols());
    for(auto cc{0}; cc<tot_matrix.cols(); ++cc)
    {
        Eigen::VectorXcd column{tot_matrix.col(cc)};
        Ez_tot_meas[cc].setLocations(probe_points);
        Ez_tot_meas[cc].setVals(column);
    }
}

void Chamber::buildDataGreen(void)
{
    mesh.buildDataGreen(
        G_b_data,
        k2_b,
        probe_points
    );
}

void Chamber::A2Q3(void)
{
	// Annihilating the Annihilator

	// Basis is green's functions, centered at receivers

	// Build a matrix H s.t. Hij is the inner product of :
	// (g_i)* and g_j

	// For this, we need a mesh, and the antennas.

	Eigen::MatrixXcd H;
	H.resize(antennas.size(),antennas.size());
	
	Eigen::MatrixXcd gg(mesh.centroids.rows(),antennas.size());
	for(auto tt{0};tt<antennas.size();++tt){
		Eigen::VectorXcd gt;
		antennas[tt].getEz(
		gt,
		mesh.centroids,
		frequency);
		gt *= k2_b;
		gg.col(tt) = gt;
	}
	

	for(auto ee{0};ee<mesh.centroids.rows();++ee){
		double area_for_ele{mesh.areas(ee)};
		for(auto ii{0};ii<gg.cols();++ii){
			std::complex<double> gistar_on_ele{
				gg(ee,ii).adjoint()
			};
			for(auto jj{0};jj<gg.cols();++jj){
				std::complex<double> gj_on_ele{
					gg(ee,jj)
				};
				H(ii,jj) += area_for_ele*gistar_on_ele*gj_on_ele;
			}
		}
	}
	
	Eigen::PartialPivLU<Eigen::MatrixXcd> H_LU;
	H_LU.compute(H);
	Eigen::MatrixXcd alphas(mesh.centroids.rows(),antennas.size());
	for(auto tt{0};tt<antennas.size();++tt){
		Eigen::VectorXcd Ezs_t;
		Ez_sct_d.getVals(Ezs_t);
		alphas.col(tt) = H_LU.solve(Ezs_t);
	}
	
}

void Chamber::A2Q5(void)
{
    auto ntx{antennas.size()};
    auto ntri{mesh.areas.size()};

    // Ez_sct_meas.resize(Ez_tot_meas.rows(),Ez_tot_meas.cols());
    Ez_sct_meas.resize(ntx);

    if(!Ez_inc_d.size())
    {
        calcDataEzInc();
    }

    Eigen::MatrixXcd Ez_sct_meas_mat(G_b_data.rows(),ntx);
    for(auto tt{0}; tt<ntx; ++tt)
    {
        Eigen::VectorXcd tot_t;
        Ez_tot_meas[tt].getVals(tot_t);
        Eigen::VectorXcd inc_t;
        Ez_inc_d[tt].getVals(inc_t);
        Eigen::VectorXcd sct_t{tot_t-inc_t};
        Ez_sct_meas[tt].setVals(sct_t);
        Ez_sct_meas_mat.col(tt) = sct_t;
    }

    // Ez_sct_meas = Ez_tot_meas - Ez_inc_d;

    if(!G_b_data.rows() || !G_b_data.cols())
    {
        std::cerr << "Annihilator needs G_b_data. Building it now..." << std::endl;
        buildDataGreen();
    }

    Eigen::MatrixXcd PhP{G_b_data.adjoint()*G_b_data};
    Eigen::MatrixXcd I{
        Eigen::MatrixXcd::Identity(
            PhP.rows(),PhP.cols()
        )
    };

    std::cerr << "PhP is " << PhP.rows() << " by " << PhP.cols() << std::endl;

    std::vector<bool> optimum_found_for_tx(ntx,false);
    std::vector<std::vector<double> > tikh_distances_for_tx(ntx);
    Eigen::MatrixXcd Phd_all{G_b_data.adjoint()*Ez_sct_meas_mat};
    std::vector<double> optimal_lambdas(ntx);
    Eigen::MatrixXcd optimal_alphas(G_b_data.cols(),ntx);
    double lambda_exp{0.0};
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

        double lambda_big{std::pow(10,lambda_exp)};
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
                // Solve system
                Eigen::VectorXcd alpha_big{LU_big.solve(Phd_all.col(itx))};
                Eigen::VectorXcd alpha_sml{LU_sml.solve(Phd_all.col(itx))};

                double res_norm_big{(Ez_sct_meas_mat.col(itx)-G_b_data*alpha_big).norm()};
                double alpha_norm_big{alpha_big.norm()};

                double res_norm_sml{(Ez_sct_meas_mat.col(itx)-G_b_data*alpha_sml).norm()};
                double alpha_norm_sml{alpha_sml.norm()};

                double tikh_dist_big{
                    std::sqrt(
                        std::pow(
                            std::log(
                                res_norm_big
                            ),
                            2
                        ) + std::pow(
                            std::log(
                                alpha_norm_big
                            ),
                            2
                        )
                    )
                };

                double tikh_dist_sml{
                    std::sqrt(
                        std::pow(
                            std::log(
                                res_norm_sml
                            ),
                            2
                        ) + std::pow(
                            std::log(
                                alpha_norm_sml
                            ),
                            2
                        )
                    )
                };

                tikh_distances_for_tx[itx].insert(
                    tikh_distances_for_tx[itx].begin(),
                    tikh_dist_sml
                );

                tikh_distances_for_tx[itx].push_back(
                    tikh_dist_big
                );

                // Check tikh_distances_for_tx[itx] for a dip.
                std::size_t ndists{tikh_distances_for_tx[itx].size()};
                if(2 < ndists)
                {
                    for(int id = 0; id < (ndists-2); ++id)
                    {
                        double dp0{tikh_distances_for_tx[itx][id+0]};
                        double dp1{tikh_distances_for_tx[itx][id+1]};
                        double dp2{tikh_distances_for_tx[itx][id+2]};
                        if((dp1<dp0) && (dp1<dp2))
                        {
                            optimum_found_for_tx[itx] = true;
                            double lambda_opt{
                                lambda_sml*std::pow(
                                    10,
                                    id+1
                                )
                            };
                            std::cerr << "\t\tOptimal lambda for tx " << itx << " is " << lambda_opt << std::endl;
                            optimal_lambdas[itx] = lambda_opt;
                            std::cerr << "\t\tAssigned!" << std::endl;
                        }
                    }
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
                optimal_alphas.col(ii) = LU_opt.solve(Phd_all.col(ii));
                std::cerr << "\tSolved and assigned!" << std::endl;
                alpha_is_calculated[ii] = true;
            }
        }


    }
    WriteMatrixToFile(
        "Alphas.txt",
        optimal_alphas
    );
    // Now I have all of my alphas. These alpha_i is a vector of  the
    // basis coefficients of the source which satisfies the scattered
    // data generated by transmitter_i.

    // From alpha, I can get a u^T, by:
    // u^T(r) = u^I(r) + G_b_domain*alpha
    if(!Ez_inc.size())
    {
        calcDomainEzInc();
    }
    if(!G_b_domain.rows() || !G_b_domain.cols())
    {
        std::cerr << "Annihilator need G_b_domain. Building it now..." << std::endl;
        mesh.buildDomainGreen(
            G_b_domain,
            k2_b
        );
    }
    std::cerr << "Making the matrix of total fields in imaging domain..." << std::endl;
    // std::cerr << "Ez_inc is " << Ez_inc.rows() << " by " << Ez_inc.size() << std::endl;
    // std::cerr << "G_b_domain is " << G_b_domain.rows() << " by " << G_b_domain.cols() << std::endl;
    // std::cerr << "optimal_alphas is " << optimal_alphas.rows() << " by " << optimal_alphas.cols() << std::endl;
    Eigen::MatrixXcd Ez_tot_opt(ntri,ntx);

    for(auto cc{0}; cc<Ez_tot_opt.cols(); ++cc)
    {
        Eigen::VectorXcd inc_c;
        Ez_inc[cc].getVals(inc_c);
        Ez_tot_opt.col(cc) = inc_c + G_b_domain*optimal_alphas.col(cc);
    }
    std::cerr << "success!" << std::endl;

    WriteMatrixToFile(
        "Ez_tot_opt.txt",
        Ez_tot_opt
    );
    Eigen::VectorXcd chi_opt(ntri);
    chi_opt.setZero();

    for(auto ii{0}; ii < ntri; ++ii)
    {
        //std::cerr << "Calculating contrast number " << ii << " of " << ntri << std::endl;
        std::complex<double> numerator{0.0},denominator{0.0};
        //std::cerr << "num and den start as " << numerator << "," << denominator << std::endl;
        for(auto jj{0}; jj<ntx; ++jj)
        {
            std::complex<double> uij{Ez_tot_opt(ii,jj)};
            std::complex<double> uijstar{std::conj(uij)};
            std::complex<double> aij{optimal_alphas(ii,jj)};
            //std::cerr << "\tPulling out " << uij << " and " << uijstar << std::endl;
            //std::cerr << "\tAlso using alpha = " << aij << std::endl;
            //std::cerr << "\tGonna add " << uijstar*aij << " to " << numerator << std::endl;
            //std::cerr << "\tGonna add " << uijstar*uij << " to " << denominator << std::endl;
            numerator += uijstar*aij;
            denominator += uijstar*uij;
            //std::cerr << "\tIt worked!!" << std::endl;
        }
        //std::cerr << "Now, num = " << numerator << ", and den = " << denominator << std::endl;
        //std::cerr << "Divide them, and we get " << numerator/denominator << std::endl;
        chi_opt(ii) = numerator/denominator;

    }
    WriteVectorToFile(
        "Chi.txt",
        chi_opt
    );

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

void Chamber::setFrequency(double freq)
{
    frequency = freq;
    double omega = 2*M_PI*freq;
    k2_b = omega*omega/CNAUGHT/CNAUGHT;
}

void Chamber::calcDomainEzInc(void)
{
    Ez_inc.resize(antennas.size());
    for(auto aa{0}; aa<antennas.size(); ++aa)
    {
        Field field_from_antenna;
        Eigen::VectorXcd ez_from_antenna;
        antennas[aa].getEz(
            ez_from_antenna,
            mesh.centroids,
            frequency
        );
        field_from_antenna.setLocations(mesh.centroids);
        field_from_antenna.setVals(ez_from_antenna);
        Ez_inc[aa] = field_from_antenna;
    }
}

void Chamber::calcDataEzInc(void)
{
    Ez_inc_d.resize(antennas.size());
    for(auto aa{0}; aa<antennas.size(); ++aa)
    {
        Field field_from_antenna_at_probes;
        Eigen::VectorXcd ez_from_antenna_at_probes;
        antennas[aa].getEz(
            ez_from_antenna_at_probes,
            probe_points,
            frequency
        );
        field_from_antenna_at_probes.setLocations(probe_points);
        field_from_antenna_at_probes.setVals(ez_from_antenna_at_probes);
        Ez_inc_d[aa] = field_from_antenna_at_probes;
    }
}

void Chamber::calcDomainEzTot(void)
{
    auto entry_point{0};
    switch(entry_point)
    {
        case 0:
        {
            // Calculate incident fields.
            //Ez_inc.resize(mesh.areas.size(),antennas.size());
            calcDomainEzInc();

            // Build domain green matrix
            std::cerr << "Building domain green..." << std::endl;
            mesh.buildDomainGreen(
                G_b_domain,
                k2_b
            );
            std::cerr << "G:(" << G_b_domain.rows() << ",";
            std::cerr << G_b_domain.cols() << ")" << std::endl;
            std::cerr << "Building Chi..." << std::endl;

            // Build contrast matrix
            Eigen::VectorXcd eps_r_from_target;
            target.eps_r.getVals(eps_r_from_target);
            Chi.resize(eps_r_from_target.size(),eps_r_from_target.size());
            Chi.setZero();
            for(int dd = 0; dd < Chi.cols();dd++)
            {
                Chi(dd,dd) = k2_b*eps_r_from_target(dd)-k2_b;
            }

            // Build domain L operator
            std::cerr << "Building L_domain" << std::endl;
            // L_domain.resize(G_b_domain.rows(),G_b_domain.cols());
            L_domain = -G_b_domain;

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
            Ez_sct.resize(antennas.size());
            for(auto ss{0}; ss < Ez_sct.size(); ss++)
            {
                Ez_sct[ss].setLocations(mesh.centroids);
            }
            // Calculate all the scattered fields
            for(int jj = 0; jj < Ez_inc.size(); jj++)
            {
                Eigen::VectorXcd inc_vec,sct_vec;
                Ez_inc[jj].getVals(inc_vec);
                sct_vec = LU_L.solve(
                    G_b_domain*(
                        Chi*inc_vec
                    )
                );
                Ez_sct[jj].setVals(sct_vec);
            }
        }
        case 2:
        {
            // Add scattered to indident to get total fields
            Ez_tot.resize(antennas.size());
            for(auto tt{0}; tt<antennas.size(); ++tt)
            {
                Eigen::VectorXcd inc_vec,sct_vec,tot_vec;
                Ez_inc[tt].getVals(inc_vec);
                Ez_sct[tt].getVals(sct_vec);
                tot_vec = inc_vec+sct_vec;
                Ez_tot[tt].setVals(tot_vec);
            }
        }
    }
}

void Chamber::calcDataEzTot(void)
{
    if(!(Ez_tot.size() == antennas.size()))
    {
        std::cerr << "Ez_tot was not built. Getting it now.";
        std::cerr << std::endl;
        calcDomainEzTot();
    }
    if(!G_b_data.rows() || !G_b_data.cols())
    {
        std::cerr << "Data green was not built. Getting it now";
        std::cerr << std::endl;
        mesh.buildDataGreen(
            G_b_data,
            k2_b,
            probe_points
        );
        // G_b_data_ready = true;
    }
    std::cerr << "Resizing d-field holders to ";
    std::cerr << probe_points.rows() << " by ";
    std::cerr << antennas.size() << std::endl;

    // Ez_inc_d.resize(probe_points.rows(),antennas.size());
    calcDataEzInc();
    Ez_sct_d.resize(antennas.size());
    Ez_tot_d.resize(antennas.size());

    // for(int aa = 0; aa < antennas.size(); aa++)
    // {
    //     std::cerr << "Calculating d-inc for antenna ";
    //     std::cerr << aa << std::endl;
    //     Eigen::VectorXcd Ez_inc_a;
    //     antennas[aa].getField(
    //         Ez_inc_a,
    //         probe_points
    //     );
    //     Ez_inc_d.col(aa) = Ez_inc_a;
    // }
    for(int tt = 0; tt < Ez_tot.size(); tt++)
    {
        std::cerr << "Calculating d-sct for antenna ";
        std::cerr << tt << std::endl;
        Eigen::VectorXcd Ez_tot_t;
        Ez_tot[tt].getVals(Ez_tot_t);
        Eigen::VectorXcd Ez_sct_d_t{G_b_data*(Chi*(Ez_tot_t))};
        Ez_sct_d[tt].setVals(Ez_sct_d_t);
        Eigen::VectorXcd Ez_inc_d_t;
        Ez_inc_d[tt].getVals(Ez_inc_d_t);
        Eigen::VectorXcd Ez_tot_d_t{Ez_sct_d_t+Ez_inc_d_t};
        Ez_tot_d[tt].setVals(Ez_tot_d_t);
    }
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
            Eigen::VectorXd dxyz;
            dxyz = (rxlocations.row(rr)-centroids.row(ee)).transpose();
            dxyz = dxyz.array().pow(2);
            double d_re = std::sqrt(dxyz.array().sum());
            double a_e = std::sqrt(areas(ee)/M_PI);

            std::complex<double> J1 = boost::math::cyl_bessel_j(
                1,
                k_b*a_e
            );
            std::complex<double> H02 = boost::math::cyl_hankel_2(
                0,
                k_b*d_re
            );
            std::complex<double> G_re = -j*M_PI*a_e*J1*H02/2.0/k_b;
            G(rr,ee) = G_re;
        }
    }
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
