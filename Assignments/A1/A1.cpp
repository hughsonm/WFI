#include <eigen3/Eigen/Dense>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>

int main (int argc, char **argv)
{
    gmsh::initialize();
    gmsh::open("TestMesh.msh");
    std::vector<std::string> model_names;
    gmsh::model::list(model_names);
    std::cout << "Listing model names: " << std::endl;
    for(int ii = 0; ii < model_names.size(); ii++)
    {
        std::cout << model_names[ii] << std::endl;
    }

    gmsh::vectorpair physical_groups;
    gmsh::model::getPhysicalGroups(physical_groups);
    std::cout << "These are my physical groups:" << std::endl;
    for (int ii = 0; ii < physical_groups.size(); ii++)
    {
        std::cout << "Dimension: " << physical_groups[ii].first << ", Tag: " << physical_groups[ii].second << std::endl;
    }
    std::vector<int> ele_types;
    std::vector<std::vector<std::size_t> > ele_tags;
    std::vector<std::vector<std::size_t> > node_tags;
    gmsh::model::mesh::getElements(ele_types,ele_tags,node_tags);
    for (int jj = 0; jj < ele_types.size(); jj++)
    {
        std::cout << "Type No. " << jj << " is " << ele_types[jj] << std::endl;
        std::string type_name;
        int type_dim;
        int type_ord;
        int type_n_nodes;
        std::vector<double> type_coords;
        gmsh::model::mesh::getElementProperties(ele_types[jj], type_name, type_dim, type_ord, type_n_nodes, type_coords);
        std::cout << "Name:\t" << type_name << "\n\tDimension:\t" << type_dim << "\n\tOrder:\t" << type_ord << "\n\tNodes:\t" << type_n_nodes << "\n\tReference Coordinates:"<< std::endl;
        for (int nn = 0; nn < type_n_nodes; nn++)
        {
            std::cout << "\t";
            for (int dd = 0; dd < type_dim; dd++)
            {
                std::cout << "\t" << type_coords[nn*type_dim+dd];
            }
            std::cout << std::endl;
        }
        for (int kk = 0; kk < ele_tags[jj].size(); kk++)
        {
            std::cout << "\t" << ele_tags[jj][kk] << std::endl;
        }
    }
    Eigen::MatrixXd m(2,2);
    m(0,0) = 1;
    m(0,1) = 2;
    m(1,0) = 3;
    m(1,1) = 4;
    std::cout << "My matrix is:\n" << m << std::endl;
    gmsh::finalize();
    return(0);
}
