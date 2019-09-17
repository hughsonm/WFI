#include <eigen3/Eigen/Dense>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

void BuildMesh(
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
            std::cout << "Physical Group Number: " << phys_idx;
            std::cout << std::endl;
        }
        int phys_dim = dimTags[phys_idx].first;
        int phys_tag = dimTags[phys_idx].second;
        if(phys_dim != 2)
        {
            continue;
        }
        if(verbose)
        {
            std::cout << "\tDimension: " << phys_dim << std::endl;
            std::cout << "\tTag: " << phys_tag << std::endl;
        }
        std::string phys_name;
        gmsh::model::getPhysicalName(
            phys_dim,
            phys_tag,
            phys_name
        );

        if(verbose)
        {
            std::cout << "\tGetting entities for this physical group's";
            std::cout << " dimension" << std::endl;
        }
        gmsh::vectorpair entities_dimtags_for_phys;
        gmsh::model::getEntities(
            entities_dimtags_for_phys,
            phys_dim
        );
        if(verbose)
        {
            std::cout << "\tPhysical Group Name: ";
            std::cout << phys_name << std::endl;
        }
        std::vector<int> ent_tags;
        gmsh::model::getEntitiesForPhysicalGroup(
            phys_dim,
            phys_tag,
            ent_tags
        );
        if(verbose)
        {
            std::cout << "\tEntity Count: " << ent_tags.size();
            std::cout << std::endl;
        }

        std::vector<std::vector<int> > node_tags_for_ele_for_phys;

        for(int ent_idx = 0; ent_idx < ent_tags.size(); ent_idx++)
        {
            int ent_for_phys_dim = entities_dimtags_for_phys[ent_idx].first;
            int ent_for_phys_tag = entities_dimtags_for_phys[ent_idx].second;

            if(verbose)
            {
                std::cout << "\t\tThis physical group has the entity";
                std::cout << " with dim = ";
                std::cout << ent_for_phys_dim;
                std::cout << " and tag = ";
                std::cout << ent_for_phys_tag;
                std::cout << std::endl;
            }

            std::string ent_for_phys_type;
            gmsh::model::getType(
                ent_for_phys_dim,
                ent_for_phys_tag,
                ent_for_phys_type
            );

            if(verbose)
            {
                std::cout << "\t\tGetting elements and nodes for this";
                std::cout << " entity" << std::endl;
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
                std::cout << "\t\tGot " << ele_for_ent_types.size();
                std::cout << " different types of elements" << std::endl;
            }
            for (int type_idx = 0; type_idx < ele_for_ent_types.size(); type_idx++)
            {
                if(verbose)
                {
                    std::cout << "\t\tType no. " << type_idx;
                    std::cout << " is " << ele_for_ent_types[type_idx] << std::endl;
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
                    std::cout << "\t\tType Name: " << type_name << std::endl;
                    std::cout << "\t\tDimension: " << type_dim << std::endl;
                    std::cout << "\t\tOrder    : " << type_order << std::endl;
                    std::cout << "\t\tNum Nodes: " << type_num_nodes << std::endl;
                    std::cout << "\t\tCount    : " << ele_for_ent_tags[type_idx].size() << std::endl;
                    std::cout << "\t\tTags     : " << std::endl;
                }
                // How many elements we got? A lot.
                int n_new_elements = ele_for_ent_tags[type_idx].size();
                Eigen::Index cur_n_rows = tri.rows();
                if(verbose)
                {
                    std::cout << "\t\tTri has " << cur_n_rows << " rows";
                    std::cout << " and we add " << n_new_elements;
                    std::cout << " for a total of ";
                    std::cout << (cur_n_rows + n_new_elements);
                    std::cout << std::endl;
                }


                tri.conservativeResize(
                    cur_n_rows + n_new_elements,
                    Eigen::NoChange
                );
                for(int ele_idx = 0; ele_idx < n_new_elements; ele_idx++)
                {
                    //std::cout << "\t\t" << ele_for_ent_tags[type_idx][ele_idx] << std::endl;
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
            std::cout << "\tThis physical's node tags for each ele:";
            std::cout << std::endl;
            for(int ele_idx = 0; ele_idx < node_tags_for_ele_for_phys.size(); ele_idx++)
            {
                std::cout << "\t";
                for(int node_idx = 0; node_idx < node_tags_for_ele_for_phys[ele_idx].size(); node_idx++)
                {
                    std::cout << node_tags_for_ele_for_phys[ele_idx][node_idx] << "\t";
                }
                std::cout << std::endl;
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

int main (int argc, char **argv)
{
    assert(argc==2);

    std::string bar = "";
    bar += "-";
    for(int bar_idx = 0; bar_idx < 80; bar_idx++) bar += "-";

    gmsh::initialize();

    Eigen::MatrixXd tri;
    Eigen::MatrixXd points;

    BuildMesh(
        argv[1],
        tri,
        points
    );

    std::cout << bar << std::endl;
    std::cout << "Tri: " << std::endl;
    std::cout << tri << std::endl;

    std::cout << bar << std::endl;
    std::cout << "Points: " << std::endl;
    std::cout << points << std::endl;

    gmsh::finalize();
    return(0);
}
