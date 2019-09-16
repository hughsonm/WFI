#include <eigen3/Eigen/Dense>
#include <gmsh.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

int main (int argc, char **argv)
{
    assert(argc==2);

    std::string bar = "";
    bar += "-";
    for(int bar_idx = 0; bar_idx < 80; bar_idx++) bar += "-";

    gmsh::initialize();
    gmsh::open(argv[1]);

    std::vector<std::string> model_names;
    gmsh::model::list(model_names);
    std::cout << "Listing model names: " << std::endl;
    for(int ii = 0; ii < model_names.size(); ii++)
    {
        std::cout << model_names[ii] << std::endl;
    }

    // Entities
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags);
    std::cout << "Entities" << std::endl;
    std::cout << "This is basically me reading you the" << std::endl;
    std::cout << "geo file" << std::endl;
    std::cout << bar << std::endl;
    for(int ent_idx = 0; ent_idx < dimTags.size(); ent_idx++)
    {

        int ent_dim = dimTags[ent_idx].first;
        int ent_tag = dimTags[ent_idx].second;
        std::cout << "Entity Number: " << ent_idx << std::endl;
        std::cout << "\tDimension: " << ent_dim << std::endl;
        std::cout << "\tTag: " << ent_tag << std::endl;
        std::string ent_name;

        gmsh::model::getEntityName(
            ent_dim,
            ent_tag,
            ent_name
        );

        std::cout << "\tName: " << ent_name << std::endl;

        std::vector<int> phys_tags;
        gmsh::model::getPhysicalGroupsForEntity(
            ent_dim,
            ent_tag,
            phys_tags
        );
        std::cout << "\tPhysical Groups: ";
        for(int phys_idx = 0; phys_idx < phys_tags.size(); phys_idx++)
        {
            std::cout << phys_tags[phys_idx] << "\t";
        }
        std::cout << std::endl;
        std::string ent_type;
        gmsh::model::getType(
            ent_dim,
            ent_tag,
            ent_type
        );
        std::cout << "\tType: " << ent_type << std::endl;

        int parent_dim,parent_tag;
        gmsh::model::getParent(
            ent_dim,
            ent_tag,
            parent_dim,
            parent_tag
        );
        std::cout << "\tParent Info:" << std::endl;
        if(parent_dim==-1 || parent_tag==-1)
        {
            std::cout << "\t\tEntity has no parent" << std::endl;
        }
        else
        {
            std::cout << "\t\tDimension: " << parent_dim << std::endl;
            std::cout << "\t\tTag: " << parent_tag << std::endl;
        }


        std::cout << std::endl;
    }

    gmsh::model::getPhysicalGroups(dimTags);
    std::cout << "Physical Groups" << std::endl;
    std::cout << "This is basically me reading you the" << std::endl;
    std::cout << "important mesh information." << std::endl;
    std::cout << bar << std::endl;
    for(int phys_idx = 0; phys_idx < dimTags.size(); phys_idx++)
    {
        std::cout << "Physical Group Number: " << phys_idx;
        std::cout << std::endl;
        int phys_dim = dimTags[phys_idx].first;
        int phys_tag = dimTags[phys_idx].second;
        std::cout << "\tDimension: " << phys_dim << std::endl;
        std::cout << "\tTag: " << phys_tag << std::endl;;
        std::string phys_name;
        gmsh::model::getPhysicalName(
            phys_dim,
            phys_tag,
            phys_name
        );

        std::cout << "\tGetting entities for this physical group's";
        std::cout << " dimension" << std::endl;
        gmsh::vectorpair entities_dimtags_for_phys;
        gmsh::model::getEntities(
            entities_dimtags_for_phys,
            phys_dim
        );
        std::cout << "\tPhysical Group Name: ";
        std::cout << phys_name << std::endl;
        std::vector<int> ent_tags;
        gmsh::model::getEntitiesForPhysicalGroup(
            phys_dim,
            phys_tag,
            ent_tags
        );
        std::cout << "\tEntity Count: " << ent_tags.size();
        std::cout << std::endl;

        //std::vector<int> ele_tags_for_phys;
        std::vector<std::vector<int> > node_tags_for_ele_for_phys;

        for(int ent_idx = 0; ent_idx < ent_tags.size(); ent_idx++)
        {
            int ent_for_phys_dim = entities_dimtags_for_phys[ent_idx].first;
            int ent_for_phys_tag = entities_dimtags_for_phys[ent_idx].second;

            std::cout << "\t\tThis physical group has the entity";
            std::cout << " with dim = ";
            std::cout << ent_for_phys_dim;
            std::cout << " and tag = ";
            std::cout << ent_for_phys_tag;
            std::cout << std::endl;

            std::string ent_for_phys_type;
            gmsh::model::getType(
                ent_for_phys_dim,
                ent_for_phys_tag,
                ent_for_phys_type
            );

            std::cout << "\t\tGetting elements and nodes for this";
            std::cout << " entity" << std::endl;

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
            std::cout << "\t\tGot " << ele_for_ent_types.size();
            std::cout << " different types of elements" << std::endl;
            for (int type_idx = 0; type_idx < ele_for_ent_types.size(); type_idx++)
            {
                std::cout << "\t\tType no. " << type_idx;
                std::cout << " is " << ele_for_ent_types[type_idx] << std::endl;

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
                std::cout << "\t\tType Name: " << type_name << std::endl;
                std::cout << "\t\tDimension: " << type_dim << std::endl;
                std::cout << "\t\tOrder    : " << type_order << std::endl;
                std::cout << "\t\tNum Nodes: " << type_num_nodes << std::endl;
                std::cout << "\t\tCount    : " << ele_for_ent_tags[type_idx].size() << std::endl;
                std::cout << "\t\tTags     : " << std::endl;
                for(int ele_idx = 0; ele_idx < ele_for_ent_tags[type_idx].size(); ele_idx++)
                {
                    //std::cout << "\t\t" << ele_for_ent_tags[type_idx][ele_idx] << std::endl;
                    //ele_tags_for_phys.push_back(ele_for_ent_tags[type_idx][ele_idx]);
                    std::vector<int> this_eles_node_tags(type_num_nodes);
                    for(int cc = 0; cc < type_num_nodes; cc++)
                    {
                        this_eles_node_tags[cc] = node_for_ent_tags[type_idx][node_tag_ptr++];
                    }
                    node_tags_for_ele_for_phys.push_back(this_eles_node_tags);
                }
            }
        }

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



    // // Get elements and nodes so we can see where stuff sits.
    // std::vector<std::size_t> node_tags;
    // std::vector<double> node_coords_flat;
    // std::vector<double> node_para_coords_flat;
    // std::vector<std::vector<double> > node_coords;
    // gmsh::model::mesh::getNodes(
    //     node_tags,
    //     node_coords_flat,
    //     node_para_coords_flat
    // );
    //
    // // Want to sort the nodes into a vector such that I can access
    // // nodes using their tags and their indices
    // std::size_t max_tag = *(max_element(node_tags.begin(),node_tags.end()));
    // std::cout << "Nodes" << std::endl;
    // std::cout << "\tMax. Tag: " << max_tag << std::endl;
    //
    // node_coords.resize(max_tag+1);
    // std::cout << "\tTag\tX\tY\tZ" << std::endl;
    // for(int tag_idx = 0; tag_idx < node_tags.size(); tag_idx++)
    // {
    //     std::vector<double> node(3);
    //     node[0] = node_coords_flat[3*tag_idx+0];
    //     node[1] = node_coords_flat[3*tag_idx+1];
    //     node[2] = node_coords_flat[3*tag_idx+2];
    //     node_coords[node_tags[tag_idx]] = node;
    //     std::cout << "\t" << node_tags[tag_idx];
    //     std::cout << "\t" << node[0];
    //     std::cout << "\t" << node[1];
    //     std::cout << "\t" << node[2];
    //     std::cout << std::endl;
    // }
    //
    //
    //
    //
    // std::vector<int> ele_types;
    // std::vector<std::vector<std::size_t> > ele_tags;
    // std::vector<std::vector<std::size_t> > ele_node_tags;
    // gmsh::model::mesh::getElements(
    //     ele_types,
    //     ele_tags,
    //     ele_node_tags,
    //     2
    // );
    // std::cout << "Elements" << std::endl;
    // std::cout << "Now that we know the physical groups," << std::endl;
    // std::cout << "let's build some elements and" << std::endl;
    // std::cout << "calculate some stuff" <<  std::endl;



    // A matrix
    Eigen::MatrixXd m(2,2);
    m(0,0) = 1;
    m(0,1) = 2;
    m(1,0) = 3;
    m(1,1) = 4;
    std::cout << "My matrix is:\n" << m << std::endl;
    gmsh::finalize();
    return(0);
}
