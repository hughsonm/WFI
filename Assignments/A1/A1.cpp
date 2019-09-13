#include <eigen3/Eigen/Dense>
#include <gmsh/Gmsh.h>
#include <stdlib.h>
#include <iostream>

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
        std::cout << std::endl;
        std::string ent_name;

        gmsh::model::setEntityName(
            ent_dim,
            ent_tag,
            ent_name);

        std::cout << "\tName: " << ent_name << std::endl;

        std::vector<int> phys_tags;
        gmsh::model::getPhysicalGroupsForEntity(
            ent_dim,
            ent_tag,
            phys_tags);
        std::cout << "\tPhysical Groups: ";
        for(int phys_idx = 0; phys_idx < phys_tags.size(); phys_idx++)
        {
            std::cout << phys_tags[phys_idx] << "\t";
        }
        std::cout << std::endl;
        std::string ent_type;
        gmsh::model::getType(ent_dim,ent_tag,ent_type);
        std::cout << "\tType: " << ent_type << std::endl;

        int parent_dim,parent_tag;
        gmsh::model::getParent(ent_dim,ent_tag,parent_dim,parent_tag);
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
            phys_name);
        std::cout << "\tPhysical Group Name: ";
        std::cout << phys_name << std::endl;

    }

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
