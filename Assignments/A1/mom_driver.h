#ifndef MOM_DRIVER_H
#define MOM_DRIVER_H

#define EPSNAUGHT 8.8541848128E-12
#define MUNAUGHT 1.25663706212E-6
#define CNAUGHT 299792508.7882675

enum AntennaStyle_t
{
    LineSource,
    PlaneWave,
    Patch
};



class Mesh
{
public:
    void buildTriangulation(
        std::string filename,
        bool verbose=false
    );
    void buildDomainGreen(
        Eigen::MatrixXcd G,
        std::complex<double> k2_b
    );
    Eigen::MatrixXd points;
    Eigen::MatrixXd tri;
    Eigen::VectorXd centroids;
    Eigen::VectorXd areas;
};

class Antenna
{
public:
    AntennaStyle_t style;
    Eigen::Vector3d location;
    Eigen::Vector3d direction;
    double frequency;
    std::complex<double> coefficient;
    void getField(
        Eigen::VectorXcd & Ez,
        const Eigen::MatrixXd & points
    );
};

class Chamber
{
private:
    double frequency;
    Eigen::MatrixXcd Ez_inc,Ez_tot,Ez_sct;
    bool Ez_inc_ready,Ez_tot_ready,Ez_sct_ready,G_b_domain_ready;
    Eigen::MatrixXcd G_b_domain;
    Eigen::MatrixXcd L_domain;
    Eigen::MatrixXcd Chi;
    Eigen::PartialPivLU<Eigen::MatrixXcd> LU_L;
    void getEzInc(Eigen::MatrixXcd & Ezdest);
public:
    Chamber(std::string meshfile);
    void addTarget(std::string targetfile);
    void setupAntennas(std::string antennafile);
    void setFrequency(double freq);
    void getDomainEzTot(Eigen::MatrixXcd & Ezdest);
    std::vector<Antenna> antennas;
    Mesh mesh;
    Eigen::VectorXcd eps_r;
    Eigen::VectorXcd k2_f;
    std::complex<double> k2_b;
};

void WriteMatrixToFile(
    std::string filename,
    Eigen::MatrixXcd matrix,
    bool append=false
);

void WriteMatrixToFile(
    std::string filename,
    Eigen::MatrixXd matrix,
    bool append=false
);


void WriteVectorToFile(
    std::string filename,
    Eigen::VectorXd vec,
    bool append=false
);


void WriteVectorToFile(
    std::string filename,
    Eigen::VectorXcd vec,
    bool append=false
);


void ReadAntennaFile(
    Eigen::MatrixXd & locations,
    Eigen::VectorXcd &  coefficients,
    const std::string filename
);


void BuildTriangulation(
    std::string filename,
    Eigen::MatrixXd & tri,
    Eigen::MatrixXd & points,
    bool verbose = false
);


void CalculateTriAreas(
    Eigen::VectorXd & areas,
    const Eigen::MatrixXd & tri,
    const Eigen::MatrixXd & points,
    bool verbose = false
);

void CalculateTriCentroids(
    Eigen::MatrixXd & centroids,
    const Eigen::MatrixXd & tri,
    const Eigen::MatrixXd & points
);



void AssignRelativeConstitutives(
    Eigen::VectorXcd & eps_r,
    const Eigen::MatrixXd & tri,
    const std::string constfilename
);

void BuildDomainGreen(
    Eigen::MatrixXcd & G,
    const Eigen::MatrixXd & centroids,
    const Eigen::VectorXd & areas,
    double k2_b
);

void BuildDataGreen(
    Eigen::MatrixXcd & G,
    const Eigen::MatrixXd & centroids,
    const Eigen::VectorXd & areas,
    const Eigen::MatrixXd & rxlocations,
    double k2_b
);
#endif
