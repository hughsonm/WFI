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
private:
    void CalculateTriAreas();
    void CalculateTriCentroids();
public:
    void buildTriangulation(
        std::string filename,
        bool verbose=false
    );
    void buildDomainGreen(
        Eigen::MatrixXcd & G,
        std::complex<double> k2_b
    );
    void buildDataGreen(
        Eigen::MatrixXcd & G,
        std::complex<double> k2_b,
        const Eigen::MatrixXd & locations
    );
    Eigen::MatrixXd points;
    Eigen::MatrixXd tri;
    Eigen::MatrixXd centroids;
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
    Eigen::MatrixXcd G_b_domain;
    Eigen::MatrixXcd G_b_data;
    Eigen::MatrixXcd L_domain;
    Eigen::MatrixXcd Chi;
    Eigen::PartialPivLU<Eigen::MatrixXcd> LU_L;
public:
    Chamber(std::string meshfile);
    void addTarget(std::string targetfile);
    void setupAntennas(std::string antennafile);
    void setupProbes(std::string probefile);
    void setFrequency(double freq);
    void calcDomainEzTot(void);
    void calcDataEzTot(void);
    std::vector<Antenna> antennas;
    Eigen::MatrixXd probe_points;
    Mesh mesh;
    Eigen::VectorXcd eps_r;
    Eigen::VectorXcd k2_f;
    std::complex<double> k2_b;
    Eigen::MatrixXcd Ez_inc,Ez_tot,Ez_sct;
    Eigen::MatrixXcd Ez_inc_d,Ez_tot_d,Ez_sct_d;
    bool Ez_inc_ready,Ez_tot_ready,Ez_sct_ready;
    bool G_b_domain_ready,G_b_data_ready;
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
