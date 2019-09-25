#ifndef MOM_DRIVER_H
#define MOM_DRIVER_H

enum AntennaStyle_t
{
    LineSource,
    PlaneWave,
    Patch
};

class Antenna
{
    AntennaStyle_t style;
    Eigen::Vector3d location;
    Eigen::Vector3d direction;
    double frequency;
    std::complex<double> coefficient;
private:
public:
    Antenna();
    void getField(
        Eigen::MatrixXcd & Ez,
        const Eigen::MatrixXd & points
    );
    ~Antenna();
};

void EvaluateIncidentField(
    Eigen::VectorXcd & obs_Ez,
    const Eigen::VectorXd source_loc,
    const std::complex<double> source_coeff,
    const double frequency,
    const Eigen::MatrixXd & obs_pts
);

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
