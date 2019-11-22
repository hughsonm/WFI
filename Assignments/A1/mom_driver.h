#ifndef MOM_DRIVER_H
#define MOM_DRIVER_H
#include <eigen3/Eigen/Eigen>
void WriteMatrixToFile(
    std::string filename,
    const Eigen::MatrixXcd & matrix,
    bool append=false
);

void WriteMatrixToFile(
    std::string filename,
    const Eigen::MatrixXd & matrix,
    bool append=false
);


void WriteVectorToFile(
    std::string filename,
    const Eigen::VectorXd & vec,
    bool append=false
);


void WriteVectorToFile(
    std::string filename,
    const Eigen::VectorXcd & vec,
    bool append=false
);

void ReadMatrixFromFile(
    std::string filename,
    Eigen::MatrixXcd & matrix
);

void ReadMatrixFromFile(
    std::string filename,
    Eigen::MatrixXd & matrix
);

void ReadVectorFromFile(
    std::string filename,
    Eigen::VectorXcd & vec
);

void ReadVectorFromFile(
    std::string filename,
    Eigen::VectorXd & vec
);


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
    std::string meshfilename;
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
    bool ready{false};
};

class Field{
private:
    Eigen::MatrixXd * locs;
    Eigen::VectorXcd vals;
public:
    void setLocations(
        Eigen::MatrixXd & set_locs
    ){
        locs = &set_locs;
    };
    void setVals(
        const Eigen::VectorXcd & set_vals
    ){
        vals.resize(set_vals.size());
        vals = set_vals;
    }
    Eigen::MatrixXd * getLocations(void){
        return locs;
    };
    void getVals(Eigen::VectorXcd & dest){
        dest.resize(vals.size());
        dest = vals;
    };
    const Eigen::VectorXcd & getValRef(void){
      return vals;
    }
    void erase(void){
        vals.resize(locs->rows());
        vals.setOnes();
    };
    void WriteValsToFile(std::string filename){WriteVectorToFile(filename,vals);};
    bool ready{false};
};

void FormFieldMatrix(
    Eigen::MatrixXcd & M,
    std::vector<Field> & V
);

class Antenna
{
public:
    AntennaStyle_t style;
    Eigen::Vector3d location;
    Eigen::Vector3d direction;
    std::complex<double> coefficient;
    void getEz(
        Eigen::VectorXcd & Ez,
        const Eigen::MatrixXd & points,
        double frequency
    );
    bool ready{false};
};

class Probe
{
public:
    Eigen::Vector3d location;
    std::complex<double> getMeasurement(
        const Field & u,
        const Eigen::VectorXcd & data_green_row
    );
    bool ready{false};
};

class Target
{
public:
    Field eps_r;
    void getSqWaveNumberField(
        Field & k2,
        double frequency
    );
    bool ready{false};
};




class Chamber
{
public:
    Chamber(std::string meshfile);
    void setTarget(std::string targetfile);
    void setupAntennas(std::string antennafile);
    void setupProbes(std::string probefile);
    void setFrequency(double freq);
    void calcDomainEzInc(void);
    void calcDomainEzTot(void);
    void calcDataEzInc(void);
    void calcDataEzTot(void);
    void buildDataGreen(void);
    void A2Q3(
        Eigen::MatrixXcd & w_calc,
        Eigen::VectorXcd & X_calc
    );
    void A2Q5(
        Eigen::MatrixXcd & w_calc,
        Eigen::MatrixXcd & u_calc,
        Eigen::VectorXcd & X_calc,
        std::vector<std::vector<Eigen::Vector2d> > & curves,
        bool tikhonov
    );
    void A3P3(
        Eigen::MatrixXd& kpts,
        Eigen::VectorXcd& kvals,
        Eigen::VectorXcd& chi
    );
    void A3P4(std::string freqfile);
    void readMeasuredData(
        std::string datafile,
        double noise_pct
    );
    double frequency;
    std::complex<double> k2_b{0.0};
    Eigen::MatrixXcd G_b_domain;
    Eigen::MatrixXcd G_b_data;
    // L = (I+G*Chi);
    Eigen::MatrixXcd L_domain;
    Eigen::MatrixXcd Chi;
    Eigen::FullPivLU<Eigen::MatrixXcd> LU_L;
    std::vector<Antenna> antennas;
    std::vector<Probe> probes;
    Eigen::MatrixXd probe_points;
    Mesh mesh;
    Target target;

    std::vector<Field> Ez_inc,Ez_tot,Ez_sct;
    std::vector<Field> Ez_inc_d,Ez_tot_d,Ez_sct_d;
    std::vector<Field> Ez_tot_meas,Ez_sct_meas;
};


#endif
