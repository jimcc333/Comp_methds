#include <iostream>
#include <map>

#include <eigen3/Eigen/Core>
#include <vector>
#include <math.h>

using namespace std;

class IsoInfo {
public:
    IsoInfo(unsigned int egroups, unsigned int f_order, unsigned int s_order, string iso_name);

    string name;
    unsigned int egroups;
    unsigned int f_order;
    unsigned int s_order;

    Eigen::VectorXf total;
    Eigen::MatrixXf ffactor; // ffactor[order][energy]
    Eigen::VectorXf chi;
    Eigen::VectorXf nufission;
    vector< Eigen::MatrixXf > skernel; // skernel[order][from][to]

    void Print();
    void Read(string data_path);
};

class RegionInfo {
public:
    float thickness;
    float dx;
    map<string,float> NumDens;

    vector<float> total;
    vector< Eigen::MatrixXf > skernel; // skernel[order][from][to]

    void Print();
};

class ParamsHolder {
public:
    // Global values
    unsigned int egroups = 10;
    unsigned int f_order = 6;
    unsigned int s_order = 9;
    unsigned int ordinates = 8;
    string data_path = "./Data/";
    string input_path = "./Input/input.txt";
    string output_name = "output.txt";
    float conv_tol = 0.002;

    // Problem objects
    vector<RegionInfo> region;

    // Isotope names
    vector<string> manifest;

    // Ordinate info
    const float mu2[2] = {0.5774, -0.5774};
    const float we2[2] = {1, 1};
    const float leg2[2][9] = {
                                {1,	0.5774, 8.61E-05, -0.384850438,	-0.38893673, -0.096349372, 0.222121708,	0.320769457, 0.152916538,},
                                {1, -0.5774, 8.61E-05, 0.384850438, -0.38893673, 0.096349372, 0.222121708, -0.320769457, 0.152916538}};

    float mu[8] = {0.9603, 0.7967, 0.5255, 0.1834, -0.1834, -0.5255, -0.7967, -0.9603};
    float we[8] = {0.1012, 0.2224, 0.3137, 0.3627, 0.3627, 0.3137, 0.2224, 0.1012};
    float leg[8][9] = {
                            {1E+00,	9.6030000E-01,	8.8326414E-01,	7.7346425E-01,	6.3737790E-01,	4.8296180E-01,	3.1913015E-01,	1.5517401E-01,	1.6161993E-04,},
                            {1E+00,	7.9670000E-01,	4.5209634E-01,	6.9175250E-02,	-2.4262639E-01,	-4.0328100E-01,	-3.8685029E-01,	-2.2670874E-01,	-1.6634017E-04,},
                            {1E+00,	5.2550000E-01,	-8.5774625E-02,	-4.2545761E-01,	-3.2693048E-01,	3.1122542E-02,	3.0242605E-01,	2.6846975E-01,	-9.6184884E-05,},
                            {1E+00,	1.8340000E-01,	-4.4954666E-01,	-2.5967810E-01,	2.5381631E-01,	2.9153232E-01,	-1.1349071E-01,	-2.8853978E-01,	8.2754120E-05,},
                            {1E+00,	-1.8340000E-01,	-4.4954666E-01,	2.5967810E-01,	2.5381631E-01,	-2.9153232E-01,	-1.1349071E-01,	2.8853978E-01,	8.2754120E-05,},
                            {1E+00,	-5.2550000E-01,	-8.5774625E-02,	4.2545761E-01,	-3.2693048E-01,	-3.1122542E-02,	3.0242605E-01,	-2.6846975E-01,	-9.6184884E-05,},
                            {1E+00,	-7.9670000E-01,	4.5209634E-01,	-6.9175250E-02,	-2.4262639E-01,	4.0328100E-01,	-3.8685029E-01,	2.2670874E-01,	-1.6634017E-04,},
                            {1E+00,	-9.6030000E-01,	8.8326414E-01,	-7.7346425E-01,	6.3737790E-01,	-4.8296180E-01,	3.1913015E-01,	-1.5517401E-01,	1.6161993E-04}};


    // Distributed source strength
    float init_source;

    void Print();
    void ReadIP();
    void BuildReg(vector<IsoInfo> &isos);
};

class Phi {
public:
    Phi(ParamsHolder &params);

    unsigned int tot;

    vector< vector < vector<float> > > flux; // [mesh][ordinate][energy]
    vector<float> distance;         // mesh-to-distance mapper
    vector<unsigned int> itoreg;    // mesh-to-region mapper, region indexed from zero
    vector< vector < vector<float> > > source; // [mesh][ordinate][energy]

    void PrintFlux();
    void PrintFlux(unsigned int ordinate, unsigned int group);
    void PrintSource();
    void SweepLR(ParamsHolder &params);
    void SweepRL(ParamsHolder &params);
    void CalcSource(ParamsHolder &params);
    void AddFlux(vector< vector < vector<float> > > addedflux);
    bool ConvCheck(vector< vector < vector<float> > > &total, float tolerance);
};

















