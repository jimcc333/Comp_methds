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
    vector< Eigen::MatrixXf > skernel; // skernel[from][to]

    void Print();
    void Read(string data_path);
};

class RegionInfo {
public:
    float thickness;
    float dx;
    map<string,float> NumDens;

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

    // Problem objects
    vector<RegionInfo> region;

    // Isotope names
    vector<string> manifest;

    // Ordinate info
    const float mu2[2] = {0.5774, -0.5774};
    const float we2[2] = {1, 1};
    const float mu8[8] = {0.9603, 0.7967, 0.5255, 0.1834, -0.1834, -0.5255, -0.7967, -0.9603};
    const float we8[8] = {0.1012, 0.2224, 0.3137, 0.3627, 0.3627, 0.3137, 0.2224, 0.1012};

    // Distributed source strength
    float source;

    void Print();
    void ReadIP();
};


class Phi {
public:
    Phi(ParamsHolder &params, vector<IsoInfo> &isos);

    vector< vector < vector<float> > > flux; // [mesh][ordinate][energy]
    vector<float> distance;



};
















