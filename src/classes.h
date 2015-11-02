#include <iostream>
#include <map>

#include <eigen3/Eigen/Core>
#include <vector>

using namespace std;

class IsoInfo {
public:
    IsoInfo(unsigned int egroups, unsigned int f_order, unsigned int s_order, string iso_name);

    string name;

    Eigen::VectorXf total;
    Eigen::MatrixXf ffactor; // ffactor[order][energy]
    Eigen::VectorXf chi;
    Eigen::VectorXf nufission;
    vector< Eigen::MatrixXf > skernel; // skernel[from][to]

    void Print();
};

class RegionInfo {
public:
    float thickness;
    map<string,float> NumDens;
    //IsoInfo smear_iso;

    void Print();
};

class ParamsHolder {
public:
    // Global values
    unsigned int egroups = 10;
    unsigned int f_order = 6;
    unsigned int s_order = 9;
    string data_path = "./Data/";
    string input_path = "./Input/input.txt";

    // Problem objects
    vector<RegionInfo> region;

    // Isotope names
    vector<string> manifest;

    void Print();
    void ReadIP();
};






