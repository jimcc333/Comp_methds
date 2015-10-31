#include <iostream>
#include<map>

#include <eigen3/Eigen/Core>
#include <vector>

using namespace std;

class RegionInfo {
public:
    float thickness;
    map<string,float> NumDens;

    void Print();
};

class ParamsHolder {
public:
    // Global values
    unsigned int egroups = 10;
    string data_path = "./Data/";
    string input_path = "./Input/input.txt";

    // Problem objects
    RegionInfo region;

    void Print();
    void ReadIso();
};


