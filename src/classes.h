#include <iostream>

#include <eigen3/Eigen/Core>
#include <vector>

using namespace std;

class ParamsHolder {
public:
    unsigned int egroups = 10;
    string data_path = "./Data/";
    string input_path = "./Input/input.txt";

    void Print();
    void ReadIso();
};


class Isotope {
    Isotope();
};
