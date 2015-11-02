#include <iostream>

#include "classes.h"

void RegionInfo::Print() {
    cout << "Region thickness: " << thickness << ". Isos and number densities:" << endl;

    for(map<string, float>::const_iterator it = NumDens.begin(); it != NumDens.end(); ++it) {
        cout << it->first << " " << it->second << endl;
    }
}

IsoInfo::IsoInfo(unsigned int egroups, unsigned int f_order, unsigned int s_order, string iso_name) {
    name = iso_name;

    total.resize(egroups,1);
    total.setZero();

    ffactor.resize(egroups, f_order);
    ffactor.setZero();

    chi.resize(egroups,1);
    chi.setZero();

    nufission.resize(egroups,1);
    nufission.setZero();

    Eigen::MatrixXf temp_matrix;
    temp_matrix.resize(egroups, egroups);
    temp_matrix.setZero();
    for(unsigned int i = 0; i < s_order; i++) {
        skernel.push_back(temp_matrix);
    }

}


void IsoInfo::Print() {
    cout << "Total:" << endl << total << endl;
    cout << "F factor:" << endl << ffactor << endl;
    cout << "Chi:" << endl << chi << endl;
    cout << "NuFission:" << endl << nufission << endl;
    cout << "Skernel[0]:" << endl << skernel[0] << endl;
}

// Prints what it's holding to terminal
void ParamsHolder::Print() {
    cout << " Database folder: " << data_path << endl
         << " Input file:      " << input_path << endl
         << " Energy groups:   " << egroups << endl;
}

// Reads isotope databases
void ParamsHolder::ReadIP() {
    ifstream input(input_path);
    if(!input.is_open()) {cout << "Couldn't open input file!" << endl; return;}

    string line;
    string name;
    float value;
    unsigned int counter = 0;

    while(getline(input, line)) {
        if(!line.compare("Manifest")) {
            cout << "..Reading manifest" << endl;
                while(getline(input, line), !line.empty()) {
                    // Read a list of all the isos in problem

                    istringstream iss(line);
                    iss >> name;

                    manifest.push_back(name);

                }
        }

        if(!line.compare("Region")) {
                counter++;
                cout << "..Reading region " << counter << endl;

                RegionInfo temp_region;

                while(getline(input, line), line[0] != 'X') {
                    // Read each line until the thickness (X)

                    istringstream iss(line);
                    iss >> name >> value;

                    temp_region.NumDens[name] = value;

                }

                // Store the thickness value
                istringstream iss(line);
                iss >> name >> value;

                temp_region.thickness = value;

                region.push_back(temp_region);
        }

    }
    for(int i = 0; i < region.size(); i++) {
        cout << "Region " << i+1 << " ";
        region[i].Print();
    }
}
