#include <iostream>

#include "classes.h"

IsoInfo::IsoInfo(unsigned int egroups, unsigned int f_order, unsigned int s_order, string iso_name) {
    name = iso_name;
    this->egroups = egroups;
    this->f_order = f_order;
    this->s_order = s_order;

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

void IsoInfo::Read(string data_path) {
    // Open isotope file
    ifstream input(data_path + name);
    if(!input.is_open()) {
        cout << endl << "(IsoInfo)Couldn't open input file " << data_path + name << endl;
        exit(1);
    }

    string line;
    int group, togroup;
    float value;
    float value1, value2, value3, value4, value5, value6, value7, value8, value9;

    while(getline(input, line)) {

        if(!line.compare("TOTAL")) {
            while(getline(input, line), !line.empty()) {

                istringstream iss(line);
                iss >> group >> value;

                total(egroups - group) = value;
            }
        }


        if(!line.compare("FFACTOR")) {
            while(getline(input, line), !line.empty()) {

                istringstream iss(line);
                iss >> group >> value1 >> value2 >> value3 >> value4 >> value5 >> value6;

                ///TODO not order 6 support
                ffactor(egroups - group, 0) = value1;
                ffactor(egroups - group, 1) = value2;
                ffactor(egroups - group, 2) = value3;
                ffactor(egroups - group, 3) = value4;
                ffactor(egroups - group, 4) = value5;
                ffactor(egroups - group, 5) = value6;
            }
        }

        if(!line.compare("CHI")) {
            while(getline(input, line), !line.empty()) {

                istringstream iss(line);
                iss >> group >> value;

                chi(egroups - group) = value;
            }
        }

        if(!line.compare("NUFISSION")) {
            while(getline(input, line), !line.empty()) {

                istringstream iss(line);
                iss >> group >> value;

                nufission(egroups - group) = value;
            }
        }

        if(!line.compare("SKERNEL")) {
            while(getline(input, line)) {

                istringstream iss(line);
                iss >> group >> togroup >> value1 >> value2 >> value3
                    >> value4 >> value5 >> value6 >> value7 >> value8 >> value9;

                ///TODO not order 8 support
                skernel[0](egroups - group, egroups - togroup) = value1;
                skernel[1](egroups - group, egroups - togroup) = value2;
                skernel[2](egroups - group, egroups - togroup) = value3;
                skernel[3](egroups - group, egroups - togroup) = value4;
                skernel[4](egroups - group, egroups - togroup) = value5;
                skernel[5](egroups - group, egroups - togroup) = value6;
                skernel[6](egroups - group, egroups - togroup) = value7;
                skernel[7](egroups - group, egroups - togroup) = value8;
                skernel[8](egroups - group, egroups - togroup) = value9;

            }
        }

    }
}

void IsoInfo::Print() {
    cout << "Total:" << endl << total << endl << endl;
    cout << "F factor:" << endl << ffactor << endl << endl;
    cout << "Chi:" << endl << chi << endl << endl;
    cout << "NuFission:" << endl << nufission << endl << endl;
    cout << "Skernel[0]:" << endl << skernel[0] << endl;
}

void RegionInfo::Print() {
    cout << "  Region thickness: " << thickness << ". Isos and number densities:" << endl;

    for(map<string, float>::const_iterator it = NumDens.begin(); it != NumDens.end(); ++it) {
        cout << "  " << it->first << " " << it->second << endl;
    }
}

// Prints what it's holding to terminal
void ParamsHolder::Print() {
    cout << " Database folder: " << data_path << endl
         << " Input file:      " << input_path << endl
         << " Energy groups:   " << egroups << endl
         << " F order:         " << f_order << endl
         << " S order:         " << s_order << endl;
}

// Reads input file
void ParamsHolder::ReadIP() {
    // Open file
    cout << "..Opening input file" << endl;
    ifstream input(input_path);
    if(!input.is_open()) {
        cout << endl << "(ParamsHolder)Couldn't open input file " << input_path << endl;
        exit(1);
    }

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

        if(!line.compare(0,9,"ordinates")) {
            istringstream iss(line);
            iss >> name >> value;

            ordinates = value;
            if(value != 8 && value != 2) {
                cout << "Code does not support " << value << " ordinates. Exiting :(" << endl;
                exit(1);
            }
            cout << "..Number of ordinates: " << value << endl;
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

                // Store delta
                getline(input, line);
                istringstream iss2(line);
                iss2 >> name >> value;

                if(name.compare("dx")) {
                    cout << " Error parsing region " << counter << " thickness." << endl;
                    exit(1);
                } else {
                    temp_region.dx = value;
                }

                region.push_back(temp_region);
        }

    }
    /*
    for(int i = 0; i < region.size(); i++) {
        cout << "Region " << i+1 << " ";
        region[i].Print();
    }*/
}






