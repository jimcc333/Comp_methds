#include <iostream>

#include "classes.h"

void RegionInfo::Print() {
    cout << "Region thickness: " << thickness << ". Isos and number densities:" << endl;

    for(map<string, float>::const_iterator it = NumDens.begin(); it != NumDens.end(); ++it) {
        cout << it->first << " " << it->second << endl;
    }
}


// Prints what it's holding to terminal
void ParamsHolder::Print() {
    cout << " Database folder: " << data_path << endl
         << " Input file:      " << input_path << endl
         << " Energy groups:   " << egroups << endl;
}

// Reads isotope databases
void ParamsHolder::ReadIso() {
    ifstream input(input_path);
    if(!input.is_open()) {cout << "Couldn't open input file!" << endl; return;}

    string line;
    string name;
    float value;

    while(getline(input, line)) {
        if(!line.compare("Input begin")) {
                cout << "Reading region 1" << endl;

                while(getline(input, line), line[0] != 'X') {
                    // Read each line until the thickness (X)

                    istringstream iss(line);
                    iss >> name >> value;

                    region.NumDens[name] = value;

                }

                // Store the thickness value
                istringstream iss(line);
                iss >> name >> value;

                region.thickness = value;
        }
        ///TODO add next region parsing method

    }
    region.Print();
}
