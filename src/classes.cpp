#include <iostream>

#include "classes.h"

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
    while(getline(input, line)) {
        if(!line.compare("NUMBER DENSITIES")) {cout << "YAY" << endl;}

        cout << line[0] << endl;
    }

}
