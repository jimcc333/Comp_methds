// Computation methods class project
// Cem Bagdatli
#include <iostream>
#include <fstream>
#include "classes.cpp"
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
    // Defaults
    ParamsHolder params;


    // Handle inputs
    if(argc == 1) {
        cout << "No inputs. Assuming:" << endl;
        params.Print();


    } else {
        unsigned int argument_count = argc - 1;

        for(int arg = 1; arg < argc; arg++) {
            if(string(argv[arg]) == "groups") {
                ///TODO check to see if next value exists
                params.egroups = stoi(argv[++arg]);
            }

            if(string(argv[arg]) == "dpath") {
                params.data_path = string(argv[++arg]);
            }

            if(string(argv[arg]) == "ipath") {
                params.input_path = string(argv[++arg]);
            }

            if(string(argv[arg]) == "f") {
                params.f_order = stoi(argv[++arg]);
            }

            if(string(argv[arg]) == "s") {
                params.s_order = stoi(argv[++arg]);
            }

        }
        cout << "Data for case: " << endl;
        params.Print();
    }

    // Read input data
    params.ReadIP();

    // Initiate isotope data object
    for(unsigned int i = 0; i < params.manifest.size(); i++) {

    }
    IsoInfo isos(params.egroups, params.f_order, params.s_order, "bla");


    //isos.Print();

    return 0;
}

















