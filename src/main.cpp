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

        }
        cout << "Data for case: " << endl;
        params.Print();
    }

    // Read in all isotope data
    params.ReadIso();

    return 0;
}
