// Computation methods class project
// Cem Bagdatli
#include <iostream>
#include <fstream>
#include "classes.cpp"
#include <vector>

using namespace std;

/***

params: hold problem variables and region info
isos:   isotope database vector



***/

int main(int argc, char* argv[]) {
    // Defaults
    ParamsHolder params;


    // Handle console inputs
    if(argc == 1) {
        cout << "No inputs. Assuming:" << endl;
        params.Print();


    } else {
        unsigned int argument_count = argc - 1;

        for(int arg = 1; arg < argc; arg++) {
            ///TODO check to see if next value exists
            if(string(argv[arg]) == "groups") {
                if(stoi(argv[1+arg]) > 10) {
                    cout << "Error! Max 10 groups supported!" << endl;
                    return 1;
                }
                params.egroups = stoi(argv[++arg]);
            }

            if(string(argv[arg]) == "dpath") {
                params.data_path = string(argv[++arg]);
            }

            if(string(argv[arg]) == "ipath") {
                params.input_path = string(argv[++arg]);
            }

            if(string(argv[arg]) == "f") {
                if(stoi(argv[1+arg]) > 6) {
                    cout << "Error! Max 6 f order supported!" << endl;
                    return 1;
                }
                params.f_order = stoi(argv[++arg]);
            }

            if(string(argv[arg]) == "s") {
                if(stoi(argv[1+arg]) > 9) {
                    cout << "Error! Max 9 s order supported!" << endl;
                    return 1;
                }
                params.s_order = stoi(argv[++arg]);
            }

        }
        cout << "Data for case: " << endl;
        params.Print();
    }

    // Read input data
    params.ReadIP();

    // Initiate isotope data object
    cout << "..Reading database." << endl;
    vector<IsoInfo> isos;
    for(unsigned int i = 0; i < params.manifest.size(); i++) {
        cout << "...Reading " << params.manifest[i] << " input file." << endl;
        IsoInfo temp_iso(params.egroups, params.f_order, params.s_order, params.manifest[i]);
        temp_iso.Read(params.data_path);
        isos.push_back(temp_iso);
    }
    cout << "..Database read." << endl;

    // Build flux vector
    Phi phi1(params, isos);

    //isos.Print();

    return 0;
}

















