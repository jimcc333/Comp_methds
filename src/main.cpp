// Computation methods class project
// Cem Bagdatli
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <boost/thread.hpp>

#include "classes.cpp"


using namespace std;

/***

params: hold problem variables and region info
isos:   isotope database vector



***/

void OutputGen(Phi &phi, ParamsHolder &params) {
    ifstream infile(params.output_name);
    if(infile.good()) {
        cout << params.output_name << " already exists. Overwriting old one." << endl;
    }

    float tot_flux[params.egroups];
    float flux = 0;
    float tot_left = 0;
    float tot_right = 0;
    ofstream output;
    output.open(params.output_name);

    output << "Transport code output file." << endl << endl;

    output << "Total source generations: " << params.tot_iter << endl;

    for(int g = 0; g < params.egroups; g++) {
        for(int n = 0; n < params.ordinates/2; n++) {
            tot_right += phi.flux[phi.tot-2][n][g] * abs(params.mu[n]) / 2; // last mesh in half point
        }

        for(int n = params.ordinates/2; n < params.ordinates; n++) {
            tot_left += phi.flux[1][n][g] * abs(params.mu[n]) / 2; // first mesh is half point
        }
    }

    output << "Right leakage: " << tot_right << " neutrons, " << tot_right/params.tot_source*100 << "%" << endl;
    output << "Left leakage : " << tot_left << " neutrons, " << tot_left/params.tot_source*100 << "%" << endl;

    output << endl << "_____Total Flux_____" << endl;
    output << "dist | flux" << endl;
    for(int i = 1; i < phi.tot; i+=2) {
        flux = 0;
        for(int g = 0; g < params.egroups; g++) {
            for(int n = 0; n < params.ordinates; n++) {
                flux += phi.flux[i][n][g];
            }
        }
        output << phi.distance[i] << " " << flux << endl;
    }

    output << endl << "_____Absorption rate_____" << endl;
    output << "dist | absorption" << endl;
    for(int i = 1; i < phi.tot; i+=2) {
        flux = 0;
        for(int g = 0; g < params.egroups; g++) {
            for(int n = 0; n < params.ordinates; n++) {
                flux += phi.flux[i][n][g] * params.region[phi.itoreg[i]].total[g];
            }
        }
        output << phi.distance[i] << " " << flux << endl;
    }

    output << endl << endl << "___Angle int. flux___" << endl;
    output << "dist | flux_g1 | flux_g2 | ..." << endl;
    for(int i = 1; i < phi.tot; i+=2) {
        output << phi.distance[i];
        for(int g = 0; g < params.egroups; g++) {
            tot_flux[g] = 0;
            for(int n = 0; n < params.ordinates; n++) {
                tot_flux[g] += phi.flux[i][n][g] * params.we[n] / 2;
            }
            output << " " << tot_flux[g];
        }
        output << endl;
    }

    output << endl << endl << "_____Psi_____" << endl;
    output << "dist" << endl << "matrix[ordinate X group]" << endl << endl;
    for(int i = 1; i < phi.flux.size(); i+=2) {
        output << phi.distance[i] << endl;
        for(int n = 0; n < phi.flux[i].size(); n++) {
            for(int g = 0; g < phi.flux[i][n].size(); g++) {
                output << phi.flux[i][n][g] << " ";
            }
            output << endl;
        }
        output << endl;
    }

    output << endl << endl << "Inputs for the case: " << endl;
    output << " Database folder: " << params.data_path << endl
         << " Input file:      " << params.input_path << endl
         << " Energy groups:   " << params.egroups << endl
         << " F order:         " << params.f_order << endl
         << " S order:         " << params.s_order << endl << endl;

    for(int r = 0; r < params.region.size(); r++) {
        output << "Region " << r+1 << ":" << endl;
        output << "  Thickness: " << params.region[r].thickness << ". Isos and number densities:" << endl;

        for(map<string, float>::const_iterator it = params.region[r].NumDens.begin(); it != params.region[r].NumDens.end(); ++it) {
            output << "  " << it->first << " " << it->second << endl;
        }
    }

    output.close();

    cout << "Completed writing results." << endl;

    return;

}

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

            if(string(argv[arg]) == "ipath" || string(argv[arg]) == "i") {
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

    // Builds the necessary parameters for each region
    params.BuildReg(isos);

    // Build flux vector
    Phi phi1(params);
    Phi total(params);

    // Determine total source
    params.tot_source = 0;
    for(int i = 1; i < phi1.source.size(); i+=2) {
        for(int j = 0; j < phi1.source[i].size(); j++) {
            for(int k = 0; k < phi1.source[i][j].size(); k++) {
                params.tot_source += phi1.source[i][j][k];
            }
        }
    }
    cout << "tot source: " << params.tot_source << endl;

    // First sweep using given source
    phi1.SweepLR(params);
    phi1.SweepRL(params);

    total.AddFlux(phi1.flux);

    unsigned int counter;
    cout << "Starting solution..." << endl;
    for(counter = 0; counter < 200; counter++) {
        // Progress output
        cout << "Iteration: " << counter << "\r";
        cout.flush();

        // Calculate new source
        phi1.CalcSource(params);
        // Sweep
        phi1.SweepLR(params);
        phi1.SweepRL(params);
        // Check convergence
        if(phi1.ConvCheck(total.flux, params.conv_tol)) {
            cout << endl << "Calculation complete after " << counter << " iterations!" << endl;
            params.tot_iter = counter;
            counter = 500;
        }
        // Add to total
        total.AddFlux(phi1.flux);
    }

    cout << "Generating output file " << params.output_name << endl;
    OutputGen(total, params);

    return 0;
}

















