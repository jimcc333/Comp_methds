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


void ThreadFunc(vector< vector < vector<float> > > &flux, vector< vector < vector<float> > > &source,
                ParamsHolder params, vector<unsigned int> &itoreg);

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
    cout << endl << "___Cem's Super-Awesome Transport Code___" << endl << endl;
    cout << "     .---------------------------." << endl;
    cout << "    /,--..---..---..---..---..--. `." << endl;
    cout << "   //___||___||___||___||___||___\\_|" << endl;
    cout << "   [j__ ######################## [_|" << endl;
    cout << "      \\============================|" << endl;
    cout << "   .==|  |'''||'''||'''||'''| |'''||" << endl;
    cout << "  /======'---''---''---''---'=|  =||" << endl;
    cout << "  |____    TRANSPORT    ____  | ==||" << endl;
    cout << "  //  \\\\         BUS   //  \\\\ |===||" << endl;
    cout << "  \\\\__//---------------\\\\__//-+---+'" << endl << endl;

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

            if(string(argv[arg]) == "threads" || string(argv[arg]) == "t") {
                if(stoi(argv[1+arg]) > 100) {
                    cout << "Error! Max 100 threads supported!" << endl;
                    return 1;
                }
                params.threads = stoi(argv[++arg]);
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
    cout << ".Total source strength: " << params.tot_source << endl;

    // First sweep using given source
    phi1.SweepLR(params);
    phi1.SweepRL(params);

    total.AddFlux(phi1.flux);

    unsigned int counter;
    cout << "Starting solution..." << endl;
    cout << "..Software threads: " << params.threads << endl;
    cout << "..Mesh points: " << phi1.flux.size() << endl;

    // Generate thread flux vectors
    vector< vector< vector < vector<float> > > > t_flux, t_source; // [thread][mesh][ordinate][energy]
    vector < vector<unsigned int> > t_itoreg;

    int mesh_size = floor(phi1.flux.size()/params.threads);
    if(mesh_size/2 != 0) {mesh_size--;} // Make sure its even

    for(int thread = 0; thread < params.threads - 1; thread++) {
        vector< vector < vector<float> > >::iterator it1 = phi1.flux.begin() + mesh_size*thread;
        vector< vector < vector<float> > >::iterator it2 = phi1.flux.begin() + mesh_size*(thread+1);
        vector< vector < vector<float> > > temp_flux(it1, it2);
        t_flux.push_back(temp_flux);
        t_source.push_back(temp_flux);

        vector<unsigned int>::iterator it3 = phi1.itoreg.begin() + mesh_size*thread;
        vector<unsigned int>::iterator it4 = phi1.itoreg.begin() + mesh_size*(thread+1);
        vector<unsigned int> temp_itoreg(it3, it4);
        t_itoreg.push_back(temp_itoreg);
    }
    vector< vector < vector<float> > >::iterator it1 = phi1.flux.begin() + mesh_size*(params.threads-1);
    vector< vector < vector<float> > >::iterator it2 = phi1.flux.end();
    vector< vector < vector<float> > > temp_flux(it1, it2);
    t_flux.push_back(temp_flux);
    t_source.push_back(temp_flux);

    vector<unsigned int>::iterator it3 = phi1.itoreg.begin() + mesh_size*(params.threads-1);
    vector<unsigned int>::iterator it4 = phi1.itoreg.end();
    vector<unsigned int> temp_itoreg(it3, it4);
    t_itoreg.push_back(temp_itoreg);

    // Start iteration
    for(counter = 0; counter < 200; counter++) {
        // Progress output
        cout << "Iteration: " << counter << "\r";
        cout.flush();

        // Update thread fluxes
        unsigned int thread = 0;
        unsigned int next = t_flux[0].size();
        unsigned int subtr = 0;
        for(unsigned int i = 1; i < phi1.flux.size(); i+=2) {
            if(i >= next) {
                thread++;
                next += t_flux[thread].size();
                subtr = i-1;
            }
            for(unsigned int n = 0; n < phi1.flux[0].size(); n++) {
                for(unsigned int g = 0; g < phi1.flux[0][n].size(); g++) {
                    //cout << "thread: " << thread+1 << " mesh: " << i << " subtr: " << subtr << " thread_i: "<< i-subtr << " next: " << next << " region: " << t_itoreg[thread][i-subtr] <<  endl;
                    t_flux[thread][i-subtr][n][g] = phi1.flux[i][n][g];
                }
            }
        }

        // Run in parallel to calculate source
        boost::thread_group threads;
        for(int thread = 0; thread < params.threads; thread++) {
            threads.create_thread(boost::bind(&ThreadFunc, boost::ref(t_flux[thread]), boost::ref(t_source[thread]),
                                          params, boost::ref(t_itoreg[thread])));
        }

        threads.join_all(); // calculated source and total threads

        // Update source from thread source
        thread = 0;
        next = t_flux[0].size();
        subtr = 0;
        for(unsigned int i = 1; i < phi1.flux.size(); i+=2) {
            if(i >= next) {
                thread++;
                next += t_flux[thread].size();
                subtr = i-1;
            }
            for(unsigned int n = 0; n < phi1.flux[0].size(); n++) {
                for(unsigned int g = 0; g < phi1.flux[0][n].size(); g++) {
                    //cout << "thread: " << thread+1 << " mesh: " << i << " subtr: " << subtr << " thread_i: " << i-subtr << " next: " << next << endl;
                    phi1.source[i][n][g] = t_source[thread][i-subtr][n][g];
                }
            }
        }

        //phi1.CalcSource(params);
        //phi1.PrintSource();
        //phi1.PrintFlux();

        // Sweep using new source
        phi1.SweepLR(params);
        phi1.SweepRL(params);

        total.AddFlux(phi1.flux);

        //total.PrintFlux();

        // Check convergence
        if(phi1.ConvCheck(total.flux, params.conv_tol)) {
            cout << endl << "Calculation complete after " << counter << " iterations!" << endl;
            params.tot_iter = counter;
            counter = 500;
        }
    }

    cout << "Generating output file " << params.output_name << endl;
    OutputGen(total, params);

    return 0;
}

// Updates source and total using flux
void ThreadFunc(vector< vector < vector<float> > > &flux, vector< vector < vector<float> > > &source,
                ParamsHolder params, vector<unsigned int> &itoreg) {
    // ! Assumes the first mesh point is a half-integer point

    const unsigned int I = flux.size();
    const unsigned int N = flux[0].size();
    const unsigned int G = flux[0][0].size();

    // Calculate source
    float inner = 0;
    float middle = 0;
    float outer = 0;

    for(unsigned int n = 0; n < N; n++) {
        for(unsigned int g = 0; g < G; g++) {
            for(unsigned int i = 1; i < I; i+=2) {

                for(unsigned int l = 0; l < params.s_order; l++) {
                // For legengre order l

                    for(unsigned int k = 0; k < G; k++) {
                    // For flux group k

                        for(unsigned int p = 0; p < N; p++) {
                        // For flux ordinate p
                            inner += params.we[p] * params.leg[p][l] * flux[i][p][k];
                        }
                        middle += inner * params.region[itoreg[i]].skernel[l](k,g);
                        inner = 0;
                    }
                    outer += middle * (2.*l + 1) * params.leg[n][l];
                    middle = 0;
                }
                source[i][n][g] = 0.5 * outer;
                outer = 0;
            }
        }
    }

    return;
}






















