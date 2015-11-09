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
    string name, name2;
    float value;
    unsigned int counter = 0;
    float tot_thickness = 0;

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

            ///TODO this should to be fixed lol
            if(value == 2) {
                mu[0] = mu2[0];
                mu[1] = mu2[1];
                we[0] = we2[0];
                we[1] = we2[1];
                for(int i = 0; i < 9; i++) {
                    leg[0][i] = leg2[0][i];
                    leg[1][i] = leg2[1][i];
                }
            }
            cout << "..Number of ordinates: " << value << endl;
        }

        if(!line.compare(0,6,"source")) {
            istringstream iss(line);
            iss >> name >> value;

            init_source = value;
            cout << "..Distributed source strength: " << value << endl;
        }

        if(!line.compare(0,6,"output")) {
            istringstream iss(line);
            iss >> name >> name2;
            output_name = name2 + ".txt";
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
                tot_thickness += value;

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
    cout << ".Input file read. Total thickness: " << tot_thickness << "." << endl;
    /*
    for(int i = 0; i < region.size(); i++) {
        cout << "Region " << i+1 << " ";
        region[i].Print();
    }*/
}

// Builds the region specific material data
void ParamsHolder::BuildReg(vector<IsoInfo> &isos) {
    // For each region r
    for(int r = 0; r < region.size(); r++) {
        // Resize variables
        region[r].total.resize(egroups, 0);

        Eigen::MatrixXf temp_matrix;
        temp_matrix.resize(egroups, egroups);
        temp_matrix.setZero();
        region[r].skernel.resize(s_order, temp_matrix);

        ///TODO the following two nested for loops can be merged to lower code overhead

        // Calculate total cross section
        for(int e = 0; e < egroups; e++) {
        // Combine isos for each energy group e

            for(map<string,float>::iterator it = region[r].NumDens.begin(); it != region[r].NumDens.end(); ++it) {
            // For every isotope in region r

                for(unsigned int iso = 0; iso < isos.size(); iso++) {
                // For each isotope in isos

                    if(!isos[iso].name.compare(it->first + ".xs")){

                        region[r].total[e] += it->second * isos[iso].total[e] * 1E-24;
                    }
                }
            }
            //cout << region[r].total[e] << " ";
        }
         //cout << endl;

        // Calculate scattering kernel
        for(int o = 0; o < s_order; o++) {
        // For order o

            for(map<string,float>::iterator it = region[r].NumDens.begin(); it != region[r].NumDens.end(); ++it) {
            // For every isotope in region r

                for(unsigned int iso = 0; iso < isos.size(); iso++) {
                // For each isotope in isos

                    if(!isos[iso].name.compare(it->first + ".xs")){
                        region[r].skernel[o] += isos[iso].skernel[o] * it->second * 1E-24;
                    }
                }
            }
        }
    }
}


Phi::Phi(ParamsHolder &params) {
    // Create empty ordinate and energy flux to push, all values zero
    const vector< vector<float> > init(params.ordinates, vector<float>(params.egroups, 0));

    // Determine number of mesh points
    unsigned int tot_mesh = 0;
    unsigned int reg_mesh;

    //  First point is at distance zero
    distance.push_back(0);
    flux.push_back(init);
    itoreg.push_back(0);

    // Goes through each region to initiate phi and distance
    for(int r = 0; r < params.region.size(); r++) {
        reg_mesh = floor(params.region[r].thickness / params.region[r].dx);
        reg_mesh *= 2; // To include half-points

        const float dx2 = params.region[r].dx / 2;
        const unsigned int start = distance.size();

        distance.resize(distance.size() + reg_mesh);
        flux.resize(flux.size() + reg_mesh, init);
        itoreg.resize(itoreg.size() + reg_mesh, r);

        for(int i = start; i < reg_mesh + start; i++) {
            distance[i] = distance[i-1]+dx2;
        }
    }
    tot = distance.size();

    // Initiate source
    source.resize(tot, init);
    for(int i = 0; i < tot; i++) {
        for(int n = 0; n < params.ordinates; n++) {
            source[i][n][0] = params.init_source * params.we[n] / 2.;
        }
    }

/*
    for(int i = 0; i < distance.size(); i++){
        cout << i << " " << itoreg[i] << " " << distance[i] << " " << flux[i][0][0] << endl;
    }*/
}

void Phi::PrintFlux() {
    cout << "------- Phi -------" << endl;
    for(unsigned int i = 0; i < tot; i++) {
    // For mesh point i
        cout << "Mesh " << i << endl;
        for(unsigned int n = 0; n < flux[0].size(); n++) {
        // For right-pointed ordinate n

            for(unsigned int g = 0; g < flux[0][0].size(); g++) {
                    cout << flux[i][n][g] << " ";
            }
            cout << endl;
        }
    }
}

void Phi::PrintFlux(unsigned int ordinate, unsigned int group) {
    cout << "Flux of ordinate " << ordinate << " and group " << group << ":" << endl;
    for(unsigned int i = 0; i < tot; i++) {
        cout << flux[i][ordinate][group] << endl;
    }
}

void Phi::PrintSource() {
    cout << "------- Source -------" << endl;
    for(unsigned int i = 0; i < tot; i++) {
    // For mesh point i
        cout << "Mesh " << i << endl;
        for(unsigned int n = 0; n < flux[0].size(); n++) {
        // For right-pointed ordinate n

            for(unsigned int g = 0; g < flux[0][0].size(); g++) {
                    cout << source[i][n][g] << " ";
            }
            cout << endl;
        }
    }
}


// Sweeps from left to right
void Phi::SweepLR(ParamsHolder &params) {
    for(unsigned int i = 1; i < tot; i++) {
    // For mesh point i

        for(unsigned int n = 0; n < params.ordinates/2; n++) {
        // For right-pointed ordinate n

            for(unsigned int g = 0; g < params.egroups; g++) {
            // For group g
                if(i%2 == 0) {
                    flux[i][n][g] = 2.*flux[i-1][n][g] - flux[i-2][n][g];
                } else {
                ///todo check source being multiplied by delta
                    flux[i][n][g] = (params.region[itoreg[i]].dx * source[i][n][g] * params.we[n] / 2. + 2.*params.mu[n]*flux[i-1][n][g])
                                    / (2.*params.mu[n] + params.region[itoreg[i]].dx*params.region[itoreg[i]].total[g]);
                }
            }

        }
    }
}

// Sweeps from right to left
void Phi::SweepRL(ParamsHolder &params) {
    for(int i = tot-2; i >= 0; i--) {
    // For mesh point i

        for(unsigned int n = params.ordinates/2; n < params.ordinates; n++) {
        // For left-pointed ordinate n

            for(unsigned int g = 0; g < params.egroups; g++) {
            // For group g
                if(i%2 == 0) {
                    flux[i][n][g] = 2.*flux[i+1][n][g] - flux[i+2][n][g];
                } else {
                ///todo check source being multiplied by delta
                    flux[i][n][g] = (params.region[itoreg[i]].dx * source[i][n][g] * params.we[n] / 2. - 2.*params.mu[n]*flux[i+1][n][g])
                                    / (-2.*params.mu[n] + params.region[itoreg[i]].dx*params.region[itoreg[i]].total[g]);
                }
            }
        }
    }
}


void Phi::CalcSource(ParamsHolder &params) {
    float inner = 0;
    float middle = 0;
    float outer = 0;

    for(int i = 1; i < tot; i+=2) {
    // For mesh point i

        for(unsigned int n = 0; n < params.ordinates; n++) {
        // For source ordinate n

            for(unsigned int g = 0; g < params.egroups; g++) {
            // For source group g

                for(unsigned int l = 0; l < params.s_order; l++) {
                // For legengre order l

                    for(unsigned int k = 0; k < params.egroups; k++) {
                    // For flux group k

                        for(unsigned int p = 0; p < params.ordinates; p++) {
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
}


void Phi::AddFlux(vector< vector < vector<float> > > addedflux) {
    for(int i = 0; i < addedflux.size(); i++) {
        for(int j = 0; j < addedflux[i].size(); j++) {
            for(int k = 0; k < addedflux[i][j].size(); k++) {
                flux[i][j][k] += addedflux[i][j][k];
            }
        }
    }
}




































