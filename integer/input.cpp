#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <map>
#include <iomanip>
using namespace std;
void fprint(vector<double>, ofstream& );
// Input Parameters

const string dimension = "2d";              // Simulation dimension
const string species[3] = {"e","Ar+","Ar"}; // Species name
const double weight = 1e-6;                 // macroparticle weight, n = Ns/weight/dV, weight *= n_0*V_0
const int Z = 1;                            // |ion charge/e-|
const double i_m = 6.63352090e-26;          // mass of ion(Ar) [kg]
const double n_0 = 1e15;                    // initial electron density
const double T_e0 = 3000.;                  // initial electron temperature [K]

const double L[6] = {0., 0.01, 0., 0.01, 0., 1.}; // Simulation size {xlo,xhi,ylo,yhi,zlo,zhi}
const int Ng[3] = {100, 100, 1};              // Number of grid points
const string conductor[2] = {"rectangle", "real"};
        
const double frequency = 13.56e6;           // RF frequency [Hz]
const int nstep_rf = 1000;                 // time steps in one RF cycle
const double T_in_rf = 1.0/frequency;       // total time in one cycle
const double DT_e = T_in_rf/nstep_rf;      // electron time step length
const int NCYCLE = 1000;                    // Total RF cycles  

const double Pressure = 10.;                // gas pressure [Pa]
const double Temperature = 300.;            // temperature [K]



const double CS_DE = 0.001;                 // energy interval [eV]
const int CS_TOT = 100000;                 // Total entries in cross section arrays
const double e_EXC_TH = 11.5;               // e + Ar -> e + Ar* (excitation threshold)
const double e_ION_TH = 15.8;               // e + Ar -> e + e + Ar+ (Ionization threshold)

// Basic constants

const double e = 1.60217662e-19;            // electron charge [C]
const double e_m = 9.10938356e-31;          // electron mass
const double k_B = 1.38064852e-23;          // Boltzmann constant
const double EPSILON0 = 8.85418781e-12;     // Vacuum dielectric constant
const double eV_TO_J = e;                   // eV <-> Joule conversion factor
const string seed = "default";

int main(int argc, char *argv[])
{   // Init non-dimension parameters
    const double lambda_De0 = pow(EPSILON0*k_B*Temperature/n_0/e/e, 0.5);
    const double omega_pe0 = pow(n_0*e*e/EPSILON0/e_m, 0.5);
    const int NSUB = static_cast<int> (pow(i_m/e_m, 0.5)*Z);                        // dti/dte = omega_pe0/omega_pi0
    const double n_gas0 = Pressure/k_B/Temperature;                 // Background neutral gas number density
    const double kTe0 = k_B * T_e0 / eV_TO_J;
    const double vte0 = pow(k_B*T_e0/e_m, 0.5);

    // Set Non-dimension forms

    cout << "Calculate the non-dimension form of Basic Parameters:... " << endl;
    cout << "   Processing Geometry parameters:... " << endl;
    double _n_gas = n_gas0/n_0;
    double _n_e = 1.;                                               // non dimension electron density n_e=n_0/n_0
    double _width_x, _width_y, _volinv, nd_L[6];
    for(int i=0; i<5; i++) { nd_L[i] = L[i]/lambda_De0; }
    nd_L[5] = L[5];
    _width_x = nd_L[0]/Ng[0], _width_y = nd_L[1]/Ng[1];
    _volinv = 1/(nd_L[0]*L[1]);
    size_t np_init = static_cast<int> (_n_e * _volinv * weight);    // Initial superparticle of e-/Ar+
    cout << "   Processing time step set:... " << endl;
    double _dt_i, _dt_e;
    _dt_e = DT_e * omega_pe0;
    _dt_i = NSUB * _dt_e;
    cout << "   Processing Total Cross Section(0.001,1000)[eV]:..." << endl;

    enum COLL {e_ELA, e_EXC, e_ION, i_ISO, i_BACK} ;
    vector<double> e_energy(CS_TOT), i_energy(CS_TOT);
    map<COLL, vector<double>> sigma;
    e_energy[0] = CS_DE;
    i_energy[0] = 2 * CS_DE;

    // Calculate electron cross sections(Data from Phelps A,PSST,1999)

    auto qela = [](auto en){ return 1e-20*(fabs(6.0 / pow(1.0 + (en/0.1) + pow(en/0.6,2.0), 3.3)
        - 1.1 * pow(en, 1.4) / (1.0 + pow(en/15.0, 1.2)) / sqrt(1.0 + pow(en/5.5, 2.5) + pow(en/60.0, 4.1)))
        + 0.05 / pow(1.0 + en/10.0, 2.0) + 0.01 * pow(en, 3.0) / (1.0 + pow(en/12.0, 6.0))); };
    auto qexc= [](const auto &en){ if(en>e_EXC_TH){ return 1e-20*(0.034 * pow(en-11.5, 1.1) * (1.0 + pow(en/15.0, 2.8)) / (1.0 + pow(en/23.0, 5.5)) 
        + 0.023 * (en-11.5) / pow(1.0 + en/80.0, 1.9)); } else { return 0.0;} };
    auto qion = [](const auto &en){ if(en>e_ION_TH){ return 1e-20*(970.0 * (en-15.8) / pow(70.0 + en, 2.0) +
                0.06 * pow(en-15.8, 2.0) * exp(-en/9)); } else {return 0.0;} };
    generate(e_energy.begin()+1, e_energy.end(), [i=1]()mutable { return CS_DE*(++i); });
    transform(e_energy.begin(),e_energy.end(), back_inserter(sigma[COLL::e_ELA]), qela);
    transform(e_energy.begin(),e_energy.end(), back_inserter(sigma[COLL::e_EXC]), qexc);
    transform(e_energy.begin(),e_energy.end(), back_inserter(sigma[COLL::e_ION]), qion);
    

    // Calculate ion cross sections(Data from Phelps A, JAP, 1994)

    auto qiso = [](const auto &e_lab){ return 2e-19 * pow(e_lab,-0.5) / (1.0 + e_lab) +
        3e-19 * e_lab / pow(1.0 + e_lab / 3.0, 2.0); };
    auto qmom= [](const auto &e_lab){ return 1.15e-18 * pow(e_lab,-0.1) * pow(1.0 + 0.015 / e_lab, 0.6); };
    auto qback = [&](const auto &x){ return (qmom(x)-qiso(x))/2.0; };
    generate(i_energy.begin()+1, i_energy.end(), [i=1]()mutable { return 2.0*CS_DE*(++i); });
    transform(i_energy.begin(),i_energy.end(), back_inserter(sigma[COLL::i_ISO]), qiso);
    transform(i_energy.begin(),i_energy.end(), back_inserter(sigma[COLL::i_BACK]), qback);

    // Non dimension energy array [] 

    const double sigma0 = n_0*vte0/omega_pe0;
    const double kTe0_inv = 1/kTe0;
    auto ndenergy = [&](auto& x) { x *= kTe0_inv; };
    auto ndcsection = [&](auto& x) { x *= sigma0; };

    for_each(e_energy.begin(), e_energy.end(), ndenergy);
    for_each(i_energy.begin(), i_energy.end(), ndenergy);
    for_each(sigma[COLL::e_ELA].begin(), sigma[COLL::e_ELA].end(), ndcsection);
    for_each(sigma[COLL::e_EXC].begin(), sigma[COLL::e_EXC].end(), ndcsection);
    for_each(sigma[COLL::e_ION].begin(), sigma[COLL::e_ION].end(), ndcsection);
    for_each(sigma[COLL::i_ISO].begin(), sigma[COLL::i_ISO].end(), ndcsection);
    for_each(sigma[COLL::i_BACK].begin(), sigma[COLL::i_BACK].end(), ndcsection);

    // Calculate e/ion total cross sections, neglect neutral collision
    const vector<string> name{"e", "p"};
    const vector<int> reaction{3, 2};
    const int react_num = reaction.size();
    ofstream f1("csection.in");
    f1 << "total " << reaction.size() << endl;
    for (int i = 0; i < react_num; ++i)
        f1 << "pairs " << name[i]
           << " reaction " << reaction[i] << endl;

    ofstream f2("e_Ar.dat");
    size_t arr_size = e_energy.size();
    double delta_E = e_energy[1] - e_energy[0];
    f2 << "bin " << arr_size << " de " << delta_E << endl;
    for (size_t i = 0; i < arr_size; ++i) 
        f2 << e_energy[i] << " " << sigma[COLL::e_ELA][i]
           << " " << sigma[COLL::e_EXC][i] << " " 
           << sigma[COLL::e_ION][i] << endl; 
    
    ofstream f3("i_Ar.dat");
    arr_size = i_energy.size();
    delta_E = i_energy[1] - i_energy[0];
    f3 << "bin " << arr_size << " de " << delta_E << endl;
    for (size_t i = 0; i < arr_size; ++i) 
        f2 << i_energy[i] << " " << sigma[COLL::i_ISO][i]
           << " " << sigma[COLL::i_BACK][i] << endl; 
    

    return 0;
}

void fprint(vector<double> vec, ofstream& ff){
    for(auto item : vec)
        ff << item << " ";
    ff << endl;
}

// cout << "Generate control.in:... " << endl;
//     const string space = "          ";
//     ofstream f1("control.in");
//     f1 << "dimension " << dimension << endl;
//     f1 << "seed      " << seed << endl;
//     f1 << "steps     " << NCYCLE*nstep_rf << endl;
//     f1 << "timestep  " << _dt_e << endl;
//     f1 << "screen    " << "100" << endl;
//     f1 << "model     " << "full" << endl;
//     f1 << "solver    " << "poisson lu" << endl;
//     f1 << "diag      " << "field/matlab 100000 100000 100000 &" << endl;
//     f1 << space << "potential number_density species all &" << endl;
//     f1 << space << "velocity species p components x y z &" << endl;
//     f1 << space << "temperature species e components x y z" << endl;
//     f1 << "read_restart false" << endl;
//     f1 << "write_restart " << "400000" << endl;
//     f1.close();
//     cout << "Generate mesh.in:... " << endl; 
//     ofstream f2("mesh.in");
//     f2 << "domain    " << nd_L[0] << " " << nd_L[1] << " "
//                        << nd_L[2] << " " << nd_L[3] << " "
//                        << nd_L[4] << " " << nd_L[5] << endl;
//     f2 << "num_cells " << Ng[0] << " " << Ng[1] << " " << Ng[2] << endl;
//     f2 << "tile      " << "8 8 1" << endl;
//     f2 << "field_bc type  n n s n p p &" << endl;
//     f2 << "         value 0 0 0 0 0 0" << endl;
//     f2 << "part_bc type   " << "v v r v p p" << endl;
//     f2 << "conductor rectangle type real &" << endl;
//     f2 << " position -2 1e-8 0. 10. 0. 0. &" << endl;
//     f2 << " potential fixed 0." << endl;
//     f2.close();
//     cout << "Generate particle.in:... " << endl;
//     ofstream f3("particle.in");
//     f3 << "species " << species[0] << " " << round(e_m/e_m*10)/10 << " -1.0 " << weight << endl;
//     f3 << "species " << species[1] << " " << round(i_m/e_m*10)/10 << " 1.0 " << weight << endl;
//     f3 << "Ground  " << species[2] << " " << round(i_m/e_m*10)/10 << " " << n_gas0 << endl;
//     f3.close();

// f4 << "e_energy ";
//     fprint(e_energy, f4);
//     f4 << "e_ela " ;
//     fprint(sigma[COLL::e_ELA], f4);
//     f4 << "e_exc " ;
//     fprint(sigma[COLL::e_EXC], f4);
//     f4 << "e_ion " ;
//     fprint(sigma[COLL::e_ION], f4);

//     f4 << "i_energy ";
//     fprint(i_energy, f4);
//     f4 << "i_iso " ;
//     fprint(sigma[COLL::i_ISO], f4);
//     f4 << "i_back " ;
//     fprint(sigma[COLL::i_BACK], f4);
//     f4.close(); 
