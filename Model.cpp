#include "Model.h"

bool Model::parseArguments(int argc, char* argv[]) {
    for(int i=1; i < argc; ++i ) {
        string arg = string(argv[i]);
        if (arg == "-h" || arg == "--help") {
            printHelp();
        }
        else if (arg == "-v" || arg == "--verbose") {
            Model::verbose = true;
        }
        else if (arg == "-p" || arg == "--physics") {
            try {
                double* physical_parameters[] = {&(Model::ax), &(Model::ay), &(Model::b), &(Model::c)};
                for (auto &physical_parameter : physical_parameters) {
                    i++;
                    *physical_parameter = stod(string(argv[i]));
                }
            }
            catch(const exception& e) {
                cerr << "Error: " << e.what() << endl;
                cerr << "-p/--physics option requires four arguments of type double (ax, ay, b, c)." << endl;
                return false;
            }
        }
        else if (arg == "-t" || arg == "--time") {
            try {
                i++;
                Model::T = stod(string(argv[i]));
            }
            catch(const exception& e) {
                cerr << "Error: " << e.what() << endl;
                cerr << "-t/--time option requires one argument of type double (T in seconds)." << endl;
                return false;
            }
        }
        else if (arg == "-g" || arg == "--geometry") {
            try {
                i++;
                Model::Lx = Model::Ly = stod(string(argv[i]));
            }
            catch(const exception& e) {
                cerr << "Error: " << e.what() << endl;
                cerr << "-g/--geometry option requires one argument of type double (L in meters)." << endl;
                return false;
            }
        }
        else if (arg == "-n" || arg == "--numeric") {
            try {
                unsigned int* numerical_parameters[] = {&(Model::Nx), &(Model::Ny), &(Model::Nt)};
                for (auto &numerical_parameter : numerical_parameters) {
                    i++;
                    *numerical_parameter = stoi(string(argv[i]));
                }
            } catch(const exception& e) {
                cerr << "Error: " << e.what() << endl;
                cerr << "-n/--numeric option requires three positive arguments of type int (Nx, Ny, Nt)." << endl;
                return false;
            }
        }
    }
    return true;
}

void Model::validateParameters() {
    Model::valid = true; // TODO
}

void Model::printHelp() {
    Model::help = true;
    cerr << "Usage: " << fname << "ax ay b c <option(s)>\n"
         << "\t-p,--physics <ax> <ay> <b> <c> \tSet physics parameters (default: ax=10, ay=0, b=0, c=0)\n"
         << "\t-g,--geometry <Lx> <Ly> <x0> <y0> \tSet length and initial position of velocity field considered (default: Lx=10, Ly=10, x0=0, y0=0)\n"
         << "\t-t,--time <end time in seconds> \tSet simulation end time (default: 1s)\n"
         << "\t-n,--numeric <Nx> <Ny> <Nt> \tSet simulation numeric parameters (default: Nx=20, Ny=20, Nt=20)\n"
         << "\t-h,--help\t\t\t\tShow this help message\n"
         << "\t-v,--verbose\t\t\t\tRun in verbose mode\n"
         << endl;
}

void Model::printParameters() {
    cout << "The overall model is " << (Model::valid ? "valid." : "not valid.") << "\n\n"
         << "Physics Parameters:\n"
         << "\tax\t" << ":" << "\t" << Model::ax << "\n"
         << "\tay\t" << ":" << "\t" << Model::ay << "\n"
         << "\tb\t" << ":" << "\t" << Model::b << "\n"
         << "\tc\t" << ":" << "\t" << Model::c << "\n"
         << "Geometric Parameters:\n"
         << "\tLx\t" << ":" << "\t" << Model::Lx << "\n"
         << "\tLy\t" << ":" << "\t" << Model::Ly << "\n"
         << "\tx0\t" << ":" << "\t" << Model::x0 << "\n"
         << "\ty0\t" << ":" << "\t" << Model::y0 << "\n"
         << "Numerical Parameters:\n"
         << "\tNx\t" << ":" << "\t" << Model::Nx << "\n"
         << "\tNy\t" << ":" << "\t" << Model::Ny << "\n"
         << "\tNt\t" << ":" << "\t" << Model::Nt << "\n"
         << "\tT\t"  << ":" << "\t" << Model::T << "\n"
         << "\tdx\t" << ":" << "\t" << Model::dx << "\n"
         << "\tdy\t" << ":" << "\t" << Model::dy << "\n"
         << "\tdt\t" << ":" << "\t" << Model::dt << "\n"
         << endl;
}
