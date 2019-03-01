#ifndef CLASS_MODEL
#define CLASS_MODEL

#include <iostream>
#include <vector>
#include <string>
using namespace std;

//#include "mpi.h"

class Model {
	public:
		Model(int argc, char* argv[]);
//        ~Model();

        void printParameters();

        // Getters
        bool            isValid()   const { return valid; }
		bool            isVerbose() const { return verbose; }
        bool            isHelp()    const { return help; }
        double          getX0()     const { return x0; }
        double          getY0()     const { return y0; }
        double          getLx()     const { return Lx; }
        double          getLy()     const { return Ly; }
        double          getT()      const { return T; }
        unsigned int    getNx()     const { return Nx; }
        unsigned int    getNy()     const { return Ny; }
        unsigned int    getNt()     const { return Nt; }
        double          getDx()     const { return dx; }
        double          getDy()     const { return dy; }
        double          getDt()     const { return dt; }
        double          getAx()     const { return ax; }
        double          getAy()     const { return ay; }
        double          getB()      const { return b; }
        double          getC()      const { return c; }

	private:
        bool parseArguments(int argc, char* argv[]);
        void printHelp();
        void validateParameters();

        // Meta Parameters
	    bool            verbose;
        bool            help;
        bool            valid = false;
        string          fname;

	    // Numerics
	    double          x0   =   0;
	    double          y0   =   0;
	    double          Lx   =   10;
	    double          Ly   =   10;
	    double          T    =   1;
	    unsigned int    Nx   =   20;
	    unsigned int    Ny   =   20;
	    unsigned int    Nt   =   40;
	    double          dx   =   Lx/(Nx-1);
	    double          dy   =   Ly/(Ny-1);
	    double          dt   =   T/Nt;

	    // Physics
	    double          ax   =   0;
	    double          ay   =   0;
	    double          b    =   0;
	    double          c    =   0;
};

Model::Model(int argc, char* argv[]) {
    if(!Model::parseArguments(argc, argv)) return;
    if (Model::verbose) Model::printParameters();
};

bool Model::parseArguments(int argc, char* argv[]) {
    for(int i=1; i < argc; ++i ) {
        string arg = string(argv[i]);
        if (arg == "-h" || arg == "--help") {
            printHelp();
            return false;
        }
        else if (arg == "-v" || arg == "--verbose") {
            Model::verbose = true;
        }
        else if (arg == "-p" || arg == "--physics") {
            try {
                double* physical_parameters[] = {&(Model::ax), &(Model::ay), &(Model::b), &(Model::c)};
                for (int j = 1; j < 5; ++j) {
                    i++;
                    *physical_parameters[j-1] = stod(string(argv[i]));
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
                unsigned int* physical_parameters[] = {&(Model::Nx), &(Model::Ny), &(Model::Nt)};
                for (int j = 1; j < 4; ++j) {
                    i++;
                    *physical_parameters[j-1] = stod(string(argv[i]));
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

void Model::printHelp() {
    help = true;
    cerr << "Usage: " << fname << "ax ay b c <option(s)>" << endl
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

#endif