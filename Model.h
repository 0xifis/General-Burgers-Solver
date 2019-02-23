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

        void PrintParameters() {};

        bool IsValid();

        void PrintHelp();

        // Getters
		bool   IsVerbose() const { return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetLx()     const { return Lx; }
        double GetLy()     const { return Ly; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }

        // Add any other getters here...

	private:
        void ParseParameters(vector <string> argv);

        void ValidateParameters();

	    bool verbose;
        bool help;
        bool valid;
        string fname;

	    // Numerics
	    double x0;
	    double y0;
	    double Lx;
	    double Ly;
	    double T;
	    int    Nx;
	    int    Ny;
	    int    Nt;
	    double dx;
	    double dy;
	    double dt;

	    // Physics
	    double ax;
	    double ay;
	    double b;
	    double c;

        // Add any additional parameters here...
};

Model::Model(int argc, char* argv[]) {
    vector <string> parameters_physics;
    for(int i=1; i < argc; ++i ) {
        string arg = argv[i];
        if (arg == "-h" || arg =="--help") {
            PrintHelp();
        } else if (arg == "-v" || arg == "-verbose") {
            verbose = true;
        } else {
            parameters_physics.push_back(arg);
        }
    }
    if (!help) {
        ParseParameters(parameters_physics);
        Lx = 10;
        Ly = 10;
        T = 0.0;

        if (verbose) {
            PrintParameters();
        }
    }
};

void Model::ParseParameters(vector <string> argv) {
    int argc = argv.size();
    if (argc != 4) {
        cout << "Incorrect number of arguments (" << argc << " instead of 4)" << endl;
        PrintHelp();
        valid = false;
    } else {
        ax = stod(argv[0]);
        ay = stod(argv[1]);
        b = stod(argv[2]);
        c = stod(argv[3]);
    }
}

void Model::PrintHelp() {
    help = true;
    cerr << "Usage: " << fname << "ax ay b c <option(s)>" << endl
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
         << "\t-v,--verbose\t\tVerbose Output\n"
         << endl;
}

#endif