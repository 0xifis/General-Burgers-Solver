#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <string>
using namespace std;

//#include "mpi.h"

class Model {
	public:
		Model(int argc, char* argv[]) {
            if(!Model::parseArguments(argc, argv)) return;
            Model::validateParameters();
		};
//        ~Model(); //TODO

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

	    // Geometric
	    double          x0   =   0;
	    double          y0   =   0;
	    double          Lx   =   10;
	    double          Ly   =   10;

	    double          T    =   1;

	    // Numeric
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

#endif