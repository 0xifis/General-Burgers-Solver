#ifndef BURGERS_H
#define BURGERS_H

#include "Model.h"

class Burgers {
    public:
        explicit Burgers(Model *m_);
        void initializeVelocityField();
        void printVelocityField();
        void integrateVelocityField();
        double fieldEnergy();

    private:
        Model* m;
        double* u;
        double* v;
        unsigned int Nx, Ny;
        double x(int col), y(int row);
        void serializeMatrix(double* m, ofstream* dataFile);
        double r_thresh = 1;
        unsigned int lbound = Nx;
        unsigned int rbound = 0;
        unsigned int tbound = Ny;
        unsigned int bbound = 0;
        void adjustBounds(unsigned int Nx, unsigned int Ny);
    
    void rollbackBounds();
};

#endif