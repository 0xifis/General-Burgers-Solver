#ifndef BURGERS_H
#define BURGERS_H

#include "Model.h"
#include "MyMPI.h"

class Burgers {
    public:
        explicit Burgers(Model *m_, MyMPI *myMPI_);
        void initializeVelocityField();
        void printVelocityField();
        void integrateVelocityField();
        double fieldEnergy();

    private:
        void serializeMatrix(double* m, const char filename[]);
        void adjustBounds(unsigned int row, unsigned int col);
        void splitDomain();
        Model* m;
        MyMPI* myMPI;
        unsigned int world_rank = myMPI->rank();
        double* u;
        double* v;
        unsigned int Nx = m->getNx();
        unsigned int Ny = m->getNy();
        unsigned int Px = m->getPx();
        unsigned int Py = m->getPy();
        unsigned int locNx, locNy;
        double x(int col), y(int row);
        double r_thresh = 1;
        unsigned int lbound = Ny;
        unsigned int rbound = 0;
        unsigned int tbound = Nx;
        unsigned int bbound = 0;
        void rollbackBounds();
        bool verbose = m->isVerbose();
};

#endif