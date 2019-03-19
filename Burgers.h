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
        double getFieldEnergy();
        void calculateFieldEnergy();

    private:
        void serializeMatrix(double* m, string filename);
        void adjustBounds(unsigned int row, unsigned int col);
        void splitDomain();
        void exchangePadding();
        void sendAndReceiveCols();
        void sendAndReceiveRows();
        int getRank(int rankx, int ranky);
        Model* m;
        double* u;
        double* v;
        unsigned int Nx = m->getNx();
        unsigned int Ny = m->getNy();
        double x(int col), y(int row);
        double r_thresh = 1;
        unsigned int lbound = Ny;
        unsigned int rbound = 0;
        unsigned int tbound = Nx;
        unsigned int bbound = 0;
        void rollbackBounds();
        bool verbose = m->isVerbose();
    
        double worldEnergy = 0.0;
        double energy = 0.0;
    
        MyMPI* myMPI;
        unsigned int world_rank = myMPI->rank();
        unsigned int rankx;
        unsigned int ranky;
        unsigned int Px = m->getPx();
        unsigned int Py = m->getPy();
        unsigned int locNx, locNy;
        unsigned int worldx = 0;
        unsigned int worldy = 0;
        unsigned int worldRef = 0;
};

#endif