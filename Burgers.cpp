#include "Burgers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>


using namespace std;

Burgers::Burgers(Model *m_) : m(m_), Nx(m->getNx()), Ny(m->getNy()) {
    u = new double[Nx*Ny];
    v = new double[Nx*Ny];
}

double Burgers::x(int col) {
    const double dx = m -> getDx();
    const double Lx = m -> getLx();
    return (dx * col - Lx / 2);
}

double Burgers::y(int row) {
    const double dy = m -> getDy();
    const double Ly = m -> getLy();
    return (dy * row - Ly / 2);
}

void Burgers::initializeVelocityField() {
    
    double r, rb;
    for(unsigned int row=0; row < Nx; ++row ) {
        for(unsigned int col=0; col < Ny; ++col) {
            r = sqrt(pow(x(col),2)+pow(y(row),2));
            if (r < r_thresh) {
                rb = 2*pow(1-r,4)*(4*r+1);
                adjustBounds(col,row);
            } else {
                rb = 0.0;
            }
            u[row*Ny+col] = rb;
            v[row*Ny+col] = rb;
        }
    }
}

void Burgers::printVelocityField() {
    const char filename[] = "velocity_u.csv";
    cout << "Writing velocity field data to file - " << filename << endl;
    ofstream dataFile (filename, fstream::trunc);
    serializeMatrix(u, &dataFile);
    dataFile.close();
    cout << "\nDone :)" << endl;
}

void Burgers::serializeMatrix(double *m, ofstream* dataFile) {
    for( size_t col=0; col<Nx; ++col) {
        for( size_t row=0; row<Nx; ++row) {
            *dataFile << (col==0 ? ' ' : ',') << m[row*Ny+col];
        }
        *dataFile << '\n';
    }
}

void Burgers::integrateVelocityField() {
    const double Nt = m->getNt();
    const double dx = m->getDx();
    const double dy = m->getDy();
    const double dt = m->getDt();
    const double ax  = m->getAx();
    const double ay  = m->getAy();
    const double b  = m->getB();
    const double c  = m->getC();
    
    const double p_xi_j = c / dx / dx + ax / dx;
    const double p_ix_j = c / dx / dx;
    const double p_i_j  = -2.0*c*(1/dx/dx + 1/dy/dy) - ax/dx - ay/dy + 1/dt;
    const double p_i_xj = c / dy / dy + ay / dy;
    const double p_i_jx = c / dy / dy;
    
    int i_j, ix_j, xi_j, i_jx, i_xj;
    double tempu, tempv;
    
	typedef std::chrono::high_resolution_clock hrc;
    hrc::time_point loop_start = hrc::now();
	
    for (int t = 1; t <= Nt; ++t) {
        for(unsigned int row = lbound; row <= rbound; ++row) {
            for(unsigned int col = tbound; col <= bbound; ++col) {
                i_j = row*Ny+col;
                i_jx= i_j + 1;
                i_xj= i_j - 1;
                ix_j= i_j + Ny;
                xi_j= i_j - Ny;
                
                tempu = p_xi_j * u[xi_j] + p_ix_j * u[ix_j] + p_i_j * u[i_j] + p_i_xj * u[i_xj] + p_i_jx * u[i_jx];
                tempu+= -1.0 * b / dx * u[i_j] * (u[i_j]- u[xi_j]) - b / dy * v[i_j] * (u[i_j] - u[i_xj]);
    
                tempv = p_xi_j * v[xi_j] + p_ix_j * v[ix_j] + p_i_j * v[i_j] + p_i_xj * v[i_xj] + p_i_jx * v[i_jx];
                tempv+= (-1.0 * b / dx) * u[i_j] * (v[i_j] - v[xi_j]) - (b / dy * v[i_j]) * (v[i_j] - v[i_xj]);
    
                u[i_j] = dt * tempu;
                v[i_j] = dt * tempv;
            }
        }
        rollbackBounds();
		
		if (t%10 == 0)
			cout << "Time Step: " << t << " of " << Nt
				 << "\tRunning time: " << chrono::duration_cast<chrono::milliseconds>(hrc::now() - loop_start).count()
				 << "ms"
				 << endl;
    }
    cout << "\n\nDone\n\n";
}

double Burgers::fieldEnergy() {
    double energy = 0.0;
    const double dx = m->getDx();
    const double dy = m->getDy();
    for (unsigned int i = 0; i < Nx * Ny; ++i) {
        energy += 0.5 * (u[i] * u[i] + v[i] * v[i]) * dx * dy;
    }
    return energy;
}

void Burgers::adjustBounds(unsigned int col, unsigned int row) {
    if(lbound > col) lbound = col;
    if(rbound < col) lbound = col;
    if(tbound > row) tbound = row;
    if(bbound < row) bbound = row;
}

void Burgers::rollbackBounds() {
    if(lbound > 1) lbound--;
    if(rbound < Nx-1) rbound++;
    if(tbound > 0) tbound--;
    if(bbound < Ny-1) bbound++;
}
