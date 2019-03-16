#include "Burgers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "mpi.h"
#include <stdio.h>

using namespace std;

Burgers::Burgers(Model *m_, MyMPI *myMPI_) : m(m_), myMPI(myMPI_) {
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
    for(unsigned int col=0; col < Nx; ++col ) {
        for(unsigned int row=0; row < Ny; ++row) {
            r = sqrt(pow(x(col),2)+pow(y(row),2));
            if (r <= r_thresh) {
                rb = 2*pow(1-r,4)*(4*r+1);
                adjustBounds(row,col);
            } else {
                rb = 0.0;
            }
            u[col*Ny+row] = rb;
            v[col*Ny+row] = rb;
        }
    }
}

void Burgers::printVelocityField() {
    serializeMatrix(u, "velocity_u.csv");
    serializeMatrix(v, "velocity_v.csv");
}

void Burgers::serializeMatrix(double *m, const char filename[]) {
    ofstream dataFile (filename, fstream::trunc);
    cout << "Writing velocity field data to file - " << filename << endl;
    for(unsigned int row=0; row<Ny; ++row) {
        for(unsigned int col=0; col<Nx; ++col) {
            dataFile << (col==0 ? ' ' : ',') << m[col*Ny+row];
        }
        dataFile << '\n';
    }
    dataFile.close();
}

void Burgers::integrateVelocityField() {
    splitDomain();
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
    const double p_i_j  = -2.0*c*(1.0/dx/dx + 1.0/dy/dy) - ax/dx - ay/dy + 1.0/dt;
    const double p_i_xj = c / dy / dy + ay / dy;
    const double p_i_jx = c / dy / dy;
    
    int i_j, ix_j, xi_j, i_jx, i_xj;
    double tempu, tempv;
    
	typedef std::chrono::high_resolution_clock hrc;
    hrc::time_point loop_start = hrc::now();
    
    auto* un = new double[Nx*Ny];
    auto* vn = new double[Nx*Ny];
    for (int t = 1; t <= Nt; ++t) {
        #pragma omp parallel for private(i_j, ix_j, xi_j, i_jx, i_xj, tempu, tempv) ordered
        for(unsigned int col = lbound; col <= rbound; ++col) {
            for(unsigned int row = tbound; row <= bbound; ++row) {
                i_j = col*Ny+row;
                i_jx= i_j + 1;
                i_xj= i_j - 1;
                ix_j= i_j + Ny;
                xi_j= i_j - Ny;
                
                tempu = p_xi_j * u[xi_j] + p_ix_j * u[ix_j] + p_i_j * u[i_j] + p_i_xj * u[i_xj] + p_i_jx * u[i_jx];
                tempu+= -1.0 * b / dx   * u[i_j] * (u[i_j] - u[xi_j]) - b / dy * v[i_j] * (u[i_j] - u[i_xj]);
    
                tempv = p_xi_j * v[xi_j] + p_ix_j * v[ix_j] + p_i_j * v[i_j] + p_i_xj * v[i_xj] + p_i_jx * v[i_jx];
                tempv+= (-1.0 * b / dx) * u[i_j] * (v[i_j] - v[xi_j]) - (b / dy ) * v[i_j] * (v[i_j] - v[i_xj]);
    
                un[i_j] = dt * tempu;
                vn[i_j] = dt * tempv;
                if (fabs(u[i_j]) > 1e-8 || fabs(v[i_j]) > 1e-8) adjustBounds(row, col);
            }
        }
        swap(un, u);
        swap(vn, v);
        
//        rollbackBounds();

		if (verbose && t%100 == 0)
			cout << "Time Step: " << t << " of " << Nt
				 << "\tRunning time: " << chrono::duration_cast<chrono::milliseconds>(hrc::now() - loop_start).count()
				 << "ms"
				 << "\tl: " << lbound << " r: " << rbound << " t: " << tbound << " b: " << bbound
				 << endl;
    }
    delete[] vn;
    delete[] un;
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

void Burgers::adjustBounds(unsigned int row, unsigned int col) {
    if(col <  Nx-2 && col > 1) {
        if (lbound >= col) lbound = col - 1;
        if (rbound <= col) rbound = col + 1;
    }
    if(row <  Ny-2 && row > 1) {
        if (tbound >= row) tbound = row - 1;
        if (bbound <= row) bbound = row + 1;
    }
}

void Burgers::rollbackBounds() {
    if(lbound > 1) lbound--;
    if(rbound < Nx-2) rbound++;
    if(tbound > 1) tbound--;
    if(bbound < Ny-2) bbound++;
}

void Burgers::splitDomain() {
    rankx = world_rank / Py;
    ranky = world_rank - rankx * Py;
    auto dvx = div((int)Nx-1,Px);
    for(int i=0; i < rankx; ++i) {
        worldx+= dvx.quot + ( i < dvx.rem? 1 : 0) + 1;
    }
    locNx = dvx.quot + ( rankx < dvx.rem? 1 : 0) + 1;
    
    auto dvy = div((int)Ny-1,Py);
    for(int i=0; i < ranky; ++i) {
        worldy+= dvy.quot + (i < dvy.rem ? 1 : 0) + 1;
    }
    locNy = dvy.quot + ( ranky < dvy.rem? 1 : 0) + 1;
    
    worldRef = worldx*Ny + worldy;
    cout    << "P" << world_rank
            << " worldref: " << worldRef
            << " worldx: " << worldx
            << " worldy: " << worldy
            << " rankx: " << rankx
            << " ranky: " << ranky
            << " locNx: " << locNx
            << " locNy: " << locNy << endl;
}
