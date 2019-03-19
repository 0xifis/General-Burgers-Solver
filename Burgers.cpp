#include "Burgers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "mpi.h"
#include <stdio.h>
#include <algorithm>

using namespace std;

Burgers::Burgers(Model *m_, MyMPI *myMPI_) : m(m_), myMPI(myMPI_) {
    u = new double[Nx*Ny];
    v = new double[Nx*Ny];
    splitDomain();
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
    for(unsigned int col=worldx; col < worldx+locNx; ++col ) {
        for(unsigned int row=worldy; row < worldy+locNy; ++row) {
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
    cout << "\tl: " << lbound << " r: " << rbound << " t: " << tbound << " b: " << bbound
               << endl;
}

void Burgers::printVelocityField() {
    serializeMatrix(u, "velocity_u.csv");
    serializeMatrix(v, "velocity_v.csv");
}

void Burgers::serializeMatrix(double *m, string filename) {
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
    unsigned int llbound, lrbound, ltbound, lbbound;
    double tempu, tempv;
    
	typedef std::chrono::high_resolution_clock hrc;
    hrc::time_point loop_start = hrc::now();
    
    auto* un = new double[Nx*Ny];
    auto* vn = new double[Nx*Ny];
    for (int t = 1; t <= Nt; ++t) {
        llbound = max(lbound, worldx+1);
        lrbound = min(rbound, worldx+locNx-2);
        ltbound = max(tbound, worldy+1);
        lbbound = min(bbound, worldy+locNy-2);
        #pragma omp parallel for private(i_j, ix_j, xi_j, i_jx, i_xj, tempu, tempv) ordered
        for(unsigned int col = llbound; col <= lrbound; ++col) {
            #pragma omp simd
            for(unsigned int row = ltbound; row <= lbbound; ++row) {
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
                if (fabs(un[i_j]) > 1e-8 || fabs(vn[i_j]) > 1e-8) adjustBounds(row, col);
            }
        }
        swap(un, u);
        swap(vn, v);
        exchangePadding();
//        rollbackBounds();
        exchangeBounds();



		if (verbose && t%500 == 0)
			cout << "Time Step: " << t << " of " << Nt
				 << "\tRunning time: " << chrono::duration_cast<chrono::milliseconds>(hrc::now() - loop_start).count()
				 << "ms"
				 << "\tl: " << lbound << " r: " << rbound << " t: " << tbound << " b: " << bbound
				 << endl;
    }
    delete[] vn;
    delete[] un;
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
    auto dvx = div((int)Nx-2,Px);
    for(int i=0; i < rankx; ++i) {
        worldx+= dvx.quot + ( i < dvx.rem? 1 : 0);
    }
    locNx = dvx.quot + ( rankx < dvx.rem? 1 : 0) + 2;
    
    auto dvy = div((int)Ny-2,Py);
    for(int i=0; i < ranky; ++i) {
        worldy+= dvy.quot + (i < dvy.rem ? 1 : 0);
    }
    locNy = dvy.quot + ( ranky < dvy.rem? 1 : 0) + 2;
    
    worldRef = worldx*Ny + worldy;
    if(m->isVerbose()) cout << "P" << world_rank
                            << " worldref: " << worldRef
                            << " worldx: " << worldx
                            << " worldy: " << worldy
                            << " rankx: " << rankx
                            << " ranky: " << ranky
                            << " locNx: " << locNx
                            << " locNy: " << locNy << endl;
}

double Burgers::getFieldEnergy() {
    return worldEnergy;
}

void Burgers::calculateFieldEnergy() {
    const double dx = m->getDx();
    const double dy = m->getDy();
    #pragma omp parallel for reduction(+:energy)
    for(unsigned int col=worldx+1; col <= worldx+locNx-2; ++col) {
        for (unsigned int row = worldy+1; row <= worldy+locNy-2; ++row) {
            energy += u[col*Ny+row] * u[col*Ny+row] + v[col*Ny+row] * v[col*Ny+row];
        }
    }
    MPI_Reduce(&energy, &worldEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    worldEnergy *= dx*dy*0.5;
}

void Burgers::exchangePadding() {
    if (Px > 1) sendAndReceiveCols();
    if (Py > 1) sendAndReceiveRows();
}

void Burgers::sendAndReceiveCols() {
//    MPI_Send( void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
    MPI_Datatype inner_col;
    MPI_Type_vector(1,locNy-2,locNy-2,MPI_DOUBLE,&inner_col);
    MPI_Type_commit(&inner_col);
    
    int right_inner_col = worldRef+(locNx-2)*Ny+1;
    int right_padding_col = worldRef+(locNx-1)*Ny+1;
    
    int left_inner_col = worldRef+Ny+1;
    int left_padding_col = worldRef+1;
    
    if(rankx > 0 && rankx < Px-1) { // Center domains; send left & right
        auto* reqs = new MPI_Request[8];
        MPI_Isend(&u[right_inner_col], 1, inner_col, getRank(rankx+1,ranky), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[right_inner_col], 1, inner_col, getRank(rankx+1,ranky), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[right_padding_col], 1, inner_col, getRank(rankx+1,ranky), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[right_padding_col], 1, inner_col, getRank(rankx+1,ranky), 1, MPI_COMM_WORLD, &reqs[3]);
    
        MPI_Isend(&u[left_inner_col], 1, inner_col, getRank(rankx-1,ranky), 0, MPI_COMM_WORLD, &reqs[4]);
        MPI_Isend(&v[left_inner_col], 1, inner_col, getRank(rankx-1,ranky), 1, MPI_COMM_WORLD, &reqs[5]);
        MPI_Irecv(&u[left_padding_col], 1, inner_col, getRank(rankx-1,ranky), 0, MPI_COMM_WORLD, &reqs[6]);
        MPI_Irecv(&v[left_padding_col], 1, inner_col, getRank(rankx-1,ranky), 1, MPI_COMM_WORLD, &reqs[7]);
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    } else if (rankx == 0) { //Left domain; send & recv right col,
//        cout << "P" << world_rank << ": sending right col to " << worldRef+(locNx-2)*Ny+1 << endl;
//        cout << "P" << world_rank << ": sending left col to " << worldRef+Ny+1 << endl;
        auto* reqs = new MPI_Request[4];
        MPI_Isend(&u[right_inner_col], 1, inner_col, getRank(rankx+1,ranky), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[right_inner_col], 1, inner_col, getRank(rankx+1,ranky), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[right_padding_col], 1, inner_col, getRank(rankx+1,ranky), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[right_padding_col], 1, inner_col, getRank(rankx+1,ranky), 1, MPI_COMM_WORLD, &reqs[3]);
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    } else { //Right domain; send & recv left col,
        auto* reqs = new MPI_Request[4];
        MPI_Isend(&u[left_inner_col], 1, inner_col, getRank(rankx-1,ranky), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[left_inner_col], 1, inner_col, getRank(rankx-1,ranky), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[left_padding_col], 1, inner_col, getRank(rankx-1,ranky), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[left_padding_col], 1, inner_col, getRank(rankx-1,ranky), 1, MPI_COMM_WORLD, &reqs[3]);
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    }
//    TODO: pull up the waitall
}

void Burgers::sendAndReceiveRows() {
//    MPI_Send( void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
    MPI_Datatype inner_row;
    MPI_Type_vector(locNx-2,1,Ny,MPI_DOUBLE,&inner_row);
    MPI_Type_commit(&inner_row);
    
    int bot_inner_row = worldRef+Ny+locNy-2;
    int bot_padding_row= worldRef+Ny+locNy-1;
    
    int top_inner_row = worldRef+Ny+1;
    int top_padding_row = worldRef+Ny;
    
    if(ranky>0 && ranky < Py-1) { // Center domains; send left & right
        auto* reqs = new MPI_Request[8];
        MPI_Isend(&u[bot_inner_row], 1, inner_row, getRank(rankx,ranky+1), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[bot_inner_row], 1, inner_row, getRank(rankx,ranky+1), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[bot_padding_row], 1, inner_row, getRank(rankx,ranky+1), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[bot_padding_row], 1, inner_row, getRank(rankx,ranky+1), 1, MPI_COMM_WORLD, &reqs[3]);
        
        MPI_Isend(&u[top_inner_row], 1, inner_row, getRank(rankx,ranky-1), 0, MPI_COMM_WORLD, &reqs[4]);
        MPI_Isend(&v[top_inner_row], 1, inner_row, getRank(rankx,ranky-1), 1, MPI_COMM_WORLD, &reqs[5]);
        MPI_Irecv(&u[top_padding_row], 1, inner_row, getRank(rankx,ranky-1), 0, MPI_COMM_WORLD, &reqs[6]);
        MPI_Irecv(&v[top_padding_row], 1, inner_row, getRank(rankx,ranky-1), 1, MPI_COMM_WORLD, &reqs[7]);
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    } else if (ranky==0) { //Left domain; send & recv right col,
//        cout << "P" << world_rank << ": sending bot row to " << bot_inner_row << endl;
//        cout << "P" << world_rank << ": waiting for " << bot_padding_row << endl;
        auto* reqs = new MPI_Request[4];
        MPI_Isend(&u[bot_inner_row], 1, inner_row, getRank(rankx,ranky+1), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[bot_inner_row], 1, inner_row, getRank(rankx,ranky+1), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[bot_padding_row], 1, inner_row, getRank(rankx,ranky+1), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[bot_padding_row], 1, inner_row, getRank(rankx,ranky+1), 1, MPI_COMM_WORLD, &reqs[3]);
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    } else { //Right domain; send & recv left col,
        auto* reqs = new MPI_Request[4];
//        cout << "P" << world_rank << ": sending top row to " << top_inner_row << endl;
//        cout << "P" << world_rank << ": waiting for " << top_padding_row << endl;
        MPI_Isend(&u[top_inner_row], 1, inner_row, getRank(rankx,ranky-1), 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&v[top_inner_row], 1, inner_row, getRank(rankx,ranky-1), 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(&u[top_padding_row], 1, inner_row, getRank(rankx,ranky-1), 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&v[top_padding_row], 1, inner_row, getRank(rankx,ranky-1), 1, MPI_COMM_WORLD, &reqs[3]);
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
        delete[] reqs;
    }
}

int Burgers::getRank(int rankx, int ranky) {
    return rankx*Py+ranky;
}

void Burgers::exchangeBounds() {
    MPI_Allreduce(MPI_IN_PLACE, &lbound, 1, MPI_UNSIGNED, MPI_MIN , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tbound, 1, MPI_UNSIGNED, MPI_MIN , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &rbound, 1, MPI_UNSIGNED, MPI_MAX , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &bbound, 1, MPI_UNSIGNED, MPI_MAX , MPI_COMM_WORLD);
}
