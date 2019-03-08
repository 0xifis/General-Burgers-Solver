#include "Burgers.h"
#include <blaze/Math.h>
#include <fstream>

Burgers::Burgers(Model *m_) : m(m_) {
    Nx = m->getNx();
    Ny = m->getNy();
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
    // use a submatrix here
    double r, rb;
    for( int row=0; row<Nx; ++row ) {
        for( int col=0; row<Ny; ++col) {
            r = sqrt(pow(x(col),2)+pow(y(row),2));
            rb = (r < r_thresh ? 2*pow(1-r,4)*(4*r+1) : 0.0);
            u[row*Ny+col] = rb;
            v[row*Ny+col] = rb;
        }
    }
}

void Burgers::printVelocityField() {
    const auto dotStep = static_cast<const size_t> (floor(u.rows() / 10.0));
    const char filename[] = "velocity_u.csv";
    cout << "Writing velocity field data to file - " << filename << endl;
    ofstream dataFile (filename, fstream::trunc);
    for( size_t col=0; col<Nx; ++col) {
        for( size_t row=0; row<Nx; ++row) {
            dataFile << (col==0 ? ' ' : ',') << u[row*Ny+col];
        }
        dataFile << '\n';
        if (col%dotStep == 0) cout << '.';
    }
    dataFile.close();
    cout << "\nDone :)" << endl;
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"
void Burgers::integrateVelocityField() {
    const auto nsx = Nx;
    const auto nsy = Ny;
    const double Nt = m->getNt();
    const double dx = m->getDx();
    const double dy = m->getDy();
    const double dt = m->getDt();
    const double ax  = m->getAx();
    const double ay  = m->getAy();
    const double b  = m->getB();
    const double c  = m->getC();

//    auto u_ix_j = blaze::submatrix(u, 1UL, 2UL, nsx, nsy);
//    auto u_xi_j = blaze::submatrix(u, 1UL, 0UL, nsx, nsy);
//    auto u_i_j  = blaze::submatrix(u, 1UL, 1UL, nsx, nsy);
//    auto u_i_xj = blaze::submatrix(u, 0UL, 1UL, nsx, nsy);
//    auto u_i_jx = blaze::submatrix(u, 2UL, 1UL, nsx, nsy);
//
//    auto v_ix_j = blaze::submatrix(v, 1UL, 2UL, nsx, nsy);
//    auto v_xi_j = blaze::submatrix(v, 1UL, 0UL, nsx, nsy);
//    auto v_i_j  = blaze::submatrix(v, 1UL, 1UL, nsx, nsy);
//    auto v_i_xj = blaze::submatrix(v, 0UL, 1UL, nsx, nsy);
//    auto v_i_jx = blaze::submatrix(v, 2UL, 1UL, nsx, nsy);
    
    const double p_xi_j = c / dx / dx + ax / dx;
    const double p_ix_j = c / dx / dx;
    const double p_i_j  = -2.0*c*(1/dx/dx + 1/dy/dy) - ax/dx - ay/dy + 1/dt;
    const double p_i_xj = c / dy / dy + ay / dy;
    const double p_i_jx = c / dy / dy;
    
    int i_j, ix_j, xi_j, i_jx, i_xj;
    double tempu, tempv;
    
//    blaze::DynamicMatrix<double> U34 (nsx, nsy), U23 (nsx, nsy), V34 (nsx, nsy), V23 (nsx, nsy);
	
	typedef std::chrono::high_resolution_clock hrc;
	typedef std::chrono::milliseconds ms;
    hrc::time_point loop_start = hrc::now();
	
    for (int t = 1; t <= Nt; ++t) {
        for(int row = 1; row < Nx-1; ++row) {
            for(int col = 1; col < Ny-1; ++col) {
                i_j = row*Ny+col;
                i_jx= i_j + 1;
                i_xj= i_j - 1;
                ix_j= i_j + Ny;
                xi_j= i_j - Ny;
                
                tempu = p_xi_j * u[xi_j] + p_ix_j * u[ix_j] + p_i_j * u[i_j] + p_i_xj * u[i_xj] + p_i_jx * u[i_xj];
                tempu+= -1.0 * b / dx * u[i_j] * (u[i_j]- u[xi_j]) - b / dy * v[i_j] * (u[i_j] - u[i_xj]);
    
                tempv = p_xi_j * v[xi_j] + p_ix_j * v[ix_j] + p_i_j * v[i_j] + p_i_xj * v[i_xj] + p_i_jx * v[i_jx];
                tempv+= (-1.0 * b / dx) * u[i_j] * (v[i_j] - v[xi_j]) - (b / dy * v[i_j]) * (v[i_j] - v[i_xj]);
    
                u[i_j] = dt * tempu;
                v[i_j] = dt * tempv;
            }
        }
		
		if (t%10 == 0)
			cout << "Time Step: " << t << " of " << Nt
				 << "\tRunning time: " << chrono::duration_cast<chrono::milliseconds>(hrc::now() - loop_start).count()
				 << "ms"
				 << endl;
    }
    cout << "\n\nDone\n\n";
}
#pragma clang diagnostic pop

double Burgers::fieldEnergy() {
    return 0.0;
//    return 0.5*sum(blaze::sqrt(u*trans(u) + v*trans(v)))*m->getDx()*m->getDy();
}