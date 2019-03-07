#include "Burgers.h"
#include <blaze/Math.h>
#include <fstream>

Burgers::Burgers(Model *m_) : m(m_) {
    const unsigned int Nx = m->getNx();
    const unsigned int Ny = m->getNy();
    Burgers::u.resize( Nx, Ny );
    Burgers::v.resize( Nx, Ny );
}

double Burgers::x(size_t col) {
    const double dx = m -> getDx();
    const double Lx = m -> getLx();
    return (dx * col - Lx / 2);
}

double Burgers::y(size_t row) {
    const double dy = m -> getDy();
    const double Ly = m -> getLy();
    return (dy * row - Ly / 2);
}

void Burgers::initializeVelocityField() {
    // use a submatrix here
    for( size_t col=0UL; col<u.rows(); ++col ) {
        for( size_t row=0UL; row<u.columns(); ++row) {
            const double r = sqrt(pow(x(col),2)+pow(y(row),2));
            u(row,col) = (r < r_thresh ? 2*pow(1-r,4)*(4*r+1) : 0.0);
        }
    }
    v = u;
}

void Burgers::printVelocityField() {
    const auto dotStep = static_cast<const size_t> (floor(u.rows() / 10.0));
    const char filename[] = "velocity_u.csv";
    cout << "Writing velocity field data to file - " << filename << endl;
    ofstream dataFile (filename, fstream::trunc);
    for( size_t row=0UL; row<u.rows(); ++row ) {
        for( size_t col=0UL; col<u.columns(); ++col) {
            dataFile << (col==0 ? ' ' : ',') << u(row,col);
        }
        dataFile << '\n';
        if (row%dotStep == 0) cout << '.';
    }
    dataFile.close();
    cout << "\nDone :)" << endl;
}

void Burgers::integrateVelocityField() {
    const auto nsx = m->getNx() - 2;
    const auto nsy = m->getNy() - 2;
    const double Nt = m->getNt();
    const double dx = m->getDx();
    const double dy = m->getDy();
    const double dt = m->getDt();
    const double ax  = m->getAx();
    const double ay  = m->getAy();
    const double b  = m->getB();
    const double c  = m->getC();

    auto u_ix_j = blaze::submatrix(u, 1UL, 2UL, nsx, nsy);
    auto u_xi_j = blaze::submatrix(u, 1UL, 0UL, nsx, nsy);
    auto u_i_j  = blaze::submatrix(u, 1UL, 1UL, nsx, nsy);
    auto u_i_xj = blaze::submatrix(u, 0UL, 1UL, nsx, nsy);
    auto u_i_jx = blaze::submatrix(u, 2UL, 1UL, nsx, nsy);
    
    auto v_ix_j = blaze::submatrix(v, 1UL, 2UL, nsx, nsy);
    auto v_xi_j = blaze::submatrix(v, 1UL, 0UL, nsx, nsy);
    auto v_i_j  = blaze::submatrix(v, 1UL, 1UL, nsx, nsy);
    auto v_i_xj = blaze::submatrix(v, 0UL, 1UL, nsx, nsy);
    auto v_i_jx = blaze::submatrix(v, 2UL, 1UL, nsx, nsy);
    
    const double p_xi_j = c / dx / dx + ax / dx;
    const double p_ix_j = c / dx / dx;
    const double p_i_j  = -2.0*c*(1/dx/dx + 1/dy/dy) - ax/dx - ay/dy + 1/dt;
    const double p_i_xj = c / dy / dy + ay / dy;
    const double p_i_jx = c / dy / dy;
    
    blaze::DynamicMatrix<double> U34 (nsx, nsy), U23 (nsx, nsy), V34 (nsx, nsy), V23 (nsx, nsy);
    
    for (int i = 1; i <= Nt; ++i) {
        U34 = p_xi_j * u_xi_j + p_ix_j * u_ix_j + p_i_j * u_i_j + p_i_xj * u_i_xj + p_i_jx * u_i_jx;
        U23 = -1.0 * b / dx * u_i_j % (u_i_j - u_xi_j) - b / dy * v_i_j % (u_i_j - u_i_xj);
        
        V34 = p_xi_j*v_xi_j + p_ix_j*v_ix_j + p_i_j*v_i_j + p_i_xj*v_i_xj + p_i_jx*v_i_jx;
        V23 = ( -1.0 * b / dx ) * u_i_j % (v_i_j - v_xi_j) - ( b / dy * v_i_j ) % (v_i_j - v_i_xj);
    
        u_i_j = dt * (U34 + U23);
        v_i_j = dt * (V34 + V23);
    }
    cout << "\n\nDone\n\n";
}

double Burgers::fieldEnergy() {
    return 0.5*sum(blaze::sqrt(u*trans(u) + v*trans(v)))*m->getDx()*m->getDy();
}