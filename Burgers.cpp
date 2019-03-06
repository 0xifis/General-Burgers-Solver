#include "Burgers.h"
#include "math.h"
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