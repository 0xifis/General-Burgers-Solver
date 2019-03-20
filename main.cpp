#include <chrono>
#include "Model.h"
#include "Burgers.h"
#include "mpi.h"
#include "MyMPI.h"

using namespace std;

int main(int argc, char* argv[]) {
    Model m(argc, argv);
    MyMPI myMPI(argc, argv, m.getPx(), m.getPy());
    
    if (m.isVerbose() && myMPI.rank()==0) m.printParameters();
    
    if(!myMPI.isValid()) return 1;
    
    if (m.isValid()) {
        MPI_Barrier(MPI_COMM_WORLD);
        Burgers b(&m, &myMPI);

        b.initializeVelocityField();

        typedef std::chrono::high_resolution_clock hrc;
        typedef std::chrono::milliseconds ms;
        
        hrc::time_point start = hrc::now();
        
        b.integrateVelocityField();

        hrc::time_point end = hrc::now();
    
        MPI_Barrier(MPI_COMM_WORLD);
    
        cout << "Processor #" << myMPI.rank()
             << ": Integration took "
             << chrono::duration_cast<ms>(end - start).count()/1000.
             << "s.\n";

        b.calculateFieldEnergy();
        if(myMPI.rank()==0) {
            cout << "Energy: " << b.getFieldEnergy() << endl;
            b.printVelocityField();
        }

    } else if(myMPI.rank()==0) {
        cerr << "The model was not run as it is not valid with the given parameters.\n"
             << "Please choose a different set of parameters and try again or use '-h' flag for the help message.\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}


