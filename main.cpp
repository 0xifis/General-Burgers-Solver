#include <chrono>
#include "Model.h"
#include "Burgers.h"
#include "mpi.h"
#include "MyMPI.h"
#include "helper.h"

using namespace std;

int main(int argc, char* argv[]) {
    MyMPI myMPI(argc, argv);
    if(!myMPI.isValid()) return 1;
    MPI_Comm myComm;
    myMPI.createSubComm(&myComm);
    
    Model m(argc, argv);
    if (m.isVerbose() && myMPI.rank()==0) m.printParameters();
    if (m.isValid()) {
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Processor #%i is now running the model with given parameters.\n", myMPI.rank());
        Burgers b(&m);

        // Call code to initialise the problem here
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
        

        // Calculate final energy and write output
        if(myMPI.rank()==0) {
            cout << "Energy: " << b.fieldEnergy() << endl;
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


