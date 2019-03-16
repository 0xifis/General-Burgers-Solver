#ifndef COURSEWORK_HELPER_H
#define COURSEWORK_HELPER_H

#include "mpi.h"
#include "MyMPI.h"
#include <string>
#include <iostream>

using namespace std;

MyMPI::MyMPI(int argc, char* argv[]) {
    for(int i=1; i < argc; ++i ) {
        string arg = string(argv[i]);
        if (arg == "-x") {
            try {
                i++;
                Px = stoi(string(argv[i]));
                i++;
                Py = stoi(string(argv[i]));
                Np = Px * Py;
            }
            catch (const exception &e) {
                cerr << "Error: " << e.what() << endl;
                cerr << "-x option requires two arguments of type int." << endl;
            }
        }
    }

    int mpi_init_err = MPI_Init(&argc, &argv);
    int mpi_rank_err = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int mpi_size_err = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    if (mpi_init_err != MPI_SUCCESS || mpi_rank_err != MPI_SUCCESS  || mpi_size_err != MPI_SUCCESS ) {
        cerr << "An error occurred initialising MPI\n";
        return;
    }
    
    if (world_size < Np) {
        if(world_rank == 0) {
            fprintf(stderr, "Not enough processors are available to split domain as requested. Exiting now.\n");
        }
        return;
    }
    
    if (world_rank==0) {
        printf("Running integration on %i process(es).\n", Np);
    }
    valid = true;
}

void MyMPI::createSubComm(MPI_Comm* myComm) {
    MPI_Group worldGroup, myGroup;
    int* group_ranks = new int[Np];
    for(int i = 0; i < Px*Py; i++) {
        group_ranks[i] = i;
    }
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
    MPI_Group_incl(worldGroup, Px * Py, group_ranks, &myGroup);
    MPI_Comm_create(MPI_COMM_WORLD, myGroup, myComm);
}

#endif