#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include "mpi.h"
class Communicator
{

    /*--------------------------
        Communicator setups the communicaton in the group
    ------------------------------*/
private:
    int *acount; // atoms count
    int *arr;
    MPI_Comm ucomm;
    MPI_Comm gcomm;

    

public:
    static int commtype; // specfies communicaton type
    static int nsim;     // number of simulation process
    static int nana;     // number of analysis process
    static int rcvrank;  // reciver rank
    
    Communicator();
    ~Communicator();
    void communicate(int*, int); 

};

#endif
