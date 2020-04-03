#ifndef UNIVERSECOMM_H_
#define UNIVERSECOMM_H_

#include <stdio.h>
#include "mpi.h"

class UniverseComm{

    private:

        static UniverseComm * instance;
       
        UniverseComm(){
         
        }


    public:
    MPI_Comm comm;
    
    static UniverseComm* setup(MPI_Comm comm){
        if(!instance){
            instance = new UniverseComm;
            instance->comm = comm;
        }
    }

    static UniverseComm * getinstance(){
        if(!instance){
            instance = new UniverseComm;
            instance->comm = MPI_COMM_WORLD;
        }

        return instance;
    }



};



#endif