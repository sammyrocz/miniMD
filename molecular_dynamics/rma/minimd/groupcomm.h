#ifndef GROUPCOMM_H_
#define GROUPCOMM_H_

#include <stdio.h>
#include "mpi.h"

class GroupComm{

    private:

        static GroupComm * instance;
        
        GroupComm(){
        
            isAvailable = false; 
            simpc = 0;
            analpc = 0;
            gsize = 0;
        }


    public:
    MPI_Win win;
    MPI_Comm comm; // communicator for the group
    bool isAvailable = false;
    int simpc; // simulation process count 
    int analpc; // analysis processes count
    int gsize; // no of processes in the group
    static GroupComm* setup(MPI_Comm comm){
        if(!instance){
            instance = new GroupComm;
            instance->comm = comm;
           
        }
    }

    void init(int csim,int canal,int cgsize){
        simpc = csim;
        analpc = canal;
        gsize = cgsize;
    }

    static GroupComm* getinstance(){
        if(!instance){
            instance = new GroupComm;
            instance->comm = MPI_COMM_WORLD;
        }

        return instance;
    }



};



#endif