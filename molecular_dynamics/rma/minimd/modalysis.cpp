/*
    author: sanjay kumar
    module: modalysis
    comments: entry point for modalysis
*/
#include "modalysis.h"
#include "groupcomm.h"
#include "universecomm.h"
#include<stdlib.h>
#include<stdio.h>



void Modalysis::init(long long int a,char *commt){
    gcomm = GroupComm::getinstance()->comm;
    ucomm = UniverseComm::getinstance()->comm;
    atoms = a;
    commtype = commt;
    MPI_Bcast(&nglobal,1,MPI_LONG_LONG_INT,0,gcomm); // defined in ljs.cpp [minimd]
    MPI_Comm_rank(gcomm,&grank);
   	MPI_Comm_rank(ucomm,&myrank);
}