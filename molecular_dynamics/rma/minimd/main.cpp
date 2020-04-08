#include <stdio.h>
#include <stdlib.h>
#include "ljs.h"
#include "mpi.h"
#include "util.h"
#include "groupcomm.h"
#include "universecomm.h"
#include "modalysis.h"
#include "communicator.h"

UniverseComm *UniverseComm::instance = 0;
GroupComm *GroupComm::instance = 0;
int Communicator::commtype = 0;
int Communicator::nana = 0;
int Communicator::nsim = 0;
int Communicator::rcvrank = 0;

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    int me;
    int nprocs;
    int proctype;
    char *commtype = 0;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    char *pratio = NULL; // simulation to analysis processes ratio
    
    for (int i = 0; i < argc; i++)
    {

        if ((strcmp(argv[i], "-npr") == 0))
        {
            if (pratio == NULL)
                pratio = new char[256];
            if (pratio != NULL)
                strcpy(pratio, argv[++i]);
            else
                printf("Error:Ratio not specified!\n");
            continue;
        } else if(strcmp(argv[i],"-comm") == 0){
                char *temp = new char[256];
                strcpy(temp,argv[++i]);

                if(strcmp(temp,"SEND") == 0){
                    Communicator::commtype = 0;  // communication type is SEND/RECV
                } else {
                    Communicator::commtype = 1; // communcation type is RMA
                }
        }
    }

    int nsim = -1;
    int nana = -1;

    MPI_Comm intercomm; // intracommunicator
    MPI_Comm intracomm; //

    
    if (pratio != 0)
    {

        // Dividing the communicator

        calstar(pratio, nsim, nana);
        
        Communicator::nana = nana;
        Communicator::nsim = nsim;    
        Communicator::rcvrank = nsim;

        if (nprocs % (nsim + nana) != 0)
        {

            if (me == 0)
                printf("Error: Correctly specify the ratio\n");
            MPI_Finalize();
            exit(0);
        }

        //plain processes till here
        joinuniverse(me, nsim, nana, MPI_COMM_WORLD, intercomm);
        GroupComm::setup(intercomm);
        GroupComm::getinstance()->isAvailable = true;
        GroupComm::getinstance()->init(nsim, nana, (nsim + nana));
        splituniverse(me, nsim, nana, proctype, MPI_COMM_WORLD, intracomm);
        UniverseComm::setup(intracomm);
    }

    if (pratio != 0)
    {

        if (proctype == SIMULATION)
        {

            miniMDinit(argc, argv, commtype);
        }
        else
        {
            #ifdef MODALYSIS

            Modalysis mod;
            mod.init(argc,argv, commtype);
            mod.process();

            #endif
        }
    }
    else
    {
        miniMDinit(argc, argv, commtype);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}