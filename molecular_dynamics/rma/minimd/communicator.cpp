#include "communicator.h"

Communicator::Communicator()
{

    acount = new long long int[(nana + nsim)];

    rcvrank = nsim;
}

Communicator::~Communicator()
{
}

void Communicator::sendrecv(void *temp, long long int atoms, int dimenstion, int ts, int rank)
{
    MPI_Gather(&atoms, 1, MPI_LONG_LONG_INT, acount, 1, MPI_LONG_LONG_INT, rcvrank, gcomm);
    
    if (rcvrank == rank)
    {
        double **array = (double **)temp;

        for (int i = nsim; i > 0; i--)
        {
            acount[i] = acount[i-1];
        }

        acount[0] = 0;

        for(int i = 0 ; i < nsim; i++){
            MPI_Recv((array[ts] + (dimenstion * acount[i])),dimenstion*acount[i+1],MPI_DOUBLE,i,ts,gcomm,MPI_STATUS_IGNORE);
        }

        
        
    }
    else
    {   
        double *array = (double *)temp;
        MPI_Ssend(array, dimenstion * atoms, MPI_DOUBLE, rcvrank, ts, gcomm);
    }
}

void Communicator::communicate(void *temp, long long int atoms, int dimension, int ts, int rank)
{


    

    if (commtype == 0)
    {
        // sendrecv function
        sendrecv(temp, atoms, dimension, ts, rank);
    }
}