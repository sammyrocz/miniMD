#include "communicator.h"

Communicator::Communicator()
{

    acount = new long long int[(nana + nsim)];
    rcvrank = nsim;
}

void Communicator::init(int a)
{
    anum = a;
    commtime = new double[anum];

    for (int i = 0; i < anum; i++)
        commtime[i] = 0.0;
}

Communicator::~Communicator()
{
}

void Communicator::sendrecv(void *temp, long long int &atoms, int dimenstion, int ts, int rank, int aindex)
{
    MPI_Gather(&atoms, 1, MPI_LONG_LONG_INT, acount, 1, MPI_LONG_LONG_INT, rcvrank, gcomm);

    if (rcvrank == rank)
    {

        for (int i = nsim; i > 0; i--)
        {
            acount[i] = acount[i - 1];
            atoms += acount[i];
        }
        acount[0] = 0;
    }

    if (rcvrank == rank)
    {
        double **array = (double **)temp;
        for (int i = 0; i < nsim; i++)
        {
            MPI_Recv((array[ts] + (dimenstion * acount[i])), dimenstion * acount[i + 1], MPI_DOUBLE, i, ts, gcomm, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        double *array = (double *)temp;
        MPI_Ssend(array, dimenstion * atoms, MPI_DOUBLE, rcvrank, ts, gcomm);
    }
}

void Communicator::rma(void *temp, long long int &atoms, int dimenstion, int ts, int rank, int aindex)
{

    MPI_Win win;

    MPI_Gather(&atoms, 1, MPI_LONG_LONG_INT, acount, 1, MPI_LONG_LONG_INT, rcvrank, gcomm);

    if (rank == rcvrank)
    {
        MPI_Win_create(MPI_BOTTOM, 0, sizeof(double), MPI_INFO_NULL, gcomm, &win);

        for (int i = nsim; i > 0; i--)
        {
            acount[i] = acount[i - 1];
            atoms += acount[i];
        }
        acount[0] = 0;
    }
    else
    {
        double *array = (double *)temp;

        MPI_Win_create(array, atoms * dimenstion * sizeof(double), sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }


   
        
         MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rcvrank, 0, win);

    if (rank == rcvrank)
    {
        double **array = (double **)temp;

        for (int i = 0; i < nsim; i++)
        {

            MPI_Get((array[ts] + (dimenstion * acount[i])), dimenstion * acount[i + 1], MPI_DOUBLE, i, 0, dimenstion * acount[i + 1], MPI_DOUBLE, win);
        }

    }

        MPI_Win_unlock(rcvrank, win);
      
   
    
   
    MPI_Win_free(&win);


}

void Communicator::communicate(void *data, long long int &atoms, int dimension, int ts, int rank, int aindex)
{

    // data holds the buffer for both sender and receiver, atoms is the total no of atoms ( reciver will have 0 atoms at this point)
    if (commtype == 0)
    {
        // sendrecv function
        sendrecv(data, atoms, dimension, ts, rank, aindex);
    }
    else if (commtype == 1)
    {

        rma(data, atoms, dimension, ts, rank, aindex);
    }
}
