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
            acount[i] = acount[i - 1];
        }

        acount[0] = 0;

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

void Communicator::rma(void *temp, long long int atoms, int dimenstion, int ts, int rank)
{

    long long int disp = 0;
    long long int size = 0;
    MPI_Scan(&atoms, &disp, 1, MPI_LONG_LONG_INT, MPI_SUM, gcomm);
    
    MPI_Win win;
    if (rank == rcvrank)
    {
        double **array = (double **)temp;
        size = disp * dimenstion;
        MPI_Win_create(array[ts], size * sizeof(double), sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }
    else
    {
        MPI_Win_create(MPI_BOTTOM, 0, sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }

    MPI_Win_fence(0, win);

    if (rank != rcvrank)
    {

        double *array = (double *)temp;
        size = atoms * dimenstion;
        disp = (disp - atoms);
        MPI_Put(array, size, MPI_DOUBLE, rcvrank, disp * dimenstion, size, MPI_DOUBLE, win);
    }

    MPI_Win_fence(0, win);

    MPI_Win_free(&win);
}

void Communicator::communicate(void *temp, long long int atoms, int dimension, int ts, int rank)
{

    if (commtype == 0)
    {
        // sendrecv function
        sendrecv(temp, atoms, dimension, ts, rank);
    }
    else if (commtype == 1)
    {

        rma(temp, atoms, dimension, ts, rank);
    }
}