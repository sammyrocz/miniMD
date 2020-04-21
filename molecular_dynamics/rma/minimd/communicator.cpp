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
    MPI_Barrier(gcomm);
    double stime = MPI_Wtime();
    if(rcvrank == rank){
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

    stime = MPI_Wtime() - stime;
    MPI_Barrier(gcomm);
    double time;
    MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, ucomm);
    commtime[aindex] += time;
}

void Communicator::rma(void *temp, long long int &atoms, int dimenstion, int ts, int rank, int aindex)
{

    long long int disp = 0;
    long long int size = 0;
    MPI_Scan(&atoms, &disp, 1, MPI_LONG_LONG_INT, MPI_SUM, gcomm);
    double **array;
    MPI_Win win;
    if (rank == rcvrank)
    {
        double **array = (double **)temp;
        size = disp * dimenstion;
        MPI_Win_create(array[ts], size * sizeof(double), sizeof(double), MPI_INFO_NULL, gcomm, &win);
        atoms = disp; // updating the atom to be reflected in coanalysis
    }
    else
    {
        MPI_Win_create(MPI_BOTTOM, 0, sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }
    
    
    MPI_Win_fence(0, win);
    MPI_Barrier(gcomm);
    double stime = MPI_Wtime();
    if (rank != rcvrank)
    {

        double *array = (double *)temp;
        size = atoms * dimenstion;
        disp = (disp - atoms);
        MPI_Put(array, size, MPI_DOUBLE, rcvrank, disp * dimenstion, size, MPI_DOUBLE, win);
    }
    stime = MPI_Wtime() - stime;
    MPI_Barrier(gcomm);
    MPI_Win_fence(0, win);
    
    double time;
    MPI_Allreduce(&stime, &time, 1, MPI_DOUBLE, MPI_MAX, ucomm);
    commtime[aindex] += time;
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
