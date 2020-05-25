#include "communicator.h"

Communicator::Communicator()
{

    acount = new long long int[(nana + nsim)];
    pending_request = false;
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

void Communicator::pre_rma(int rank)
{

    if (commtype == 0) // Request Type SEND
        return;

    // Request is of type POST ONLY

    if (pending_request == false) // no pending requests
        return;

    if (rank != rcvrank)
    {
        MPI_Win_wait(win);
        
    } 

    MPI_Win_free(&win);
    MPI_Group_free(&rmagroup);
    MPI_Group_free(&comm_group);

    pending_request = true;
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

    MPI_Igather(&atoms, 1, MPI_LONG_LONG_INT, acount, 1, MPI_LONG_LONG_INT, rcvrank, gcomm,&request);

    if (rank == rcvrank)
    {
        MPI_Win_create(MPI_BOTTOM, 0, sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }
    else
    {
        double *array = (double *)temp;

        MPI_Win_create(array, atoms * dimenstion * sizeof(double), sizeof(double), MPI_INFO_NULL, gcomm, &win);
    }

    MPI_Comm_group(gcomm, &comm_group);

    if (rank == rcvrank)
    {

        int ranks[nsim];

        for (int i = 0; i < nsim; i++)
            ranks[i] = i;

        MPI_Group_incl(comm_group, nsim, ranks, &rmagroup);
        /* Begin the access epoch */
        MPI_Win_start(rmagroup, 0, win);

	MPI_Wait(&request,MPI_STATUS_IGNORE);

     
        for (int i = nsim; i > 0; i--)
        {
            acount[i] = acount[i - 1];
            atoms += acount[i];
        }
        acount[0] = 0;

        /* Put into rank==0 according to my rank */
        double **array = (double **)temp;

        for (int i = 0; i < nsim; i++)
        {

            MPI_Get((array[ts] + (dimenstion * acount[i])), dimenstion * acount[i + 1], MPI_DOUBLE, i, 0, dimenstion * acount[i + 1], MPI_DOUBLE, win);
        }

        /* Terminate the access epoch */
         MPI_Win_complete(win);
       
    }
    else
    {

        int ranks[] = {rcvrank};

        MPI_Group_incl(comm_group, 1, ranks, &rmagroup);
        MPI_Win_post(rmagroup, 0, win);
       
    }

   
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
