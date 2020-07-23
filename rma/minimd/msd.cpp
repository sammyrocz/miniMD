/*

    author: sanjay kumar
    module : modalysis


*/

#include <mpi.h>
#include "modalysis.h"
#include <stdio.h>
void Modalysis::compute_msd(double *arr)
{

    double dx, dy, dz;

    double vector[4] = {0};

    for (int i = 0; i < nlocal; i++)
    {

        dx = arr[i * PAD + 0] - xoriginal[i * PAD + 0];
        dy = arr[i * PAD + 1] - xoriginal[i * PAD + 1];
        dz = arr[i * PAD + 2] - xoriginal[i * PAD + 2];

        vector[0] += dx * dx;
        vector[1] += dy * dy;
        vector[2] += dz * dz;
        vector[3] += dx * dx + dy * dy + dz * dz;
    }
    MPI_Allreduce(vector, msd, 4, MPI_DOUBLE, MPI_SUM, ucomm);

  

    

    return;
}
