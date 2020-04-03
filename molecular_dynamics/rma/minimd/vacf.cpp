/*

    author: sanjay kumar
    module : modalysis


*/

#include <mpi.h>
#include <stdio.h>
#include "modalysis.h"

void Modalysis::compute_vacf(double *arr) {

	double vxsq, vysq, vzsq;
	double vector[4] = {0};

	for (int i = 0; i < nlocal; i++) {

		vxsq = arr[i*PAD+0] * voriginal[i*PAD+0];
		vysq = arr[i*PAD+1] * voriginal[i*PAD+1];
		vzsq = arr[i*PAD+2] * voriginal[i*PAD+2];
		vector[0] += vxsq;
		vector[1] += vysq;
		vector[2] += vzsq;
		vector[3] += vxsq + vysq + vzsq;

	}

	MPI_Allreduce(vector, vacf, 4, MPI_DOUBLE, MPI_SUM, ucomm);

	
	return;

}

