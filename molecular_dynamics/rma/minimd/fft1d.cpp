#include <fftw3.h>
#include <fftw3-mpi.h>
#include <stdlib.h>
#include "modalysis.h"

void Modalysis::compute_fft_1d(double *arr, double **temporalarr, int aindex, long long int natoms)
{

    int index = acurrstep[aindex] % atevery[aindex];
    if (temporalarr[index] != 0)
    {
        delete temporalarr[index];
        temporalarr[index] = 0;
    }
    temporalsize[aindex][index] = natoms;
    temporalarr[index] = new double[nlocal];
    for (int i = 0; i < nlocal; i++)
    {
        temporalarr[index][i] = arr[i];
    }

    if (index < (atevery[aindex] - 1))
    {
        return;
    }

    for (int p = 0; p < nlocal; p++)
    {
        fftw_plan plan;
        double *data;
        fftw_complex *out;
        int alloc_size, i, j, k;

        int dim = atevery[aindex];

        /* get local data size and allocate */
        data = (double *)malloc(sizeof(double) * dim);
        alloc_size = 2 * (dim / 2 + 1);
        out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * alloc_size);
        if (out == NULL || data == NULL)
        {
            printf("Memory allocation failed for out or data: %d\n", alloc_size);
            exit(1);
        }

        /* create plan for forward DFT */
        plan = fftw_plan_dft_r2c_1d(dim, data, out, FFTW_MEASURE);

#ifdef DEBUG
        double t = MPI_Wtime();
#endif

        for (int ts = 0, j = -1; ts < atevery[aindex]; ts++)
        {
            if (temporalsize[aindex][ts] < p)
            {   // only if atom is present 
                data[++j] = temporalarr[ts][p];
            }
        }

#ifdef DEBUG
        t = MPI_Wtime() - t;
        if (myrank < 3)
            printf("%d: FFT: array copy time %lf\n", myrank, t);
#endif

        /* compute transforms */
        fftw_execute(plan);

        fftw_destroy_plan(plan);
        fftw_free(data);
    }
}
