#ifndef MODALYSIS_H
#define MODALYSIS_H

/*
    author: sanjay kumar
    module: modalysis
    comments: Class definition for modalysis
*/

#include "mpi.h"
#include "types.h"
#define PAD 3
#define FILENAMELEN 64

class Modalysis
{

private:
    int myrank;
    long long int nglobal;
    long long int nlocal;
    long long int atoms;
    MPI_Comm gcomm;
    int anum;
    int *adim;
    int *atevery;
    int *atsteps;
    int *acurrstep;

    double **array;
    double ***temporalarr; // temporal data
    int **temporalsize;
    // holds the data and performs computation
    double *xoriginal;
    double *voriginal;

    double *msd;
    double *vacf;
    char *commtype;
    int grank; // group rank
    MPI_Comm ucomm;
    MPI_Offset mpifo;
    // test
    char *afname;
    char **filenames;
    MPI_File *afh;
    int *timer;
    MPI_Win win;
public:
    void init(long long int, char *);
    void coanalyze();
    void configure(); // gets configuration from config file
    void allocate();  // allocates required memory for setup
    void process();
    void checktimestep();
    void readdata();

    void compute_vacf(double *arr);
    void compute_msd(double *);
    void compute_histo(double *);
    void compute_fft_1d(double *, double **, int, long long int);
    void cleanup();
};

#endif