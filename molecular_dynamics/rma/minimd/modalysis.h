#ifndef MODALYSIS_H
#define MODALYSIS_H

/*
    author: sanjay kumar
    module: modalysis
    comments: Class definition for modalysis
*/

#include "mpi.h"
#include "types.h"
#include "communicator.h"
#define PAD 3
#define FILENAMELEN 64

class Modalysis
{

private:


    // *** ----- configure setup ------ *** // 
    char *config;
    int anum; // analysis count
    int num_steps;
    int *adim; // dimesnion of the analysis
    int *atevery; 
    int *atsteps;
    int *acurrstep;
    int *afreq;
    int *istemporal;
    char *afname;
    char *aname;
    int myrank;

    // *** ---- communicators setup ---- *** //
    Communicator transmitter;
    MPI_Comm gcomm;
    MPI_Comm ucomm;

    // *** ----- data for communication ----- *** //

    
    long long int tatoms; // total atoms
    long long int nlocal;
    long long int atoms;

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
    
    MPI_Offset mpifo;
    // test
    char **filenames;
    MPI_File *afh;
    int *timer;
    MPI_Win win;
public:

    Modalysis();
    ~Modalysis();

    void init(int argc, char **argv,char *);
    void coanalyze();
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