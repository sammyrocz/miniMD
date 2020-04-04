/*
    author: sanjay kumar
    module: modalysis
    comments: entry point for modalysis
*/
#include "modalysis.h"
#include "groupcomm.h"
#include "universecomm.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#define FILENAMELEN 64
#define ANAMELEN 16

Modalysis::Modalysis()
{

    gcomm = GroupComm::getinstance()->comm;
    ucomm = UniverseComm::getinstance()->comm;
    MPI_Comm_rank(gcomm, &grank);
    MPI_Comm_rank(ucomm, &myrank);

    atoms = 0;
    config = NULL;
}

Modalysis::~Modalysis(){
    // destructing the created objects
}

void Modalysis::allocate()
{
    // allocates required space for the communication

    adim = new int[anum];
    atevery = new int[anum];
    atsteps = new int[anum];
    acurrstep = new int[anum];
    afreq = new int[anum];
    afname = new char[anum * FILENAMELEN];
	aname = new char[anum * ANAMELEN];

    int simpc = GroupComm::getinstance()->simpc;
    array = new double *[anum];

    for (int i = 0; i < anum; i++)
    {
        array[i] = 0;
    }

    xoriginal = 0;
    voriginal = 0;
    msd = new double[4];
    vacf = new double[4];
}

void Modalysis::init(int argc, char **argv, char *commt)
{

    num_steps = 500;
    commtype = commt;
    int i;
    
    for (i = 0; i < argc; i++)
    {

        if(strcmp(argv[i], "-atoms") == 0){
            atoms = atoll(argv[++i]);
        }else if(strcmp(argv[i],"-acfg") == 0){
            config = new char[256];
            strcpy(config, argv[++i]);        
        }
    }

    if (config == NULL)
    {
        printf("Config file NULL error %d\n", ucomm);
        exit(1);
    }

    if (myrank == 0)
    {

        FILE *fp = fopen(config, "r");

        if (fp == NULL)
        {
            printf("Config file open error %d %s\n", errno, strerror(errno));
            exit(1);
        }

        fscanf(fp, "%d", &anum);
        allocate();

        i = 0;
        while (i < anum)
        {
            fscanf(fp, "%d %d %d %d %s %s", &adim[i], &atevery[i], &atsteps[i], &acurrstep[i], afname + i * FILENAMELEN, aname + i * ANAMELEN);
            ++i;
        }

        fclose(fp);
    }

    MPI_Bcast(&anum, 1, MPI_INT, 0, ucomm);

    if (afreq == NULL)
		allocate();

	if (MPI_Bcast(adim, anum, MPI_INT, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis config adim bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atevery, anum, MPI_INT, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis config atevery bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atsteps, anum, MPI_INT, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis config atsteps bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(acurrstep, anum, MPI_INT, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis config acurrstep bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(afname, anum * FILENAMELEN, MPI_CHAR, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis dump file name bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(aname, anum * ANAMELEN, MPI_CHAR, 0, ucomm) != MPI_SUCCESS)
		printf("\nAnalysis name bcast error %d %s\n", errno, strerror(errno));

    for (i = 0; i < anum; i++)
	{
		afreq[i] = num_steps / atsteps[i]; //todo: deal with integrality later
	}

        

}