/*
 * Developed at Argonne National Laboratory
 *
 * Contact: pmalakar@anl.gov, malakar.preeti@gmail.com
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <mpi.h>

#include "comm.h"
#include "atom.h"
#include "dump.h"

void Dump::initAnalysisDump(Comm &comm, char *analysiscfg)
{

	gcomm = GroupComm::getinstance()->comm;

	transmitter.gcomm = GroupComm::getinstance()->comm;
	transmitter.ucomm = UniverseComm::getinstance()->comm;
	MPI_Comm_rank(gcomm, &grank);

	int i;

	if (analysiscfg == NULL)
	{
		printf("Config file NULL error %d\n", comm.me);
		exit(1);
	}

	configFile = new char[strlen(analysiscfg)];
	strcpy(configFile, analysiscfg);
	//read from config file
	if (comm.me == 0)
	{	

		FILE *fp = fopen(analysiscfg, "r");

		if (fp == NULL)
		{
			printf("Config file open error %d %s\n", errno, strerror(errno));
			exit(1);
		}

		fscanf(fp, "%d", &anum);
		aalloc(anum);

		i = 0;
		while (i < anum)
		{
			fscanf(fp, "%d %d %d %d %d %s %s", &adim[i], &atevery[i], &atsteps[i], &acurrstep[i], &istemporal[i], afname + i * FILENAMELEN, aname + i * ANAMELEN);
			++i;
		}

		fclose(fp);
	}

	MPI_Bcast(&anum, 1, MPI_INT, 0, comm.subcomm);

	transmitter.init(anum);
	//allocate anum elements in all processes but 0 (reader)
	if (afreq == NULL)
		aalloc(anum);

	if (MPI_Bcast(adim, anum, MPI_INT, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis config adim bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atevery, anum, MPI_INT, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis config atevery bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(atsteps, anum, MPI_INT, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis config atsteps bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(acurrstep, anum, MPI_INT, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis config acurrstep bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(istemporal, anum, MPI_INT, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis config temporal bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(afname, anum * FILENAMELEN, MPI_CHAR, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis dump file name bcast error %d %s\n", errno, strerror(errno));
	if (MPI_Bcast(aname, anum * ANAMELEN, MPI_CHAR, 0, comm.subcomm) != MPI_SUCCESS)
		printf("\nAnalysis name bcast error %d %s\n", errno, strerror(errno));

	for (i = 0; i < anum; i++)
	{
		afreq[i] = num_steps / atsteps[i]; //todo: deal with integrality later
	}

	for(int i = 0 ; i < anum ; i++){
		dumper[i] = new double[adim[i]*tatoms];
	}

}

void Dump::apack(Atom &atom, Comm &comm, int n, int aindex)
{

	int ret;
	nlocal = atom.nlocal;

	arraylen = adim[aindex] * nlocal;
	array = (MMD_float *)malloc(arraylen * sizeof(MMD_float));
	if (array == NULL)
		printf("Error: Memory allocation failed");

	MPI_Scan(&nlocal, &numAtoms, 1, MPI_LONG_LONG_INT, MPI_SUM, comm.subcomm);
	MPI_Allreduce(&numAtoms, &totalAtoms, 1, MPI_LONG_LONG_INT, MPI_MAX, comm.subcomm);

	//fix
	for (long long int i = 0; i < nlocal; i++)
	{
		if (aindex == 0)
			array[i * PAD + 0] = atom.x[i * PAD + 0], array[i * PAD + 1] = atom.x[i * PAD + 1], array[i * PAD + 2] = atom.x[i * PAD + 2];
		else if (aindex == 1)
			array[i * PAD + 0] = atom.v[i * PAD + 0], array[i * PAD + 1] = atom.v[i * PAD + 1], array[i * PAD + 2] = atom.v[i * PAD + 2];
		else if (aindex == 2)
			array[i] = atom.x[i * PAD + 0];
		else if (aindex == 3)
			array[i] = atom.x[i * PAD + 1];
		else if (aindex == 4)
			array[i] = atom.v[i * PAD + 0];
		else if (aindex == 5)
			array[i] = atom.v[i * PAD + 1];
	}


	for (long long int i = 0; i < nlocal; i++)
	{
		if (aindex == 0)
			dumper[aindex][i * PAD + 0] = atom.x[i * PAD + 0], dumper[aindex][i * PAD + 1] = atom.x[i * PAD + 1], dumper[aindex][i * PAD + 2] = atom.x[i * PAD + 2];
		else if (aindex == 1)
			dumper[aindex][i * PAD + 0] = atom.v[i * PAD + 0], dumper[aindex][i * PAD + 1] = atom.v[i * PAD + 1], dumper[aindex][i * PAD + 2] = atom.v[i * PAD + 2];
		else if (aindex == 2)
			dumper[aindex][i] = atom.x[i * PAD + 0];
		else if (aindex == 3)
			dumper[aindex][i] = atom.x[i * PAD + 1];
		else if (aindex == 4)
			dumper[aindex][i] = atom.v[i * PAD + 0];
		else if (aindex == 5)
			dumper[aindex][i] = atom.v[i * PAD + 1];
	}



	transmitter.pre_rma(grank); // returns if not RMA
	transmitter.communicate(array, nlocal, adim[aindex], acurrstep[aindex], grank, aindex);
}

void Dump::aunpack()
{
	if (array == NULL)
		delete array;
}

void Dump::writeAOutput(Atom &atom, Comm &comm, int n, int aindex)
{

	if (acurrstep[aindex] >= atsteps[aindex])
	{
		return;
	}

	apack(atom, comm, n, aindex);
	aunpack();

	acurrstep[aindex]++;
}

void Dump::aalloc(int anum)
{
	adim = new int[anum];
	atevery = new int[anum];
	atsteps = new int[anum];
	acurrstep = new int[anum];
	afreq = new int[anum];
	afname = new char[anum * FILENAMELEN];
	aname = new char[anum * ANAMELEN];
	istemporal = new int[anum];

	afh = (MPI_File *)malloc(anum * sizeof(MPI_File));

	dumper = new double*[anum];
	
}

void Dump::finiAnalysisDump()
{

	delete adim;
	delete atevery;
	delete atsteps;
	delete acurrstep;
	delete afreq;
	delete afname;
	delete aname;
	free(afh);
}

char *Dump::getConfigFile()
{
	return configFile;
}
