#ifndef DUMP_H_
#define DUMP_H_

#include <stdint.h>
#include "mpi.h"
#include "atom.h"
#include "comm.h"

#define ANAMELEN 16
#define FILENAMELEN 64

class Dump 
{

	char *dumpfile; //[] = "dump.txt";
	char *posfile, *posfilename; // = "positions.txt";
	char *velfile, *velfilename; // = "velocities.txt";
	//FILE *dumpfp, *configfp;

	int output_frequency;
	int num_steps;

	MMD_int nlocal; //long long int nlocal;
	MMD_int numAtoms;
	MMD_int totalAtoms;
	MMD_int bufsize;
	int count, rcount;

	MPI_Offset mpifo;
	MPI_File posfh, velfh;
	
	MPI_Datatype dtype;

	MPI_File *afh;
	MMD_float *array; 
	MMD_int arraylen;
	char *configFile;

	int *adim;
	int *atevery;
	int *atsteps;
	int *acurrstep;
	char *afname;
	char *aname;
	double *time;
	
	

  public:

  Dump();
  ~Dump();
	double * time_to_write;
	double ttime;
	int getFreq();
	char *getConfigFile();
  void initDump(Comm &, int, int, char *, char *);
	void writeFile(Atom &, int, Comm &);
	void pack(Atom &, int, Comm &);
	void dump(Atom &, int, Comm &);
	void unpack(void);
	void finiDump(Comm &);

	void initAnalysisDump(Comm &, char *);
	void finiAnalysisDump();
	void aalloc(int);
	
	void apack(Atom &, Comm &, int, int);
	void adump(Atom &, Comm &, int, int);
	void aunpack();
	void updateConfig(Comm &);
	void writeAOutput(Atom &, Comm &, int, int);

	MMD_float *pos, *vel, *rtest;
	//char *dumpdir;
	
	int anum;
	int *afreq;

};

#endif

