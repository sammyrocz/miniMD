#ifndef UTIL_H_
#define UTIL_H_

/* 

author : Sanjay Kumar
module : util
comments : Utility functions for communicators and analysis to processes raio

*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "comm.h"
#include "mpi.h"
#define SIMULATION 1
#define ANALYSIS 0 

// computes the gcd of the two number
int gcd(int a, int b) 
{ 
    // Everything divides 0  
    if (a == 0) 
       return b; 
    if (b == 0) 
       return a; 
   
    // base case 
    if (a == b) 
        return a; 
   
    // a is greater 
    if (a > b) 
        return gcd(a-b, b); 
    return gcd(a, b-a); 
} 

// computes simulation to analysis ratio for each node
void calstar(char *inp,int& nsim,int& nanal){
    
    if(inp == 0)
    return;

    char temp[strlen(inp)];
    strcpy(temp,inp);
    char* token = strtok(temp, ":"); 
    
  
    while (token != NULL) { 
        
        if(nsim == -1){
            nsim = atoi(token);
        } else {
            nanal = atoi(token); 
        }
        token = strtok(NULL, "-");
    }
    
    int _gcd = gcd(nsim,nanal);
    nsim = nsim / _gcd;
    nanal = nanal/_gcd;
}

// creating two universe for miniMD and modalysis
void splituniverse(int me,int simr,int analyr,int &proctype,MPI_Comm oldcomm,MPI_Comm &newcomm){

    int total = simr + analyr;
    int type;
    if((me % total) >= simr){
        type = ANALYSIS; 
    } else {
        type = SIMULATION;
    }
    proctype = type;
    MPI_Comm_split(oldcomm, type, me, &newcomm);
}


// joins node local processes
// si => simulation
// an => analysis

// creating groups containing process for miniMD and modalysis 
void joinuniverse(int me, int si,int an,MPI_Comm oldcomm,MPI_Comm &comm){
    
    // every group has it's own communicator    
    int color = me/(si+an);
    MPI_Comm_split(oldcomm,color,me,&comm);

}



#endif