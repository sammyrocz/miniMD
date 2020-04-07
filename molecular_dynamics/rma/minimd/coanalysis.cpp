#include "modalysis.h"

// void Modalysis::readdata()
// {

//     int aindex;
//     nlocal = 0;
//     long long int arraylen;

//     MPI_Recv(&aindex, 1, MPI_INT, 0, 0, gcomm, MPI_STATUS_IGNORE);
//     MPI_Allreduce(&nlocal, &nglobal, 1, MPI_LONG_LONG_INT, MPI_SUM, gcomm);

//     nlocal = nglobal;

//     arraylen = adim[aindex] * nlocal;
//     array[aindex] = (double *)malloc(arraylen * sizeof(double));
//     if (strcmp(commtype, "SEND") == 0)
//     {
//         int simpc = GroupComm::getinstance()->simpc;
//         long long int natom;
//         long long int tatom = 0;

//         for (int i = 0; i < simpc; i++)
//         {
//             MPI_Recv(&natom, 1, MPI_LONG_LONG_INT, i, 0, gcomm, MPI_STATUS_IGNORE);
//             long long int count = natom * adim[aindex];
//             MPI_Recv(array[aindex] + tatom, count, MPI_DOUBLE, i, 0, gcomm, MPI_STATUS_IGNORE);
//             tatom = tatom + count;
//         }
//     }
//     else if (strcmp(commtype, "RMA") == 0)
//     {
//         long long int disp, temp;
//         disp = temp = 0;
//         MPI_Scan(&temp, &disp, 1, MPI_LONG_LONG_INT, MPI_SUM, gcomm);
//         MPI_Win_create(array[aindex], arraylen * sizeof(double), sizeof(double), MPI_INFO_NULL, gcomm, &win);

//         MPI_Win_fence(0,win);

//         MPI_Win_fence(0,win);

//         MPI_Win_free(&win);
//     }

//     if (acurrstep[aindex] == 0)
//     {

//         if (aindex == 0)
//         {
//             long long int size = adim[aindex] * atoms;
//             xoriginal = new double[size];
//             for (long long int i = 0; i < size; i++)
//             {
//                 xoriginal[i] = 0;
//             }

//             for (int i = 0; i < adim[aindex] * nlocal; i++)
//             {
//                 xoriginal[i] = array[aindex][i];
//             }
//         }
//         else if (aindex == 1)
//         {
//             long long int size = adim[aindex] * atoms;
//             voriginal = new double[size];
//             for (long long int i = 0; i < size; i++)
//             {
//                 voriginal[i] = 0;
//             }
//             for (int i = 0; i < adim[aindex] * nlocal; i++)
//             {
//                 voriginal[i] = array[aindex][i];
//             }
//         }
//     }

//     switch (aindex)
//     {
//     case 0:

//         compute_msd(array[aindex]);
//         acurrstep[aindex] = acurrstep[aindex] + atevery[aindex];
//         break;

//     case 1:
//         compute_vacf(array[aindex]);
//         acurrstep[aindex] = acurrstep[aindex] + atevery[aindex];
//         break;

//     case 2:
//     case 3:
//         compute_histo(array[aindex]);
//         acurrstep[aindex] = acurrstep[aindex] + atevery[aindex];
//         break;
//     case 4:
//     case 5:
//         compute_fft_1d(array[aindex], temporalarr[aindex], aindex, arraylen);
//         acurrstep[aindex]++;
//         break;
//     }

//     delete array[aindex];

//     MPI_Barrier(ucomm);
// }

void Modalysis::readdata(int aindex, int n)
{   
    
    
    if (acurrstep[aindex] >= atsteps[aindex]) return;
    
    if((acurrstep[aindex] % atevery[aindex]) == 0 && istemporal[aindex] == 0){
        
        transmitter.communicate(array[aindex],nlocal,adim[aindex],acurrstep[aindex],grank);
    } else if(istemporal[aindex]){
        transmitter.communicate(array[aindex],nlocal,adim[aindex],acurrstep[aindex],grank);
    }

    acurrstep[aindex]++;
}

void Modalysis::process()
{

    for (int i = 0; i <= num_steps; i++)
    {
        for (int j = 0; j < anum; j++)
            if (i % afreq[j] == 0)
            {
                readdata(j, i);
                nlocal = 0;
            }
    }
}

void Modalysis::coanalyze()
{

    process();
}