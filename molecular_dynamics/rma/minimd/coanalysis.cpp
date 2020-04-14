#include "modalysis.h"

void Modalysis::coanalyze(double **array, int aindex, int ts)
{

    switch (aindex)
    {
    case 0:
        // msd processing

        if (ts == 0)
            for (int i = 0; i < nlocal * adim[aindex]; i++)
                xoriginal[i] = array[ts][i];

        compute_msd(array[ts]);
        break;
    case 1:

        if (ts == 0)
            for (int i = 0; i < nlocal * adim[aindex]; i++)
                voriginal[i] = array[ts][i];

        compute_vacf(array[ts]);
        break;

    case 2:
        compute_histo(array[ts]);
        break;

    case 3:
        compute_histo(array[ts]);
        break;

    case 4:
        for (int k = 0; k < nlocal; k++)
            compute_fft_1d(ts - (acurrstep[aindex] - 1), ts, k, array);
        break;

    case 5:
        for (int k = 0; k < nlocal; k++)
            compute_fft_1d(ts - (acurrstep[aindex] - 1), ts, k, array);
        break;
    }
}
// reads the data from the simulation processes
void Modalysis::readdata(int aindex, int ts)
{

    if (acurrstep[aindex] >= atsteps[aindex])
        return;

    if ((acurrstep[aindex] % atevery[aindex]) == 0 && istemporal[aindex] == 0)
    {

        transmitter.communicate(array[aindex], nlocal, adim[aindex], acurrstep[aindex], grank,aindex);
        
        coanalyze(array[aindex], aindex, acurrstep[aindex]);
    }
    else if (istemporal[aindex])
    {
        transmitter.communicate(array[aindex], nlocal, adim[aindex], acurrstep[aindex], grank,aindex);

        if ((ts + 1) % atevery[aindex] == 0)
        {   
            coanalyze(array[aindex], aindex, acurrstep[aindex]);
        }
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
                nlocal = 0; //clearing up the no of atoms for new iteration
                readdata(j, i);
            }
    }
}
