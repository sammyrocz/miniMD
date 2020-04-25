#!/bin/bash
make || exit
sleep 0.5
mpirun -np 2 ./modalysis -areadp /tmp -a 131072 -acfg ../minimd/config.sample 
