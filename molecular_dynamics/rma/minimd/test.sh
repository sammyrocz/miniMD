#!/bin/bash
rm config.sample
cp ../config.sample config.sample
make -f Makefile.default || exit 1
nx=32
ny=32
nz=32
np=6
mpirun -np $np ./miniMD -acfg config.sample -nx $nx -ny $ny -nz $nz -npr 2:1 -comm RMA | tee output.txt
#du -h *.txt