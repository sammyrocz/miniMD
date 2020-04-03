#!/bin/bash
rm config.sample
cp ../config.sample config.sample
make -f Makefile.default || exit 1
mpirun -np 6 ./miniMD -acfg config.sample -npr 2:1 -atom 131072 -comm RMA | tee output.txt
#du -h *.txt