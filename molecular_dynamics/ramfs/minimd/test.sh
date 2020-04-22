#!/bin/bash

cp ../config.sample .
make -f Makefile.default || exit
mpirun -np 4 ./miniMD -acfg config.sample -aoutp /tmp
