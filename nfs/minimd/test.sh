#!/bin/bash

function strict_execute() {

    arrs=$1  #placement for simulation
    arra=$2  #placement for analysis
    atoms=$3 # -nx -ny -nz

    atoms=(${atoms[0]})
    # generating comma seprated string for simulation array
    # temp=${arrs[*]}
    # places=${temp//${IFS:0:1}/,}

    # generating comma seprated string for analysis array
    # temp=${arra[*]}
    # placea=${temp//${IFS:0:1}/,}

    # getting the atom count from -nx -ny -nz
    atomcount=1
    for i in ${atoms[@]}; do
        atomcount=$(expr $atomcount \* $i)
    done
    atomcount=$(expr $atomcount \* 4)

    cp ../config.sample .

    ofprefix=$(echo $arrs + $arra | bc)
    otfname="${ofprefix}.txt"
    atfname="${ofprefix}t.txt"

    cd ../modalysis
    mpirun -np $arra -f aipaddress ./modalysis -atoms $atomcount -acfg config.sample &
    cd ../minimd
    mpirun -np $arrs -f sipaddress ./miniMD -acfg config.sample -nx ${atoms[0]} -ny ${atoms[1]} -nz ${atoms[2]} | tee $otfname

    python3 readtime.py $otfname >>$atfname

}

option=$1

if [ $option -eq 3]; then
    
    rm -rf *t.txt
    rm -rf 3

    mkdir -p 3

    arrs=9
    arra=3
    atoms=(32 64 64)
    head -n 3 hostss3 > sipaddress
    head -n 3 hostsa3 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"

    arrs=18
    arra=6
    atoms=(32 64 64)
    head -n 6 hostss3 > sipaddress
    head -n 6 hostsa3 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"


    arrs=36
    arra=6
    atoms=(32 64 64)
    head -n 12 hostss3 > sipaddress
    head -n 12 hostsa3 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"

    mv *t.txt 3
    

fi


if [ $option -eq 5]; then
    
    rm -rf *t.txt
    rm -rf 5
    mkdir -p 5

    arrs=10
    arra=2
    atoms=(32 64 64)
    head -n 2 hostss5 > sipaddress
    head -n 2 hostsa5 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"

    arrs=20
    arra=4
    atoms=(32 64 64)
    head -n 4 hostss3 > sipaddress
    head -n 4 hostsa3 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"


    arrs=40
    arra=8
    atoms=(32 64 64)
    head -n 8 hostss3 > sipaddress
    head -n 8 hostsa3 > aipaddress
    strict_execute $arrs $arrs "$(echo ${atoms[@]})"

    mv *t.txt 5
    

fi

