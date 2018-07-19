#!/bin/sh

#  simulation_loop.sh
#  
#
#  Created by james roach on 8/25/15.
#


#get current directory


#echo $1
#echo $2
#echo $3
cwd=$(pwd)

#compile code
icpc "${cwd}/network_dyn.cpp" -O2 -o run1.out
icpc "${cwd}/AMD_v4.cpp" -O2 -o run2.out

#C=$(awk 'BEGIN{for(i=0.0;i<=1.5;i+=0.25)print i}')
C=$(awk 'BEGIN{for(i=0.216;i<=0.25;i+=0.24)print i}')
D=$(awk 'BEGIN{for(i=20.0;i<=10020.05;i+=1000)print i}')
for k in $C
do
    for l in $D
    do
        for o in {1..5}
        do
            #sleep `expr $RANDOM % 10`

            #make excecution directory
            newdir=${cwd}"/gks_"${k}"net_"${l}"run_"${o}
            mkdir "$newdir"
            #cd to execution directory
            cd "$newdir"
            cp "${cwd}/curr.txt"    "${newdir}/curr.txt"
            cp "${cwd}/mcond.txt"    "${newdir}/mcond.txt"
            cp "${cwd}/ifg_mat.txt"    "${newdir}/ifg_mat.txt"
            #run code  gks, freq, network, dc spread
            "${cwd}/run1.out" $k 10.0 $l 0.0 0.25 1.0 0.02 0.1
            "${cwd}/run2.out" "raster_dat.txt"
        done
    done
done
