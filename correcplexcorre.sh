#!/bin/bash

#msg=('B&B' 'C&B Greedy' 'C&B Clique' 'C&B Greedy y Clique' 'Greedy' 'Clique' 'Greedy y Clique' )
#arr=('' '-G -X' '-K -X' '-X -K -G' '-G' '-K' '-G -K')
#
#len=${#msg[@]}
#
#for file in `ls inst*.lp`
#do
#    for((i=0;i<$len;i++))
#    do
#        echo "" >> log
#        echo "Corriendo $file con ${msg[i]}" >> log
#        echo "./Tp -f $file -t ${arr[i]} >/dev/null" >> log
#        echo "***********">>log
#        ./Tp -f $file -t ${arr[i]} >/dev/null 2>> log
#        echo "***********">>log
#    done
#done

file="./miplib2010-benchmark/reblock67.mps.gz"
echo "" >> log
echo "Corriendo $file con CPLEX " >> log
echo "./Tp -t -f $file -c >/dev/null" >> log
echo "***********">>log
./Tp -t -f $file -c >/dev/null 2>> log
echo "***********">>log

echo "" >> log
echo "Corriendo $file con Greedy y Cover " >> log
echo "./Tp -t -f $file -G -K >/dev/null" >> log
echo "***********">>log
./Tp -f $file -G -K >/dev/null 2>> log
echo "***********">>log

echo "" >> log
echo "Corriendo $file con B&B " >> log
echo "./Tp -t -f $file  >/dev/null" >> log
echo "***********">>log
./Tp -t -f $file >/dev/null 2>> log
echo "***********">>log
