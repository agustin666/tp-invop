#!/bin/bash

msg=('B&B' 'C&B Greedy' 'C&B Clique' 'C&B Greedy y Clique' 'Greedy' 'Clique' 'Greedy y Clique' )
arr=('' '-G -X' '-K -X' '-X -K -G' '-G' '-K' '-G -K')

len=${#msg[@]}

for file in `ls inst*.lp`
do
    for((i=0;i<$len;i++))
    do
        echo "Corriendo $file con ${msg[i]}" >> log
        echo "./Tp -f $file -t ${arr[i]} >/dev/null" >> log
        ./Tp -f $file -t ${arr[i]} >/dev/null 2>> log
    done
done
