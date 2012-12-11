#!/bin/bash
msg=('B&B' 'C&B Greedy' 'C&B Clique' 'C&B Greedy y Clique' )
arr=('' '-G -X' '-K -X' '-X -K -G')

len=${#msg[@]}

for file in `ls inst*.lp`
do
    for((i=0;i<$len;i++))
    do
        echo "" >> log
        echo "Corriendo $file con ${msg[i]}"
        echo "Corriendo $file con ${msg[i]}" >> log
        echo "./Tp -f $file -t ${arr[i]} >/dev/null" >> log
        echo "***********">>log
        ./Tp -f $file ${arr[i]} -t >/dev/null 2>> log
        echo "***********">>log
    done
done
