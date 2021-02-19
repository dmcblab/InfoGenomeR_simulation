#!/bin/bash
declare -a arr;
k=0;
for i in 15 14 13 12 11
do
	arr[$k]="echo $i &&sleep $i";
	k=$(($k+1));
done

for z in `seq 0 22`;do
	echo $z;
done
