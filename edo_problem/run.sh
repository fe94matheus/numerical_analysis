#!/bin/bash

make

make clean
mkdir results
for i in 15 31 63 127 255
do

./exe $i 0 1 > resultados_$i.log

done
mv resultados_* results/

