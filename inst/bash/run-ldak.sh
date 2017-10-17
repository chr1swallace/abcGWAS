#!/bin/bash

d=$1
log=${d}/ldak.log

bedfile="${d}/ldak"
if [ ! -e ${bedfile}.bed ]
then
  echo " ${bedfile}.bed does not exist"
fi

echo "running ldak, see $log for logs"

ldak --cut-weights $d \
     --bfile $bedfile > $log

export ncut=`cat $d/section_details.txt|tail -1|cut -f1 -d' '`

for j in $(seq "$ncut");
do
 export iii=$j
 ldak \
 --calc-weights $d \
 --bfile $bedfile \
 --section ${iii} >> $log
done



## join the weight files
ldak \
 --join-weights $d >> $log



