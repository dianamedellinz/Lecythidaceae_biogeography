#!/bin/bash

gene=$1

t=$2

probes=$3

# transforming n for -
perl -pe 'tr/n/-/ if /^(?!>)/' ../01_alignment/$t/aligned_${probes}_${gene}_${t}.fa > $t/gaps_${probes}_${gene}_${t}.fa

#remove columns that are gaps

trimal -in $t/gaps_${probes}_${gene}_${t}.fa -out $t/noallgaps_${probes}_${gene}_${t}.fa -noallgaps

#trimming
trimal -in $t/noallgaps_${probes}_${gene}_${t}.fa -out $t/trimmed_${probes}_${gene}_${t}.fa -gappyout -sgc -sgt -htmlout $t/trimal.${probes}_${gene}_${t}.html > $t/trimal.${probes}_${gene}_${t}.log

