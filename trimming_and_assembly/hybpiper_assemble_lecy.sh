#!/usr/bin/bash
conda activate HybPiper
#call like this: nohup hybpiper_assemble [target file] [sample list] &> logs/hybpiper_assemble_[target type]_nohup.log

targetFile=$1
targetType=${targetFile%_targets.fa}
sampleList=$2

while read sampleCode; do
pairedReads="/data_harddisk/diana/a353_new/00_read_process/$sampleCode/${sampleCode}.clean_[12].fq.gz"
singleReads="/data_harddisk/diana/a353_new/00_read_process/$sampleCode/${sampleCode}.clean_S.fq.gz"
echo "Starting assembly on $sampleCode at $(date) with target file $targetFile"
echo
hybpiper assemble --cpu 50 --bwa -t_dna $targetFile --readfiles $pairedReads --unpaired $singleReads --merged --prefix $sampleCode &> logs/${sampleCode}_hybpiper_assemble_${targetType}.log
perl -nle 'print if /^\[.+\]:\s+/' logs/${sampleCode}_hybpiper_assemble_${targetType}.log
echo
echo "HybPiper assembly on $sampleCode finished at $(date)"
echo
echo "--------------------------------------------------------------------"
echo
done < $sampleList
