conda activate HybPiper

# change to the hybpiper assembly directory
cd ~/diana/a353_new/01_hybpiper/

# Get a353_new target file
# this is the name of the file a353_new/01_hybpiper/lecy_a353_new_targets.fa

# Get a list of samples
ls -d ../00_read_process/L* | cut -d'/' -f3 > sample_list.txt

# shift+G to see the end of the sample_list.txt file
# check the number of samples
# less sample_list.txt -N intro
# 96 lecy samples

# Make a logs directory
mkdir logs

# Run hybpiper assemble for a353_new targets
nohup bash scripts/hybpiper_assemble_lecy.sh lecy_a353_targets.fa sample_list.txt &> logs/hybpiper_assemble_a353_nohup.log &

# Try Hybpiper again for samples with less than 800 genes with recovered sequences
# with no read merging and --cov_cutoff 2
# make list of samples
# sed -n '3p; 6p; 13p; 19,20p; 23p; 27,28p; 35p; 39p; 81,82p' sample_list.txt > no_merge_sample_list.txt

# nohup bash ../scripts/hybpiper_assemble_no_merge.sh mb2_targets.fa &> ../logs/hybpiper_assemble_mb2_no_merge_nohup.log &
nohup bash scripts/hybpiper_assemble_no_merge.sh lecy_a353_targets.fa &> logs/hybpiper_assemble_a353_no_merge_nohup.log &

# Collect stats from mb2 target assemblies
targetFile=lecy_a353_targets.fa
hybpiper stats -t_dna $targetFile gene --seq_lengths_filename a353_exon_lengths --stats_filename a353_exons_stats sample_list.txt
#hybpiper stats -t_dna $targetFile supercontig --seq_lengths_filename a353_supercontig_lengths --stats_filename a353_supercontigs_stats sample_list.txt

hybpiper recovery_heatmap --figure_length 60 --gene_text_size 2 --heatmap_filename a353_exon_recovery_heatmap --heatmap_filetype pdf a353_exon_lengths.tsv

# Collect stats on paralogs
hybpiper paralog_retriever -t_dna $targetFile --figure_length 60 --gene_text_size 4 \
--paralog_report_filename a353_paralog_report --paralogs_above_threshold_report_filename a353_paralogs_above_threshold_report.txt \
--heatmap_filename a353_paralog_heatmap --heatmap_filetype pdf sample_list.txt

#Get a list of genes with no paralogs detected and exclude an leaf or nod gene as well as mb2_gene293, which only had sequence from 3 samples
grep "^>" mb2_targets.fa | sed -e 's/^>Albizia-//' > mb2_gene_list.txt
tail -n+139 mb2_paralogs_above_threshold_report.txt > mb2_genes_with_paralogs.txt
comm -23 <(sort mb2_gene_list.txt) <(sort mb2_genes_with_paralogs.txt) | head -576 | grep -v "mb2_gene293" > mb2_genes_without_paralogs.txt

#get targets for just genes without paralogs
grep -A1 -f mb2_genes_without_paralogs.txt  mb2_targets.fa | sed '/^--$/d' > mb2_single_copy_targets.fa

#retrieve sequences for genes without paralogs
hybpiper retrieve_sequences -t_dna lecy_a353_targets.fa --sample_names sample_list.txt --fasta_dir retrieved_exons_a353 dna


# Alignment
#make parallel lists of sample codes and tip labels to change the sequence ids for alignments
cd ~/acacia/hybpiper_assembly
cut -f2 ../read_processing/ssp01_acacia_lib_info.tsv | tail -n+2 > ../alignments/sample_code_list.txt
cut -f2 ../read_processing_outgroups/mimosoid_outgroups_lib_info.tsv | tail -n+2 >> ../alignments/sample_code_list.txt
cut -f7 ../read_processing/ssp01_acacia_lib_info.tsv | tail -n+2 > ../alignments/tip_label_list.txt
cut -f7 ../read_processing_outgroups/mimosoid_outgroups_lib_info.tsv | tail -n+2 >> ../alignments/tip_label_list.txt

cd ~/acacia/alignments
test -d mb2 || mkdir mb2
cd mb2
test -d unaligned || mkdir unaligned

#copy unaligned fastas and change seq ids from sample codes to taxon names
cp ../../hybpiper_assembly/mb2/mb2_genes_without_paralogs.txt mb2_genes_without_paralogs.txt
while read g; do
perl -pe 's/^(>(A|OG)\d+) .+$/\1/' ../../hybpiper_assembly/mb2/retrieved_single_copy_mb2/${g}.FNA | pxrls -c ../sample_code_list.txt -n ../tip_label_list.txt > unaligned/${g}.fa
done < mb2_genes_without_paralogs.txt

#align and trim
conda activate alignment
cd ~/acacia/alignments

parallel -j 50 --progress bash scripts/align_and_trim.sh {} mb2 :::: mb2/mb2_genes_without_paralogs.txt

#remove empty samples removed by trimal files
while read g; do
test -s trimmed/removed_samples_${g}.txt || rm -f trimmed/removed_samples_${g}.txt
done < mb2_genes_without_paralogs.txt

#concatenate fasta files for supermatrix
sed 's/^/trimmed\//;s/$/_trimmed.fa/' mb2_genes_without_paralogs.txt > mb2_supermatrix_aln_file_list.txt
pxcat -f mb2_supermatrix_aln_file_list.txt -p mb2_partitions_unedited.txt > mb2_concatenated_alignments.fa
sed 's/trimmed\///; s/_trimmed.fa//' mb2_partitions_unedited.txt > mb2_partitions.txt

#copy alignment file and partitions file to phylogenies dir
cd ~/acacia/phylogenies/mb2
cp ../../alignments/mb2/mb2_concatenated_alignments.fa .
cp ../../alignments/mb2/mb2_partitions.txt .

#Iqtree find best paritioning scheme and models
test -d iqtree_supermatrix || mkdir iqtree_supermatrix
cd iqtree_supermatrix
nohup iqtree -s ../mb2_concatenated_alignments.fa -p ../mb2_partitions.txt -pre mb2_partition_finding \
-m MF+MERGE -T AUTO -ntmax 50 &> iqtree_mb2_partition_finding_nohup.log &
#should use -rclusterf 10 to speed up?

#use best partitioning scheme and models for tree search
nohup iqtree -s ../mb2_concatenated_alignments.fa -p mb2_partition_finding.best_scheme.nex -pre mb2_tree_search -o Dichrostachys_cinerea -B 1000 -alrt 1000 --boot-trees -T AUTO &> iqtree_mb2_tree_search_nohup.log &

#Estimate individual gene trees after finding the optimal model for each
cd ~/acacia/phylogenies/mb2
test -d iqtree_individual_gene_trees || mkdir iqtree_individual_gene_trees
cd iqtree_individual_gene_trees
#-o Dichrostachys_cinerea to specify an outgroup doesn't work because Dichrostachys cinerea does not have a sequence for each gene. it doesn't make sense to root on anything else if it is absent
OMP_NUM_THREADS=1 nohup iqtree -s ../mb2_concatenated_alignments.fa -S ../mb2_partitions.txt -pre mb2_indiv_gene_trees \
-B 1000 -alrt 1000 -T AUTO -ntmax 50 &> iqtree_mb2_indiv_gene_trees_nohup.log &

#extract aLRT and ufboot values from ML tree and put in separate tree files
perl -pe 's/([\d\.]+)\/([\d\.]+):/\1:/g' mb2_indiv_gene_trees.treefile > mb2_indiv_gene_trees_alrt.treefile
perl -pe 's/([\d\.]+)\/([\d\.]+):/\2:/g' mb2_indiv_gene_trees.treefile > mb2_indiv_gene_trees_ufbt.treefile


#Run ASTRAL
test -d astral || mkdir astral
cd astral
java -Xmx8G -jar ~/programs/Astral/astral.5.7.8.jar -i ../mb2_indiv_gene_trees_ufbt.treefile --outgroup Dichrostachys_cinerea

#add arbitrary branch lengths to tips so R can plot the astral tree
python ~/programs/Astral/add-bl.py - mb2_astral.tre > mb2_astral_bl.tre

#prune Senegalia chundra from astral tree
pxrmt -t mb2_astral_bl.tre -n S_chundra > mb2_astral_bl_pruned.tre

############
grep -A1 "$g" ../../hybpiper_assembly/mb2/mb2_targets.fa | perl -pe '$_ = ">Albizia_julibrissin\n" if /^>/' >> ${g}.fa
mafft --auto ${g}.fa | perl -pe 'tr/n/-/ if /^(?!>)/' > ${g}_mafft.fa
trimal -in ${g}_mafft.fa -out ${g}_aligned.fa -noallgaps
trimal -in ${g}_aligned.fa -out ${g}_trimmed.fa -automated1 

while read s; do
while read g; do
echo $s/$g/$s/exonerate_stats.tsv
done < mb2_genes_without_paralogs.txt
done < <(head -1 ../full_sample_list.txt)
############


### Chloroplast genes ###
# Get chloplast genes (cpg) target file
cd ~/acacia/hybpiper_assembly
mkdir cpg
sed -e 's/-/-cp_/' ../plastome_targets/plastome_targets.fa | sed '/^>/ s/$/ foo/g' | tr -d '\n' | sed 's/ foo/\n/g; s/>/\n>/2g' > cpg/cpg_targets.fa

# Run hybpiper assemble for chloroplast gene targets
cd cpg
nohup bash ../scripts/hybpiper_assemble.sh cpg_targets.fa &> ../logs/hybpiper_assemble_cpg_nohup.log &

# Run hybpiper assemble for chloroplast gene targets for outgroups
nohup bash ../scripts/hybpiper_assemble_outgroups.sh cpg_targets.fa &> ../logs/hybpiper_assemble_cpg_outgroups_nohup.log &

# Try HybPiper again for samples with only a few recovered genes
# with no read merging and --cov_cutoff 2
#make list of samples
echo -e "OG01\nOG04" > no_merge_sample_list.txt

nohup bash ../scripts/hybpiper_assemble_no_merge.sh cpg_targets.fa &> ../logs/hybpiper_assemble_cpg_no_merge_nohup.log &


# Collect stats from cpg target assemblies
targetFile=cpg_targets.fa
hybpiper stats -t_dna $targetFile gene --seq_lengths_filename cpg_exon_lengths --stats_filename cpg_exons_stats ../full_sample_list.txt
#hybpiper stats -t_dna $targetFile supercontig --seq_lengths_filename cpg_supercontig_lengths --stats_filename cpg_supercontigs_stats ../full_sample_list.txt

hybpiper recovery_heatmap --heatmap_filename cpg_exon_recovery_heatmap --heatmap_filetype pdf cpg_exon_lengths.tsv

# Collect stats on paralogs
hybpiper paralog_retriever -t_dna $targetFile \
--paralog_report_filename cpg_paralog_report --paralogs_above_threshold_report_filename cpg_paralogs_above_threshold_report \
--heatmap_filename cpg_paralog_heatmap --heatmap_filetype pdf ../full_sample_list.txt

#Using all genes even if they had paralog warnings
#make a list of cp genes
grep "^>" cpg_targets.fa | perl -pe 's/^>((Ssengal)|(Vnilotica))-//' | sort | uniq > cpg_genes.txt

#retrieve sequences for cp genes 
hybpiper retrieve_sequences -t_dna cpg_targets.fa --sample_names ../full_sample_list.txt --fasta_dir retrieved_cpg dna

# Alignment
#make parallel lists of sample codes and tip labels to change the sequence ids for alignments
#alread done
#cd ~/acacia/hybpiper_assembly
#cut -f2 ../read_processing/ssp01_acacia_lib_info.tsv | tail -n+2 > ../alignments/sample_code_list.txt
#cut -f2 ../read_processing_outgroups/mimosoid_outgroups_lib_info.tsv | tail -n+2 >> ../alignments/sample_code_list.txt
#cut -f7 ../read_processing/ssp01_acacia_lib_info.tsv | tail -n+2 > ../alignments/tip_label_list.txt
#cut -f7 ../read_processing_outgroups/mimosoid_outgroups_lib_info.tsv | tail -n+2 >> ../alignments/tip_label_list.txt

cd ~/acacia/alignments
test -d cpg || mkdir cpg
cd cpg
test -d unaligned || mkdir unaligned

#copy unaligned fastas and change seq ids from sample codes to taxon names
cp ../../hybpiper_assembly/cpg/cpg_genes.txt .
while read g; do
perl -pe 's/^(>(A|OG)\d+) .+$/\1/' ../../hybpiper_assembly/cpg/retrieved_cpg/${g}.FNA | pxrls -c ../sample_code_list.txt -n ../tip_label_list.txt > unaligned/${g}.fa
done < cpg_genes.txt

#align and trim
conda activate alignment
cd ~/acacia/alignments

#modify the Albizia cp cds file to have matching gene names
perl -pe 's/-(\w+)$/-cp_\1/' ../plastome_targets/Ajulibrissin_plastome_cds.fa | sed '/^>/ s/$/ foo/g' | tr -d '\n' | sed 's/ foo/\n/g; s/>/\n>/2g' > cpg/Albizia_cpg_genes.fa

parallel -j 50 --progress bash scripts/align_and_trim.sh {} cpg :::: cpg/cpg_genes.txt

#remove empty samples removed by trimal files
while read g; do
test -s trimmed/removed_samples_${g}.txt || rm -f trimmed/removed_samples_${g}.txt
done < cpg_genes.txt

#concatenate fasta files for supermatrix
sed 's/^/trimmed\//;s/$/_trimmed.fa/' cpg_genes.txt > cpg_supermatrix_aln_file_list.txt
pxcat -f cpg_supermatrix_aln_file_list.txt -p cpg_partitions_unedited.txt > cpg_concatenated_alignments.fa
sed 's/trimmed\///; s/_trimmed.fa//' cpg_partitions_unedited.txt > cpg_partitions.txt

# Estimate Phylogenies
#copy alignment file and partitions file to phylogenies dir
cd ~/acacia/phylogenies/cpg
cp ../../alignments/cpg/cpg_concatenated_alignments.fa .
cp ../../alignments/cpg/cpg_partitions.txt .

#Iqtree find best partitioning scheme and models
test -d iqtree_supermatrix || mkdir iqtree_supermatrix
cd iqtree_supermatrix
nohup iqtree -s ../cpg_concatenated_alignments.fa -p ../cpg_partitions.txt -pre cpg_partition_finding \
-m MF+MERGE -T 4 &> iqtree_cpg_partition_finding_nohup.log &

#use best partitioning scheme and models for tree search
nohup iqtree -s ../cpg_concatenated_alignments.fa -p cpg_partition_finding.best_scheme.nex -pre cpg_tree_search -o Dichrostachys_cinerea -B 1000 -alrt 1000 --boot-trees -T AUTO &> iqtree_cpg_tree_search_nohup.log &

#extract aLRT and ufboot values from ML tree and put in separate tree files
perl -pe 's/([\d\.]+)\/([\d\.]+):/\1:/g' cpg_tree_search.treefile > cpg_tree_search_alrt.treefile
perl -pe 's/([\d\.]+)\/([\d\.]+):/\2:/g' cpg_tree_search.treefile > cpg_tree_search_ufbt.treefile

#get true majority rule consensus tree and estimate branch lengths with same model used to find ML tree while collapsing near zero-length branches
iqtree -t cpg_tree_search.ufboot -con -minsup 0.5 -pre cpg_tree_search_50mrc -o Dichrostachys_cinerea -T AUTO -ntmax 4
iqtree -s ../cpg_concatenated_alignments.fa -p cpg_partition_finding.best_scheme.nex -te cpg_tree_search_50mrc.contree \
-pre cpg_tree_search_50mrc_collapsed -czb -o Dichrostachys_cinerea -T AUTO -ntmax 4

#prune Stryphnodendron from cpg tree
pxrmt -t cpg_tree_search_ufbt.treefile -n Stryphnodendron_adstringens > cpg_tree_search_ufbt_pruned.tre
pxrmt -t cpg_tree_search_50mrc_collapsed.treefile -n Stryphnodendron_adstringens > cpg_tree_search_50mrc_pruned.tre

##############
sampleCode=A004
targetFile=mb2_targets.fa
targetType=${targetFile%_targets.fa}
pairedReads="../../read_processing/$sampleCode/${sampleCode}.clean_[12].fq.gz"
singleReads="../../read_processing/$sampleCode/${sampleCode}.clean_S.fq.gz"
hybpiper assemble --cpu 50 --bwa -t_dna $targetFile --readfiles $pairedReads --unpaired $singleReads --cov_cutoff 2 --prefix $sampleCode --keep_intermediate_files

hybpiper assemble --cpu 50 --bwa -t_dna $targetFile --readfiles $pairedReads --unpaired $singleReads --cov_cutoff 2 --prefix $sampleCode --keep_intermediate_files --start_from exonerate_contigs


exonerate -m protein2genome --showalignment yes --showvulgar no -V 0 --refine full --ryo ">%ti,%qi,%qab,%qae,%pi,(%tS),%tab,%tae\\n%tcs\\n" \
mb2_gene128_target.fasta mb2_gene128_contigs_k31.fasta > exonerate_results.fasta

hybpiper assemble --cpu 50 --bwa -t_dna $targetFile --readfiles $pairedReads --unpaired $singleReads --target Albizia --prefix $sampleCode


