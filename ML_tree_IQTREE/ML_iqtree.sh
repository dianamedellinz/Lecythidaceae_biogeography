#Getting gene trees from the all_loci_ML_tree tree file (664 genes)

grep "344_AT" Parts.txt > Parts_344.txt
grep "a353_" Parts.txt > Parts_a353.txt

#merging exons and introns for 344 genes
while read g; do grep "$g" ../Parts_344.txt | sed -e 's/exon\///; s/intron_ok\///' > Part_344_${g}.txt; done < 344_unique_genes_list.txt

#MFP Extended model selection followed by gene tree
parallel -j 12 --progress "OMP_NUM_THREADS=1 iqtree2 -s ../Supermatrix.fa -p Part_344_{}.txt -pre 344_{}_tree_search -m MFP -B 1000 -alrt 1000 -T 2 &> 344_{}_tree_search_parallel.log" :::: 344_unique_genes_list.txt

# Appending the treefiles from the 344 probe set to a single file:
while read g; do cat 344_${g}_tree_search.treefile >> all_344_tree_search.treefile; done < 344_unique_genes_list.txt

#extract aLRT and ufboot values from ML trees and put them in separate tree files:
perl -pe 's/([\d\.]+)\/([\d\.]+):/\1:/g' all_344_tree_search.treefile > all_344_tree_search_alrt.treefile
perl -pe 's/([\d\.]+)\/([\d\.]+):/\2:/g' all_344_tree_search.treefile > all_344_tree_search_ufbt.treefile

#MFP Extended model selection followed by gene tree for a353
nohup parallel -j 12 "OMP_NUM_THREADS=1 iqtree2 -s ../Supermatrix.fa -p Part_a353_{}.txt -pre a353_{}_tree_search -m MFP -B 1000 -alrt 1000 -T 2 &> a353_{}_tree_search_parallel.log" :::: a353_unique_genes_list.txt &> a353_parallel_nohup_output.log &

# Appending the treefiles from a353 to a single file:
while read g; do cat a353_${g}_tree_search.treefile >> all_a353_tree_search.treefile; done < a353_unique_genes_list.txt

#extract aLRT and ufboot values from ML trees and put them in separate tree files:
perl -pe 's/([\d\.]+)\/([\d\.]+):/\1:/g' all_a353_tree_search.treefile > all_a353_tree_search_alrt.treefile
perl -pe 's/([\d\.]+)\/([\d\.]+):/\2:/g' all_a353_tree_search.treefile > all_a353_tree_search_ufbt.treefile

#copying and appending the two treefiles (31_overlapped_genes, a353 because these have Marcgraviaceae/Tetrameristedaceae as the outgroup):
cp ~/diana/Lecy_2025/03_iqtree/a353_gene_trees/all_a353_tree_search_ufbt.treefile ~/diana/Lecy_2025/04_astral
cp ~/diana/Lecy_2025/03_iqtree/31_overlapped_gene_trees/all_31_overlapped_tree_search_ufbt.treefile ~/diana/Lecy_2025/04_astral

# Appending the three treefiles:
cat all_31_overlapped_tree_search_ufbt.treefile all_a353_tree_search_ufbt.treefile > 353_plus_31_loci_gene_trees_ufbt.treefile
cat 31_overlapped_genes_list.txt a353_unique_genes_list.txt > all_353_loci_genes_list.txt

#rerooting the 353 gene trees using R script
cp ~/diana/a353_combined/06_astral/phyparts/reroot_genetrees.R .
Rscript reroot_genetrees.R 353_plus_31_loci_gene_trees_ufbt.treefile  outgroup_species_list.txt all_loci_ .lecy.treefile

#rerooting the 344 gene trees using R script
cp ~/diana/a353_combined/06_astral/phyparts/reroot_genetrees.R .
Rscript reroot_genetrees.R all_344_tree_search_ufbt.treefile 344_unique_genes_list.txt outgroup_samples_344_genes.txt all_344_ .lecy.treefile
