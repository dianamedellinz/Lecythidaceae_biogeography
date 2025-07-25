#Running ASTRAL for the 664 gene trees

java -Xmx8G -jar /home/brlab/Astral/astral.5.7.8.jar -i all_loci_gene_trees_ufbt.treefile -o Lecy_all_loci_astral.tre 2> Lecy_all_loci_astral.out.log

#colapsing branches with bootstrap support < 10%, then run Astral again:cd

#nw are utilities that fix things like nw_prune, nw_reroot, nw_ed (editing), etc.
conda activate nwUtils
nw_ed all_loci_gene_trees_ufbt.treefile 'i & b<=10' o > all_loci_gene_trees_ufbt-BS10.treefile

#dropping samples <2000bp from genes trees
pxrmt -s -t all_loci_gene_trees_ufbt-BS10.treefile -f samples_to_drop_gene_trees > all_loci_gene_trees_ufbt-BS10_clean.treefile

java -Xmx8G -jar /home/brlab/Astral/astral.5.7.8.jar -i all_loci_gene_trees_ufbt-BS10_clean.treefile -o Lecy_all_loci_astral-BS10_clean.tre 2> Lecy_all_loci_astral-BS10_clean.out.log

#rerooting the gene trees using R script
Rscript reroot_genetrees.R all_loci_gene_trees_ufbt-BS10_clean.treefile 665_genes_list.txt outgroup_species_list.txt all_loci_ .lecy.treefile

#rooting the clean astral tree
Rscript reroot_genetrees.R Lecy_all_loci_astral-BS10_clean.tre rooted_astral_tree_filename.txt outgroup_species_list.txt
sed -e 's/Root;$/;/' -e 's/NaN/1.0/g' Lecy_all_loci_astral-BS10_clean_rooted.tre > Lecy_all_loci_astral-BS10_clean_rooted.treefile

#dropping samples <2000bp from astral tree
pxrmt -s -t Lecy_all_loci_astral-BS10_clean_rooted.treefile -f samples_to_drop_gene_trees > Lecy_all_loci_astral-BS10_clean_rooted.tre
