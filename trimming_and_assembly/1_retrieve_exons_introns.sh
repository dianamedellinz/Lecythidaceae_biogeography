# Collect stats from a353 target assemblies
targetFile=lecy_a353_ALL_targets.fa
hybpiper stats -t_dna $targetFile gene --seq_lengths_filename a353_exon_lengths --stats_filename a353_exons_stats sample_list.txt
hybpiper stats -t_dna $targetFile supercontig --seq_lengths_filename a353_supercontig_lengths --stats_filename a353_supercontigs_stats sample_list.txt

hybpiper recovery_heatmap --figure_length 60 --gene_text_size 2 --heatmap_filename a353_exon_recovery_heatmap --heatmap_filetype pdf a353_exon_lengths.tsv

# Collect stats on paralogs
hybpiper paralog_retriever -t_dna $targetFile --figure_length 60 --gene_text_size 4 \
--paralog_report_filename a353_paralog_report --paralogs_above_threshold_report_filename a353_paralogs_above_threshold_report.txt \
--heatmap_filename a353_paralog_heatmap --heatmap_filetype pdf sample_list.txt

### paralogs ###
##make a list of the genes (because there is more than one copy for each gene)
grep "^>" lecy_a353_ALL_targets.fa | perl -ple 's/^>\w+-//' | sort | uniq > a353_genes_list.txt

#python3 scripts/count_warnings_by_gene.py a353_genes_list.txt sample_list.txt

#retrieve exon sequences for genes without paralogs
hybpiper retrieve_sequences -t_dna lecy_a353_ALL_targets.fa --sample_names sample_list.txt --fasta_dir retrieved_exons_a353 dna

#retrieve introns
hybpiper retrieve_sequences -t_dna lecy_a353_ALL_targets.fa --cpu 20 --sample_names sample_list.txt --fasta_dir retrieved_introns_a353 intron
