## ALIGNMENT AND TRIMMING - Jan 2025

# 1.a. running alignment for the 31_overlapped genes EXONS
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/exon/{}_31_overlapped_exon.fasta > exon/aligned_{}_31_overlapped_exon.fa" :::: <(cut -f1 ../00_unaligned/31_genes_ID.txt)

# 1.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g exon; done < <(cut -f1 ../00_unaligned/31_genes_ID.txt)

# 2.a. running alignment for a353 genes exons
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/exon/a353_{}_exon.fa > exon/aligned_a353_{}_exon.fa" :::: ../00_unaligned/a353_unique_genes_list.txt

# 2.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g exon; done < ../00_unaligned/a353_unique_genes_list.txt

# 3.a. running alignment for 344 genes exons
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/exon/{}_exon.fasta > exon/aligned_344_{}_exon.fa" :::: ../00_unaligned/344_unique_genes_list.txt

# 3.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g exon; done < ../00_unaligned/344_unique_genes_list.txt

# 4.a. running alignment for the 31_overlapped genes INTRONS
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/intron/{}_31_overlapped_intron.fasta > intron/aligned_{}_31_overlapped_intron.fa" :::: <(cut -f1 ../00_unaligned/31_genes_ID.txt)

# 4.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g intron 31_overlapped; done < <(cut -f1 ../00_unaligned/31_genes_ID.txt)

# 5.a. running alignment for a353 genes introns
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/intron/{}_intron.fa > intron/aligned_a353_{}_intron.fa" :::: ../00_unaligned/a353_unique_genes_list.txt

# 5.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g intron a353; done < ../00_unaligned/a353_unique_genes_list.txt

# 6.a. running alignment for 344 genes introns
cd ~/diana/Lecy_2025/01_alignment/
parallel -j 20 --progress "mafft --thread 2 --auto ../00_unaligned/intron/{}_intron.fasta > intron/aligned_344_{}_intron.fa" :::: ../00_unaligned/344_unique_genes_list.txt

# 6.b. running trimal
conda activate trimming
cd ~/diana/Lecy_2025/02_trimal/
while read g; do bash trimal.sh $g intron 344; done < ../00_unaligned/344_unique_genes_list.txt