### NUCLEAR PROBES COMPARISON ###

mkdir a353_compared
wget https://github.com/oscarvargash/lecy_taxo/blob/main/trimming_and_assembly/Lecythid_complete_targets
#comparing with Arapidopsis thaliana genes
#getting the Annotation_Summary.txt from the suplementary materials (Jhonson et al., 2019)
cut -f2 Annotation_Summary.txt | grep "^AT" | sort > a353_TAIR_ID.txt
#getting the seq IDs from the Lecy probe (Vargas et al., 2019)
grep "^>" Lecythid_complete_targets.txt | perl -ple 's/>DAT\d+-//' | sort | uniq > Lecy_354_TAIR_ID.txt
#common genes for both datasets:
comm -12 a353_TAIR_ID.txt Lecy_354_TAIR_ID.txt > shared_TAIR_ID.txt
