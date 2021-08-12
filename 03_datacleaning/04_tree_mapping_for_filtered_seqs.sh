#!bin/bash


#. ~/miniconda3/etc/profile.d/conda.sh
. activate qiime2-2020.11

# make wd
output='04_tree_mapping_for_filtered_seqs'
mkdir $output
### Paths #####

#fna='../02_externally_generated_data/16S_sequencing/intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito.qza'
alexFilteredSeqs='../02_externally_generated_data/16S_sequencing/intermediate_files/imported_databases/alexFilteredSeqs.qza'
mySeqs='../02_externally_generated_data/16S_sequencing/intermediate_files/dada2/rep-seqs.qza'
seqkeep='03_manual_seq_compare/downstream/subsetASVs_for_tree_list.txt'
taxonomy='../02_externally_generated_data/16S_sequencing/intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza'
#rep_align='../02_externally_generated_data/16S_sequencing/raw_data/SILVA_128_QIIME_release/rep_set_aligned/99/99_otus_aligned.fasta'

### Filter sequencess
qiime feature-table filter-seqs \
--i-data $mySeqs \
--m-metadata-file $seqkeep \
--o-filtered-data $output/mySeqs_filtered_seqs.qza

### Merge mine and alex's sequences
qiime feature-table merge-seqs \
--i-data $output/mySeqs_filtered_seqs.qza \
--i-data $alexFilteredSeqs \
--o-merged-data $output/mergedwAlex_filteredseqs.qza


## Align sequences with MAFFT
# REFERNECE: https://doi.org/10.1093/molbev/mst010

qiime alignment mafft \
--i-sequences $output/mergedwAlex_filteredseqs.qza \
--o-alignment $output/mergedwAlex_aligned.qza

# mask alignmnet 
qiime alignment mask \
--i-alignment $output/mergedwAlex_aligned.qza \
--o-masked-alignment  $output/mergedwAlex_aligned_masked.qza


## Make RAxML tree
# REFERENCE: https://doi.org/10.1093/bioinformatics/btu033
# Rapid bootsrap ref: https://dx.doi.org/10.1080/10635150802429642

qiime phylogeny raxml-rapid-bootstrap \
--i-alignment $output/mergedwAlex_aligned_masked.qza \
--p-seed 1239 \
--p-rapid-bootstrap-seed 8394 \
--p-bootstrap-replicates 100 \
--p-substitution-model GTRCAT \
--o-tree $output/raxml-cat-bootstrap-tree.qza \
--verbose

### Export tree
qiime tools export \
--input-path $output/raxml-cat-bootstrap-tree.qza \
--output-path $output/raxml_tree_filteredseqs.tre




