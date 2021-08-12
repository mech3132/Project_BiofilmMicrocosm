#!bin/bash

#. ~/miniconda3/etc/profile.d/conda.sh
. activate qiime2-2020.11

### Paths #####
mappingfile='raw_data/mapping_files/sample_prep_info.txt'
alexseqs='raw_data/alex_isolate_sequences/Full_BT_bact.fasta'
inhibIso='../../00_isolate_info/manual_chosen_isolates/Inhib_to_include.txt'
nonIso='../../00_isolate_info/manual_chosen_isolates/nonInhib_to_include.txt'
###### Step one: demultiplex ###########

#mkdir paired-end-sequences
#cp raw_data/seq_data/Undetermined* paired-end-sequences/

#mv paired-end-sequences/*R1* paired-end-sequences/forward.fastq.gz
#mv paired-end-sequences/*R2* paired-end-sequences/reverse.fastq.gz
#mv paired-end-sequences/*I1* paired-end-sequences/barcodes.fastq.gz

## Import using qiime tools
#qiime tools import \
#--type EMPPairedEndSequences \
#--input-path paired-end-sequences \
#--output-path paired-end-sequences/emp-paired-end-sequences.qza

## Demultiplex
#qiime demux emp-paired \
#--m-barcodes-file $mappingfile \
#--m-barcodes-column barcode \
#--p-rev-comp-mapping-barcodes \
#--i-seqs paired-end-sequences/emp-paired-end-sequences.qza \
#--o-per-sample-sequences paired-end-sequences/demux-full.qza \
#--o-error-correction-details paired-end-sequences/demux-details.qza

## Look at  summary
#qiime demux summarize \
#--i-data paired-end-sequences/demux-full.qza \
#--o-visualization paired-end-sequences/demux-full.qzv
#qiime tools export \
#--input-path paired-end-sequences/demux-full.qzv \
#--output-path paired-end-sequences/demux-full-vis/

## Clean up and remove intermediate files
#rm paired-end-sequences/*fastq.gz

#### TESTE WITH DIFFERENT DEMUX PARAMS DUE TO LOW READ COUNT ###

mkdir paired-end-sequences-2
#cp raw_data/seq_data/*McKenzie* paired-end-sequences-2/
cp raw_data/seq_data/Undetermined* paired-end-sequences/

mv paired-end-sequences-2/*R1* paired-end-sequences-2/forward.fastq.gz
mv paired-end-sequences-2/*R2* paired-end-sequences-2/reverse.fastq.gz
mv paired-end-sequences-2/*I1* paired-end-sequences-2/barcodes.fastq.gz

# Import using qiime tools
qiime tools import \
--type EMPPairedEndSequences \
--input-path paired-end-sequences-2 \
--output-path paired-end-sequences-2/emp-paired-end-sequences.qza

# Demultiplex
qiime demux emp-paired \
--m-barcodes-file $mappingfile \
--m-barcodes-column barcode \
--i-seqs paired-end-sequences-2/emp-paired-end-sequences.qza \
--p-rev-comp-mapping-barcodes \
--p-rev-comp-barcodes \
--o-per-sample-sequences paired-end-sequences-2/demux-full.qza \
--o-error-correction-details paired-end-sequences-2/demux-details.qza

# Look at  summary
qiime demux summarize \
--i-data paired-end-sequences-2/demux-full.qza \
--o-visualization paired-end-sequences-2/demux-full.qzv
qiime tools export \
--input-path paired-end-sequences-2/demux-full.qzv \
--output-path paired-end-sequences-2/demux-full-vis/

# For sequence submission
qiime tools export \
--input-path paired-end-sequences-2/demux-full.qza \
--output-path paired-end-sequences-2/exported_sequences

# Clean up and remove intermediate files
rm paired-end-sequences-mckenzie/*fastq.gz
################## 

mkdir intermediate_files
mkdir intermediate_files/dada2
# denoise with dada2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-sequences-2/demux-full.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table intermediate_files/dada2/table.qza \
  --o-representative-sequences intermediate_files/dada2/rep-seqs.qza \
  --o-denoising-stats intermediate_files/dada2/denoising-stats.qza
  
  
  ## Notes on DADA2:
  # dada2 does quality filtering within itself, and also checks for chimeras, does abundance filtering, etc. OUtput is PURE ASVs babyyy


qiime feature-table summarize \
  --i-table intermediate_files/dada2/table.qza \
  --o-visualization intermediate_files/dada2/table.qzv \
  --m-sample-metadata-file $mappingfile
qiime tools export \
  --input-path intermediate_files/dada2/table.qzv \
  --output-path intermediate_files/dada2/table-vis
  
#### Preparing databases ####

mkdir intermediate_files/imported_databases

# import gg database, silva database, alex, and woodhams
# GG
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ./raw_data/gg_13_8_otus/rep_set/99_otus.fasta \
--output-path ./intermediate_files/imported_databases/gg_99_repset.qza
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ./raw_data/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
--output-path ./intermediate_files/imported_databases/99_gg_ref-taxonomy.qza
# SILVA
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ./raw_data/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta \
--output-path ./intermediate_files/imported_databases/SILVA_99_repset.qza
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ./raw_data/SILVA_128_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt \
--output-path ./intermediate_files/imported_databases/99_SILVA_ref-taxonomy.qza

# Import alex
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $alexseqs \
--output-path ./intermediate_files/imported_databases/alex_imported.qza
# Import woodhams
# There's a flaw in the woodhams database; : signs. Remove all the ones with this'
tac ./raw_data/woodhams_database/Amphibian-skin_bacteria_16S_sequences.fna | sed '/:/I,+1 d' | tac > ./intermediate_files/imported_databases/woodhams_temp.fna 
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ./intermediate_files/imported_databases/woodhams_temp.fna  \
--output-path ./intermediate_files/imported_databases/woodhams_16S.qza

## Get relavent sequences from Alex
cat $inhibIso $nonIso | sort -u > ./intermediate_files/allIsolates.txt
# now, get sequence IDs
echo '#SampleID' > ./intermediate_files/imported_databases/alexKeepSames.txt
while read id; do
	grep $id $alexseqs | sed 's/>//g' >> ./intermediate_files/imported_databases/alexKeepSames.txt
done < ./intermediate_files/allIsolates.txt
sort ./intermediate_files/imported_databases/alexKeepSames.txt -u |sort -r > ./intermediate_files/imported_databases/alexKeepSamesunique.txt
rm ./intermediate_files/imported_databases/alexKeepSames.txt

# Filter alex's non-relative sequences
qiime feature-table filter-seqs \
--i-data ./intermediate_files/imported_databases/alex_imported.qza \
--m-metadata-file ./intermediate_files/imported_databases/alexKeepSamesunique.txt \
--o-filtered-data ./intermediate_files/imported_databases/alexFilteredSeqs.qza


# extract reads from gg, SILVA, Alex, and woodhams
qiime feature-classifier extract-reads \
--i-sequences ./intermediate_files/imported_databases/gg_99_repset.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-min-length 150 \
--p-max-length 400 \
--o-reads ./intermediate_files/imported_databases/gg-99-repset_extractedreads515806.qza
#qiime feature-classifier extract-reads \
#--i-sequences ./intermediate_files/imported_databases/SILVA_99_repset.qza \
#--p-f-primer GTGCCAGCMGCCGCGGTAA \
#--p-r-primer GGACTACHVGGGTWTCTAAT \
#--p-min-length 150 \
#--p-max-length 400 \
#--o-reads ./intermediate_files/imported_databases/SILVA-99-repset_extractedreads515806.qza
qiime feature-classifier extract-reads \
--i-sequences ./intermediate_files/imported_databases/woodhams_16S.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-min-length 150 \
--p-max-length 400 \
--o-reads ./intermediate_files/imported_databases/woodhams_16S_extractedreads515806.qza
qiime feature-classifier extract-reads \
--i-sequences ./intermediate_files/imported_databases/alexFilteredSeqs.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-min-length 150 \
--p-max-length 400 \
--o-reads ./intermediate_files/imported_databases/alex_extractedreads515806.qza

# Merge all reads together
mkdir intermediate_files/merged_repset
qiime feature-table merge-seqs \
--i-data ./intermediate_files/imported_databases/alex_extractedreads515806.qza ./intermediate_files/dada2/rep-seqs.qza ./intermediate_files/imported_databases/woodhams_16S_extractedreads515806.qza \
--o-merged-data ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged.qza

##### Assigning taxonomy ######

mkdir intermediate_files/bayes_classifier
# Build bayes classifer for gg
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ./intermediate_files/imported_databases/gg-99-repset_extractedreads515806.qza \
--i-reference-taxonomy ./intermediate_files/imported_databases/99_gg_ref-taxonomy.qza \
--o-classifier ./intermediate_files/bayes_classifier/99_gg_13_8_bayesclassifier.qza

qiime feature-classifier classify-sklearn \
--i-reads ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged.qza \
--i-classifier ./intermediate_files/bayes_classifier/99_gg_13_8_bayesclassifier.qza \
--o-classification ./intermediate_files/bayes_classifier/gg_classifications_biofilm_alex_woodhams_merged.qza

# Export
qiime tools export \
--input-path ./intermediate_files/bayes_classifier/gg_classifications_biofilm_alex_woodhams_merged.qza \
--output-path ./intermediate_files/bayes_classifier/gg_classifications_biofilm_alex_woodhams_merged

# The bayes classifer for SILVA is too memory-intensive to build-- so I'm going to use the existing 138 classifier on the QIIME2 resources page, that uses 515/806 sequences

qiime feature-classifier classify-sklearn \
--i-reads ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged.qza \
--p-reads-per-batch 500 \
--i-classifier ./raw_data/SILVA_138_classifier/silva-138-99-515-806-nb-classifier.qza \
--o-classification ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza

# Filter out mitochondria, chloroplast
mkdir ./intermediate_files/filtering_table
qiime taxa filter-table \
--i-table intermediate_files/dada2/table.qza \
--i-taxonomy ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza \
--p-exclude mitochondria,chloroplast,archaea \
--o-filtered-table ./intermediate_files/filtering_table/00_table_withoutchloromito.qza

qiime taxa filter-seqs \
--i-sequences ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged.qza \
--i-taxonomy ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza \
--p-exclude mitochondria,chloroplast,archaea \
--o-filtered-sequences ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito.qza

# Export
qiime tools export \
--input-path ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza \
--output-path ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged

# Convert out seq repset file, for blastn later
qiime tools export \
--input-path ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito.qza \
--output-path ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito

## Quickly checking to see how many sequences were lost after filter chloro, mito, archaea
qiime feature-table summarize --i-table ./intermediate_files/filtering_table/00_table_withoutchloromito.qza \
--o-visualization ./intermediate_files/filtering_table/00_table_withoutchloromito.qzv

## Redone after finding out demultipexing wasn't work:
# Original filter table had 5,939,614 reads
# filtered table ahd 5,938,950 reads

## The filtered table has 221 features; 118,137 reads
## The original table had 230 features; 118,234 reads
## So we lost about 100 reads; that's fine.

## Build bayes classifer for SILVA
#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads ./intermediate_files/imported_databases/SILVA-99-repset_extractedreads515806.qza \
#--i-reference-taxonomy ./intermediate_files/imported_databases/99_SILVA_ref-taxonomy.qza \
#--o-classifier ./intermediate_files/bayes_classifier/99_SILVA128_bayesclassifier.qza

#qiime feature-classifier classify-sklearn \
#--i-reads ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged.qza \
#--i-classifier ./intermediate_files/bayes_classifier/99_SILVA128_bayesclassifier.qza \
#--o-classification ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza
## Export
#qiime tools export \
#--input-path ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged.qza \
#--output-path ./intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged

## Citation
#Bokulich, N.A., Kaehler, B.D., Rideout, J.R. et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome 6, 90 (2018). https://doi.org/10.1186/s40168-018-0470-z

##### Phylogenetic tree creation ######
mkdir intermediate_files
mkdir intermediate_files/phylogeny

# Merge sequences from mine, alex, and woodhams so that I can generate a "master tree"

### NOTE: apparently short frags are not good for making a tree, so I'm using insert fragments into the gg tree
#qiime phylogeny align-to-tree-mafft-fasttree \
#--i-sequences intermediate_files/dada2/rep-seqs.qza \
#--o-alignment aligned_rep_seqs.qza \
#--o-masked-alighment masked-alignment-rep-seqs.qza \
#--o-tree unrooted-tree.qza \
#--o-rooted-tree rooted-tree.qza \

# Import SILVA 128 database 
qiime tools import \
--input-path ./raw_data/sepp-package-silva/ref \
--type SeppReferenceDatabase \
--output-path ./intermediate_files/imported_databases/99_otus_aligned_masked1977_SILVASEPP
## CITATION: smirarab/sepp-refs on github

# Insert all fragments into a reference phylogeny using SEPP
#qiime fragment-insertion sepp \
#--i-representative-sequences ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito.qza \
#--i-reference-database ./raw_data/sepp-refs-gg-13-8.qza \
#--o-tree ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_gg.qza \
#--o-placements ./intermediate_files/phylogeny/SEPP_placements_gg.qza

qiime fragment-insertion sepp \
--i-representative-sequences ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito.qza \
--i-reference-database ./intermediate_files/imported_databases/99_otus_aligned_masked1977_SILVASEPP.qza \
--o-tree ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_SILVA.qza \
--o-placements ./intermediate_files/phylogeny/SEPP_placements_SILVA.qza

# Remove features that didn't fit into the tree from my database
## gg
#qiime fragment-insertion filter-features \
#--i-table ./intermediate_files/filtering_table/00_table_withoutchloromito.qza \
#--i-tree ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_gg.qza \
#--o-filtered-table ./intermediate_files/filtering_table/01_SEPP_filtered_table_gg.qza \
#--o-removed-table ./intermediate_files/filtering_table/01_SEPP_removed_table_gg.qza 
## Check if placements are zero
#qiime feature-table summarize \
#--i-table ./intermediate_files/filtering_table/01_SEPP_removed_table_gg.qza \
#--o-visualization ./intermediate_files/filtering_table/01_SEPP_removed_table_gg.qzv
# There are actually 0 features in the removed table; all ASVs kept
## Export tree and placements
#qiime tools export \
#--input-path ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_gg.qza \
#--output-path ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_gg
#qiime tools export \
#--input-path ./intermediate_files/phylogeny/SEPP_placements_gg.qza \
#--output-path ./intermediate_files/phylogeny/SEPP_placements_gg

## SILVA
qiime fragment-insertion filter-features \
--i-table ./intermediate_files/filtering_table/00_table_withoutchloromito.qza \
--i-tree ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_SILVA.qza \
--o-filtered-table ./intermediate_files/filtering_table/01_SEPP_filtered_table_SILVA.qza \
--o-removed-table ./intermediate_files/filtering_table/01_SEPP_removed_table_SILVA.qza 
# Check if placements are zero
qiime feature-table summarize \
--i-table ./intermediate_files/filtering_table/01_SEPP_removed_table_SILVA.qza \
--o-visualization ./intermediate_files/filtering_table/01_SEPP_removed_table_SILVA.qzv
# There are actually 0 features in the removed table; all ASVs kept
# Export tree and placements
qiime tools export \
--input-path ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_SILVA.qza \
--output-path ./intermediate_files/phylogeny/SEPP_inserted_tree_rooted_SILVA
qiime tools export \
--input-path ./intermediate_files/phylogeny/SEPP_placements_SILVA.qza \
--output-path ./intermediate_files/phylogeny/SEPP_placements_SILVA

# Citation
#Phylogenetic Placement of Exact Amplicon Sequences Improves Associations with Clinical Information. Stefan Janssen, Daniel McDonald, Antonio Gonzalez, Jose A. Navas-Molina, Lingjing Jiang, Zhenjiang Zech Xu, Kevin Winker, Deborah M. Kado, Eric Orwoll, Mark Manary, Siavash Mirarab, Rob Knight. mSystems 2018. doi: https://doi.org/10.1128/mSystems.00021-18

##### OTU filtering ######
# Remove OTUs with less than 50 total reads in whole dataset and less than 10 reads per sample
#qiime feature-table filter-features \
#--i-table ./intermediate_files/filtering_table/01_SEPP_filtered_table.qza \
#--p-min-frequency 50 \
#--p-min-samples 2 \
#--o-filtered-table ./intermediate_files/filtering_table/02_filtered_otu_table_LT50_LS10.qza

#### Exporting OTU tables ######
#qiime tools export \
#--input-path ./intermediate_files/filtering_table/01_SEPP_filtered_table_gg.qza \
#--output-path ./intermediate_files/filtering_table/01_SEPP_filtered_table_gg
qiime tools export \
--input-path ./intermediate_files/filtering_table/01_SEPP_filtered_table_SILVA.qza \
--output-path ./intermediate_files/filtering_table/01_SEPP_filtered_table_SILVA
#qiime tools export \
#--input-path ./intermediate_files/filtering_table/02_filtered_otu_table_LT50_LS10.qza \
#--output-path ./intermediate_files/filtering_table/02_filtered_otu_table_LT50_LS10

# Convert to normal ASV table
#biom convert -i ./intermediate_files/filtering_table/01_SEPP_filtered_table_gg/feature-table.biom \
#--to-tsv \
#-o ./intermediate_files/filtering_table/01_otu_table_SEPPfiltered_gg.txt
biom convert -i ./intermediate_files/filtering_table/01_SEPP_filtered_table_SILVA/feature-table.biom \
--to-tsv \
-o ./intermediate_files/filtering_table/01_otu_table_SEPPfiltered_SILVA.txt
#biom convert -i ./intermediate_files/filtering_table/02_filtered_otu_table_LT50_LS10/feature-table.biom \
#--to-tsv \
#-o ./intermediate_files/filtering_table/02_otu_table_lowabundfiltered.txt

# Let's try clustering against Woodham's database at 100% similarity-- use "unfiltered" version of woodhams
# using blastn instead of vsearch I think? I need queryID

#qiime vsearch cluster-features-closed-reference \
#--i-sequences ./intermediate_files/dada2/rep-seqs.qza \
#--i-table ./intermediate_files/dada2/table.qza \
#--i-reference-sequences ./intermediate_files/imported_databases/woodhams_16S.qza \
#--p-perc-identity 1 \
#--output-dir ./intermediate_files/vsearch_against_woodhams

## Look at matched sequences
#qiime tools export \
#--input-table 


# Blastn search 
mkdir intermediate_files/blast
echo "Creating Woodhams db"
makeblastdb -in ./intermediate_files/imported_databases/woodhams_temp.fna -parse_seqids -blastdb_version 5  -title "Woodhams" -dbtype nucl -out ./intermediate_files/blast/Woodhams

echo "Creating Alex db"
qiime tools export \
--input-path ./intermediate_files/imported_databases/alexFilteredSeqs.qza \
--output-path ./intermediate_files/imported_databases/alexFilteredSeqs

makeblastdb -in ./intermediate_files/imported_databases/alexFilteredSeqs/dna-sequences.fasta -parse_seqids -blastdb_version 5  -title "AlexIsos" -dbtype nucl -out ./intermediate_files/blast/AlexIsos

# Blast against woodham's full db
blastn -db ./intermediate_files/blast/Woodhams -query ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito/dna-sequences.fasta -out ./intermediate_files/blast/blast_against_Woodhams.txt -outfmt "6 qseqid sseqid pident length mismatch stitle evalue" -perc_identity 97

blastn -db ./intermediate_files/blast/AlexIsos -query ./intermediate_files/merged_repset/biofilm_alex_woodhams_repset_merged_nochlormito/dna-sequences.fasta -out ./intermediate_files/blast/blast_against_Alex.txt -outfmt "6 qseqid sseqid pident length mismatch stitle evalue" -perc_identity 97

echo -e 'QueryID\tMatchedSequenceID\tPercentIdentity\tLength\tMismatch\ttaxaID\tevalue' > ./intermediate_files/headers_for_blast.txt
sed 's/-e //g' intermediate_files/headers_for_blast.txt > ./intermediate_files/headers_for_blast2.txt

cat ./intermediate_files/headers_for_blast2.txt ./intermediate_files/blast/blast_against_Woodhams.txt > intermediate_files/blast/blast_against_Woodhams_wheader.txt
cat ./intermediate_files/headers_for_blast2.txt ./intermediate_files/blast/blast_against_Alex.txt > intermediate_files/blast/blast_against_Alex_wheader.txt

rm intermediate_files/headers_for_blast*
rm intermediate_files/blast/blast_against_Woodhams.txt
rm intermediate_files/blast/blast_against_Alex.txt



