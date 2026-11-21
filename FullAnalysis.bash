#!/bin/bash

# Run from the repository folder: bash FullAnalysis.bash
# Tools and R packages are in conda_environment.yml (conda activate zebra_finch)

###########################################
###########################################
##  Make data structure & download data  ##
###########################################
###########################################

# Download genomic data
mkdir Genomic_data
cd Genomic_data
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 "4A" "1A" "1B" "LG2" "LG5" "LGE22" "MT" "Z"
do
	wget ftp://ftp.ensembl.org/pub/release-98/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.dna.chromosome.$i.fa.gz
done
cd ..
pigz -d Genomic_data/*


# TF analysis
mkdir JASPAR_pfms_meme
cd JASPAR_pfms_meme
wget https://jaspar2018.genereg.net/download/data/2018/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.zip --user-agent "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.22 (KHTML, like Gecko) Ubuntu Chromium/25.0.1364.160 Chrome/25.0.1364.160 Safari/537.22"
unzip *
rm *.zip
cd ..

mkdir JASPAR_counts
mkdir JASPAR_output



##############
##############
##  RNAseq  ##
##############
##############


######################
# Age and Sex effect #
######################

# The RNA-seq count table (Plastoscine.txt, Plastoscine_annot.txt) is part of this repository.
# To start from the raw data instead: download the FASTQ files from GEO (GSE147974) and run
# bash fastq2counts.bash (conda environment zebra_finch_preprocessing, paths in environment.yml)


# Differential expression analysis
mkdir GeneLists
mkdir boxplots
Rscript DE_analysis.R

# General statistics & Plots for paper
Rscript DE_analysis_statisticsPaper.R
Rscript Circos.R

# Gene set analysis: done outside this repository (Genomatix GePS and Metascape) on the gene lists in GeneLists/

# TF analysis
mkdir fimo_RNAseq_out
Rscript ExtractSeq_TFsearching_RNAseq.R
N=40
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RNAseq_out/fimo_$file JASPAR_pfms_meme/$file seq_RNAseq_TFsearching.fa &
done
wait
)

Rscript Parse_FimoAll_RNAseq.R

for coi in "1dph" "20dph" "65dph" "adult" "age"
do
	Rscript JASPAR_TF_Pairwise.R $coi,nonsig &
done
wait


####################
# Cluster analysis #
####################

# Make clustering
mkdir Trends
Rscript RNAseq_ClusterAnalysis.R

# Gene set analysis: done outside this repository (Genomatix GePS and Metascape)

# TF analysis
for coi in $(ls GeneLists | grep patt | cut -d. -f 1 | cut -d_ -f 2,3,4,5)
do
	Rscript JASPAR_TF_Pairwise.R $coi &
done
wait



############
############
##  RRBS  ##
############
############

# RRBS pre-processing (Trim Galore, Bismark) is not part of this repository.
# The Bismark coverage files are available from GEO (GSE147974); set their folders in environment.yml

# Differential methylation analysis
Rscript RRBS_analyse.R
Rscript RRBS_analyse_hydroxy.R
Rscript RRBS_analyse_nonCpG.R

# General statistics & Plots for paper
Rscript RRBS_statisticsPaper.R
Rscript RRBS_plots.R

# Gene set analysis: done outside this repository (Genomatix GePS and Metascape)

# TF analysis
mkdir fimo_RRBS_out
Rscript ExtractSeq_TFsearching_RRBS.R
N=40
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RRBS_out/fimo_$file JASPAR_pfms_meme/$file seq_RRBS_TFsearching.fa &
done
wait
)

Rscript Parse_FimoAll_RRBS.R
Rscript JASPAR_TF_Pairwise_RRBS.R age

############################################
############################################
##  Correlation Methylation & Expression  ##
############################################
############################################


Rscript Expression_methylation_cor.R



#####################
#####################
##  Dosage Effect  ##
#####################
#####################

Rscript Sexually_dimorphic_genes.R

Rscript DosageEffect.R

# TF analysis
Rscript JASPAR_TF_Pairwise.R dosage,nodosage &
Rscript JASPAR_TF_Pairwise.R dosage_bis,nodosage &
Rscript JASPAR_TF_Pairwise.R MaleOver,FemaleOver &
wait


#############
#############
##  Plots  ##
#############
#############

# RNA-seq coverage tracks: sorted BAM files of fastq2counts.bash (folder dir_mapping in environment.yml)
dir_mapping=$(awk '$1=="dir_mapping:" {print $2}' environment.yml)

samtools merge -f young.bam $(for i in $(seq 1 6); do echo ${dir_mapping}/JD-${i}_sorted.bam; done)
samtools merge -f old.bam $(for i in $(seq 7 24); do echo ${dir_mapping}/JD-${i}_sorted.bam; done)

samtools index old.bam
samtools index young.bam

bamCoverage -b old.bam -o old.bw
bamCoverage -b young.bam -o young.bw
bamCoverage --scaleFactor 0.333 -b old.bam -o old_scaled.bw
