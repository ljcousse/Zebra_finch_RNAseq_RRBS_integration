#!/bin/sh

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
wget http://jaspar2018.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.zip --user-agent "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.22 (KHTML, like Gecko) Ubuntu Chromium/25.0.1364.160 Chrome/25.0.1364.160 Safari/537.22"
unzip *
rm *.zip
cd ..

mkdir JASPAR_counts
mkdir JASPAR_output

# Install Necessary R packages
Rscript InstallPackages.R



##############
##############
##  RNAseq  ##
##############
##############


######################
# Age and Sex effect #
######################

# Download RNAseq data from GEO
# mkdir RNAseq_data
# cd RNAseq_data
# wget GEO-link
# cd ..

# Preprocessing RNAseq
# mkdir RNASeq_preproc
# bash fastq2counts.bash # Adjust folders!

# OR

# Download count and annotation data from GEO
# implement after GEO submission is filed


# Differential expression analysis
mkdir GeneLists
mkdir boxplots
Rscript DE_analysis.R

# General statistics & Plots for paper
Rscript DE_analysis_statisticsPaper.R
Rscript Circos.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis 
mkdir fimo_RNAseq_out
Rscript ExtractSeq_TFsearching_RNAseq_new.R
N=40 
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RNAseq_out/fimo_$file JASPAR_pfms_meme/$file seq_RNAseq_TFsearching.fa &
done
)

Rscript Parse_FimoAll_RNAseq.R

for coi in "1dph" "20dph" "65dph" "adult" "age"
do
	Rscript JASPAR_TF_Pairwise.R $coi,nonsig &
done


####################
# Cluster analysis #
####################g

# Make clustering 
mkdir Trends
Rscript RNAseq_ClusterAnalysis.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis
for coi in $(ls GeneLists | grep patt | cut -d. -f 1 | cut -d_ -f 2,3,4,5)
do
	Rscript JASPAR_TF_Pairwise.R $coi &
done



############
############
##  RRBS  ##
############
############

# Preprocessing RRBS done by Sandra Steyaert

# Download RRBS data from GEO
# mkdir RRBSdata
# cd RRBS_data
# wget GEO-link
# cd ..
# Adjust folders in RRBS_analyse* scripts!

# Differential methylation analysis
Rscript RRBS_analyse.R 
Rscript RRBS_analyse_hydroxy.R
Rscript RRBS_analyse_nonCpG.R

# General statistics & Plots for paper
Rscript RRBS_statisticsPaper.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis 
mkdir fimo_RRBS_out
Rscript ExtractSeq_TFsearching_RRBS_new.R
N=40 
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RRBS_out/fimo_$file JASPAR_pfms_meme/$file seq_RRBS_TFsearching.fa &
done
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

Rscript DosageEffect

# TF analysis
Rscript JASPAR_TF_Pairwise.R dosage,nodosage &
Rscript JASPAR_TF_Pairwise.R dosage_bis,nodosage &
Rscript JASPAR_TF_Pairwise.R MaleOver,FemaleOver &


#############
#############
##  Plots  ##
#############
#############

samtools merge -f young.bam JD-1_sorted.bam JD-2_sorted.bam JD-3_sorted.bam JD-4_sorted.bam JD-5_sorted.bam JD-6_sorted.bam 
samtools merge -f old.bam JD-7_sorted.bam JD-8_sorted.bam JD-9_sorted.bam JD-10_sorted.bam JD-11_sorted.bam JD-12_sorted.bam JD-13_sorted.bam JD-14_sorted.bam JD-15_sorted.bam JD-16_sorted.bam JD-17_sorted.bam JD-18_sorted.bam JD-19_sorted.bam JD-20_sorted.bam JD-21_sorted.bam JD-22_sorted.bam JD-23_sorted.bam JD-24_sorted.bam 

samtools index old.bam
samtools index young.bam

bamCoverage -b old.bam -o old.bw
bamCoverage -b young.bam -o young.bw
bamCoverage --scaleFactor 0.333 -b old.bam -o old_scaled.bw