#!/bin/sh

###########################################
###########################################
##  Make data structure & download data  ##
###########################################
###########################################


# Cluster analysis
mkdir Trends

# TF analysis
mkdir fimo_RNAseq_out

mkdir JASPAR_pfms_meme
cd JASPAR_pfms_meme
wget http://jaspar2018.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.zip --user-agent "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.22 (KHTML, like Gecko) Ubuntu Chromium/25.0.1364.160 Chrome/25.0.1364.160 Safari/537.22"
cd ..
unzip *
rm *.zip

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

# Preprocessing RNAseq
# bash fastq2counts.bash

# Differential expression analysis
Rscript DE_analysis.R

# General statistics & Plots for paper
Rscript DE_analysis_statisticsPaper.R
Rscript Circos.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis 
Rscript ExtractSeq_TFsearching_RNAseq_new.R
N=40 
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RNAseq_out/fimo_$file JASPAR_pfms_meme/$file seq_RNAseq_TFsearching.fa &
done
)

for coi in "1dph" "20dph" "65dph" "adult" "age"
do
	Rscript JASPAR_TF_Pairwise.R $coi,nonsig &
done
Rscript JASPAR_TF_Pairwise.R dosage,nodosage &
Rscript JASPAR_TF_Pairwise.R dosage_bis,nodosage &
Rscript JASPAR_TF_Pairwise.R MaleOver,FemaleOver &


####################
# Cluster analysis #
####################g

# Make clustering 
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

# Differential methylation analysis
Rscript RRBS_analyse.R
Rscript RRBS_analyse_hydroxy.R
Rscript RRBS_analyse_nonCpG.R

# General statistics & Plots for paper
Rscript RRBS_statisticsPaper.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis 
Rscript ExtractSeq_TFsearching_RRBS_new.R
N=40 
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RRBS_out/fimo_$file JASPAR_pfms_meme/$file seq_RRBS_TFsearching.fa &
done
)



############################################
############################################
##  Correlation Methylation & Expression  ##
############################################
############################################


Rscript Expression_methylation_cor



#####################
#####################
##  Dosage Effect  ##
#####################
#####################


Rscript DosageEffect



#############
#############
##  Plots  ##
#############
#############


