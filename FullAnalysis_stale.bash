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
wget http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.zip
unzip *
rm *.zip

mkdir JASPAR_counts
mkdir JASPAR_output



##############
##############
##  RNAseq  ##
##############
##############


# Preprocessing RNAseq
bash fastq2counts.bash

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



########################
########################
##  Cluster analysis  ##
########################
########################

# Make clustering 
Rscript RNAseq_ClusterAnalysis.R

# Gene set analysis
## INPUT JOLIEN/WIM

# TF analysis
for coi in $(ls GeneLists | grep patt | cut -d. -f 1 | cut -d_ -f 2,3,4,5)
do
	Rscript JASPAR_TF_Pairwise.R $coi &
done



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








jaspar2meme -pfm TF > Fimo_in.txt

# RNAseq
########

Rscript ExtractSeq_TFsearching_RNAseq.R 

for coi in "1dph" "20dph" "65dph" "adult" "age"
do
	fimo -oc fimo_$coi\_out Fimo_in.txt seq_$coi\_TFsearching.fa &
done
fimo -oc fimo_RNAseq_ns_out Fimo_in.txt seq_RNAseq_ns_TFsearching.fa

for coi in "1dph" "20dph" "65dph" "adult" "age"
do
	Rscript JASPAR_TF_Pairwise.R $coi,nonsig &
done
Rscript JASPAR_TF_Pairwise.R dosage,nodosage &
Rscript JASPAR_TF_Pairwise.R dosage_bis,nodosage &
Rscript JASPAR_TF_Pairwise.R MaleOver,FemaleOver &
for coi in $(ls GeneLists | grep patt | cut -d. -f 1 | cut -d_ -f 2,3,4,5)
do
	Rscript JASPAR_TF_Pairwise.R $coi &
done


fimo -oc fimo_RNAseq_ns_out Fimo_in.txt seq_RNAseq_ns_TFsearching.fa

Rscript JASPAR_TF_RNAseq.R

# RRBS
######

Rscript ExtractSeq_TFsearching_RNAseq.R 

for coi in "age"
do
	fimo -oc fimo_$coi\_rrbs_out Fimo_in.txt seq_$coi\_rrbs_TFsearching.fa &
done
fimo -oc fimo_rrbs_ns_out Fimo_in.txt seq_rrbs_ns_TFsearching.fa

Rscript JASPAR_TF_RRBS.R

# Dosage Compensation
#####################

Rscript ExtractSeq_TFsearching_DosageComp.R 

for coi in "dosage" "nodosage" "dosage_bis"
do
	fimo -oc fimo_$coi\_out Fimo_in.txt seq_$coi\_TFsearching.fa &
done

Rscript JASPAR_TF_DosageComp.R

Rscript JASPAR_TF_Pairwise.R dosage,nodosage &
Rscript JASPAR_TF_Pairwise.R dosage_bis,nodosage &
Rscript JASPAR_TF_Pairwise.R MaleOver,FemaleOver &

# Trends
#####################

Rscript ExtractSeq_TFsearching_Trends.R 

cd Trends/

for coi in $(ls | grep genelist | cut -d. -f 1 | awk -F 'genelist_patt' '{print $2}')
do
	fimo -oc fimo_$coi\_out /home/louisc/Zebra_Finch_new/Fimo_in.txt seq_$coi\_TFsearching.fa &
	fimo -oc fimo_$coi\_ns_out /home/louisc/Zebra_Finch_new/Fimo_in.txt seq_$coi\_ns_TFsearching.fa &
done

cd ..

Rscript JASPAR_TF_DosageComp.R

# Dosage Compensation
#####################

Rscript ExtractSeq_TFsearching_M_F.R

for coi in "MaleOver" "FemaleOver"
do
	fimo -oc fimo_$coi\_out Fimo_in.txt seq_$coi\_TFsearching.fa &
done

Rscript JASPAR_TF_DosageComp.R



for coi in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
do
	fimo -oc seq_RNAseq_TFsearching/fimo_$coi\_out --thresh 1 --max-stored-scores 10000000 Fimo_in.txt seq_RNAseq_TFsearching/seq_$coi\.fa &
done

fimo -oc seq_test_out --thresh 1 --max-stored-scores 10000000 TF/MA0074.1.pfm seq_test.fa &
done

fimo -oc fimo_RNAseq_out Fimo_in_Vert.txt seq_RNAseq_TFsearching.fa


## RNAseq
N=40 
(
for file in $(ls JASPAR_pfms_meme)
do
	((i=i%N)); ((i++==0)) && wait
	fimo -oc fimo_RNAseq_out/fimo_$file JASPAR_pfms_meme/$file seq_RNAseq_TFsearching.fa &
done
)
## RRBS
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



########
# MEME #
########

meme -maxsize 10000000 -dna -w 10 -revcomp -p 20 -nmotifs 10 -o dosage_meme_out seq_dosage_TFsearching.fa
meme -maxsize 10000000 -dna -w 10 -revcomp -p 20 -nmotifs 10 -o nodosage_meme_out seq_nodosage_TFsearching.fa


#######
# End #
#######



Rscript Trend_expression.R

Rscript RRBS_analyse.R

# Correlation Exrpession & Methylation
Rscript Expression_methylation_cor.R

# TF analysis
Rscript ExtractSeq_TFsearching.R

fimo -oc fimo_out Fimo_in.txt seq_TFsearching.fa 
fimo -oc fimo_rna_age_out Fimo_in.txt seq_rna_age_TFsearching.fa
fimo -oc fimo_rna_int_out Fimo_in.txt seq_rna_int_TFsearching.fa
fimo -oc fimo_rna_sex_out Fimo_in.txt seq_rna_sex_TFsearching.fa

fimo -oc fimo_ns_out Fimo_in.txt seq_ns_TFsearching.fa 
fimo -oc fimo_rna_age_ns_out Fimo_in.txt seq_rna_age_ns_TFsearching.fa
fimo -oc fimo_rna_int_ns_out Fimo_in.txt seq_rna_int_ns_TFsearching.fa
fimo -oc fimo_rna_sex_ns_out Fimo_in.txt seq_rna_sex_ns_TFsearching.fa

Rscript JASPAR_TF.R
Rscript JASPAR_TF_RNAseq.R

# Cluster analysis

Rscript GetClustering.R

Rscript ExtractSeq_TFsearching_clusters.R

for cluster in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26; do
  fimo -oc clustering/cluster_$cluster Fimo_in.txt clustering/seq_rna_int_TFsearching_cluster$cluster.fa &
done

Rscript JASPAR_TF_clusters.R
Rscript Compare_TF.R

# Trend analysis

Rscript TrendsInExpression.R
Rscript JASPAR_TF_trends.R

Rscript Compare_TF_trends.R

# Sex analysis

Rscript DivideSexes.R

fimo -oc fimo_rna_sex_Mbias_out Fimo_in.txt seq_rna_sex_Mbias_TFsearching.fa &
fimo -oc fimo_rna_sex_Fbias_out Fimo_in.txt seq_rna_sex_Fbias_TFsearching.fa &

Rscript JASPAR_TF_DividedSexes.R