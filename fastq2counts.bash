# !/bin/sh 

###############################
##                           ##
## From Fastq to count table ##
##                           ##
###############################


#################
## Directories ##
#################

dir_fastq=/data/rundata/run174/fastq
dir_home=/data/louisc/Project_Plastoscine/

##########
## STAR ##
##########

cd $dir_fastq
for dir in $(ls)
do 
	cd $dir_fastq/$dir
	pigz -d *
	mkdir $dir_home/$dir
	for files in $(ls)
	do
		# echo $dir_home/$dir/$files
		cd $dir_home/$dir
		STAR --runThreadN 6 --genomeDir $dir_home/TAIR --readFilesIn $dir_fastq/$dir/$files --outFileNamePrefix $files --outReadsUnmapped Fastx
	done
	cd $dir_fastq/$dir
	pigz *
done
 
################
## SAM to BAM ##
################

cd $dir_fastq/
for dir in $(ls)
do 
	cd /$dir_home/$dir/
	for files in $(ls | grep .Aligned.out.sam)
	do 
		samtools view -bS -o $(echo $files | cut -d_ -f 1,3).bam $files 
	done
done

################
## Merge BAMs ##
################

cd $dir_home
for prefix in $(ls -LR | cut -d_ -f 1 | grep JD | sort -u) 
do 
	samtools merge $prefix\_merged $(ls -rt -d -1 $PWD/{*,.*}/* | grep .bam | grep $prefix\_)
done

####################################
## Summarize and make count table ##
####################################

cd $dir_home
for prefix in $(ls -LR | cut -d_ -f 1 | grep JD | sort -u) 
do 
	echo $prefix
	samtools sort $prefix\_merged $prefix\_sorted
	htseq-count -s no -f bam $prefix\_sorted.bam Taeniopygia_guttata.taeGut3.2.4.87.gtf > $prefix\_htseqobj.txt
done

pr -mTs $'\t' <(cut -d$'\t' -f2 JD-1_htseqobj.txt) <(cut -d$'\t' -f2 JD-2_htseqobj.txt) <(cut -d$'\t' -f2 JD-3_htseqobj.txt) <(cut -d$'\t' -f2 JD-4_htseqobj.txt) <(cut -d$'\t' -f2 JD-5_htseqobj.txt) <(cut -d$'\t' -f2 JD-6_htseqobj.txt) <(cut -d$'\t' -f2 JD-7_htseqobj.txt) <(cut -d$'\t' -f2 JD-8_htseqobj.txt) <(cut -d$'\t' -f2 JD-9_htseqobj.txt) <(cut -d$'\t' -f2 JD-10_htseqobj.txt) <(cut -d$'\t' -f2 JD-11_htseqobj.txt) <(cut -d$'\t' -f2 JD-12_htseqobj.txt) <(cut -d$'\t' -f2 JD-13_htseqobj.txt) <(cut -d$'\t' -f2 JD-14_htseqobj.txt) <(cut -d$'\t' -f2 JD-15_htseqobj.txt) <(cut -d$'\t' -f2 JD-16_htseqobj.txt) <(cut -d$'\t' -f2 JD-17_htseqobj.txt) <(cut -d$'\t' -f2 JD-18_htseqobj.txt) <(cut -d$'\t' -f2 JD-19_htseqobj.txt) <(cut -d$'\t' -f2 JD-20_htseqobj.txt) <(cut -d$'\t' -f2 JD-21_htseqobj.txt) <(cut -d$'\t' -f2 JD-22_htseqobj.txt) <(cut -d$'\t' -f2 JD-23_htseqobj.txt) <(cut -d$'\t' -f2 JD-24_htseqobj.txt) > Plastoscine.txt

cut -d$'\t' -f1 JD-1_htseqobj.txt > Plastoscine_annot.txt