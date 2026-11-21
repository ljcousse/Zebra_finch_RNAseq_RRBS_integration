# DNA methylation and gene expression in the developing zebra finch telencephalon

Code for the manuscript *"DNA methylation regulates transcription factor-specific neurodevelopmental but not sexually dimorphic gene expression dynamics in zebra finch telencephalon"* (Diddens, Coussement et al., Frontiers in Cell and Developmental Biology, 2021, [doi:10.3389/fcell.2021.583555](https://doi.org/10.3389/fcell.2021.583555)).

We profiled the telencephalon of male and female zebra finches (*Taeniopygia guttata*) at 1, 20 and 65 days post-hatch (dph) and in adults with RNA-seq and (oxidative) reduced representation bisulfite sequencing ((ox)RRBS). For each sex and age, three pools of three telencephalons were sequenced (24 samples); the same 24 samples were used for RNA-seq and RRBS. Differential expression (edgeR, stagewise testing) and differential methylation (limma on M-values) were tested for the age effect and for the sex effect at each age. The analysis further covers age-dependent expression patterns, transcription factor binding site enrichment (JASPAR 2018, FIMO), the association between methylation and expression (linear mixed models), sexually dimorphic expression and Z chromosome dosage compensation.

## Repository contents

| File | Purpose |
|---|---|
| `FullAnalysis.bash` | Runs the full analysis in order: downloads the genome (Ensembl release 98, taeGut3.2.4) and the JASPAR 2018 motifs, then calls the scripts below, FIMO and the RNA-seq coverage tracks (samtools, deepTools `bamCoverage`) |
| `fastq2counts.bash` | RNA-seq pre-processing: mapping (STAR), SAM → BAM and merging per sample (samtools), counting (htseq-count, Ensembl release 87 GTF) → `Plastoscine.txt` |
| `Plastoscine.txt`, `Plastoscine_annot.txt` | RNA-seq count table (24 samples, `JD-1` to `JD-24`) and the Ensembl gene ID per row |
| `DE_analysis.R` | Differential expression: edgeR GLM with omnibus test and stagewise post-hoc testing (stageR); gene lists per contrast, boxplots, `FullData_output.txt` |
| `DE_analysis_statisticsPaper.R` | Summary statistics and sample clustering of the RNA-seq data |
| `Circos.R` | Circos plot (RCircos) of the age- and sex-related genes over the genome (uses `genome.txt`) |
| `RNAseq_ClusterAnalysis.R` | Classification of genes into 12 age-dependent expression patterns |
| `ExtractSeq_TFsearching_RNAseq.R` | Promoter sequences (−2,000 to +500 bp around the TSS) of all expressed genes for FIMO |
| `Parse_FimoAll_RNAseq.R`, `Parse_FimoAll_RRBS.R` | Merge the FIMO output of all motifs into one table |
| `JASPAR_TF_Pairwise.R` | TF binding site enrichment: chi-squared test per motif, genes of interest vs an outgroup |
| `RRBS_analyse.R` | Differential CpG methylation: import of the Bismark files (BiSeq), filtering, M-values, annotation (promoter/genic/intergenic), limma per contrast, clustering |
| `RRBS_analyse_hydroxy.R` | Hydroxymethylation (RRBS − oxRRBS) |
| `RRBS_analyse_nonCpG.R` | Differential CpH (CHG and CHH) methylation, same approach as `RRBS_analyse.R` |
| `RRBS_plots.R` | Combined figure of CpG, hydroxymethylation and CpH methylation |
| `RRBS_statisticsPaper.R` | Summary statistics, correlation and linear mixed models between methylation and expression, differential methylation in DE vs non-DE genes, bedGraph tracks |
| `ExtractSeq_TFsearching_RRBS.R`, `JASPAR_TF_Pairwise_RRBS.R` | TF binding site enrichment around differentially methylated promoter CpGs (±50 bp) |
| `Expression_methylation_cor.R` | Correlation between promoter methylation and expression for the age-related CpGs (CpG and CpH) |
| `Sexually_dimorphic_genes.R` | Distribution of sex-related genes over autosomes and the Z chromosome; male- and female-overexpressed Z genes |
| `DosageEffect.R` | Dosage compensation: M:F expression ratios, dosage-compensated (from birth or acquired) and non-compensated Z genes, expression and methylation per group |
| `Info_paper.R` | Mapping and trimming statistics for the supplementary tables (run on the mapping servers) |
| `AnnotFromCSV.R` | Builds `JASPAR_annot.txt` (JASPAR 2018 motif annotation) from the JASPAR search exports |
| `JASPAR_annot.txt` | JASPAR 2018 motif annotation (ID, name, species, class, family) |
| `genome.txt` | Chromosome sizes of taeGut3.2.4 (for `Circos.R`) |
| `conda_environment.yml`, `conda_environment_preprocessing.yml` | Conda environments for the analysis and for the RNA-seq pre-processing |
| `environment_example.yml` | Template for the local paths (copy to `environment.yml`) |

## Installation

```bash
git clone https://github.com/ljcousse/PlastoscineProject.git
cd PlastoscineProject

conda env create -f conda_environment.yml
conda activate zebra_finch
```

For the RNA-seq pre-processing (`fastq2counts.bash`), create the second environment: `conda env create -f conda_environment_preprocessing.yml`.

The original analysis was run on Linux servers without conda, with R 3.6.1, edgeR 3.8.6, limma 3.22.7, UniProt.ws 2.14.0, STAR 2.5.2b, samtools 0.1.18, htseq-count 0.6.0, FastQC 0.10.1, Trim Galore, Bismark 0.9.0 (Bowtie2 mode), FIMO 4.11.3 and deepTools 3.4.3 (versions as in the manuscript). These versions come from different Bioconductor releases and cannot be installed together with conda. The conda environments therefore use R 3.6 with Bioconductor 3.10, MEME 4.11.2, deepTools 3.4.3, STAR 2.5.2b, samtools 0.1.18 and htseq-count 0.6.1. All pinned versions have Linux builds on conda-forge/bioconda; results may differ slightly from the published ones. If an environment does not solve, use `conda config --set channel_priority strict` or `mamba`.

The annotation queries the Ensembl January 2019 archive (release 95) with biomaRt and UniProt with UniProt.ws. The UniProt web service has changed since the analysis, so `UniProt.ws` will likely fail; the cached `uniprot.Rda`, `Genes.Rda` and `Genes_RRBS.Rda` files are then needed in the repository folder.

The RRBS pre-processing (Trim Galore, Bismark) is not part of this repository; the analysis starts from the Bismark coverage files.

## Configuration

Local paths are not hard-coded. Copy the template and fill in your paths:

```bash
cp environment_example.yml environment.yml
```

`environment.yml` contains:
- `dir_rrbs_cpg`: Bismark CpG coverage files, RRBS (`*.B*`) and oxRRBS (`*.O*`) per sample;
- `dir_rrbs_chg`, `dir_rrbs_chh`: Bismark CHG and CHH coverage files;
- `dir_fastq`, `dir_mapping`, `dir_star_index`, `file_gtf`: FASTQ files, mapping output folder, STAR index and Ensembl release 87 GTF (for `fastq2counts.bash`);
- `dir_rrbs_runs`: RRBS run folders with the trimming reports (for `Info_paper.R`).

Run the scripts from the repository folder; all output is written there (the output folders and files are listed in `.gitignore`). Keep `environment.yml` local (it is listed in `.gitignore`).

## Running the analysis

0. **RNA-seq pre-processing (optional):** `bash fastq2counts.bash` (FASTQ → STAR alignments → merged and sorted BAM per sample → htseq-count → `Plastoscine.txt`). The count table is already in the repository.
1. **Full analysis:** `bash FullAnalysis.bash`. It runs, in this order:
   1. downloads: genome FASTA → `Genomic_data/`, JASPAR 2018 motifs → `JASPAR_pfms_meme/`;
   2. RNA-seq: `DE_analysis.R` → `GeneLists/`, `FullData_output.txt`; `DE_analysis_statisticsPaper.R`; `Circos.R`;
   3. TF analysis (RNA-seq): `ExtractSeq_TFsearching_RNAseq.R` → FIMO per motif (40 jobs in parallel) → `Parse_FimoAll_RNAseq.R` → `JASPAR_TF_Pairwise.R` per contrast → `JASPAR_output/`;
   4. expression patterns: `RNAseq_ClusterAnalysis.R` → `GeneLists/ListGenes_patt_*.txt`, `Trends/`, then `JASPAR_TF_Pairwise.R` per pattern;
   5. RRBS: `RRBS_analyse.R`, `RRBS_analyse_hydroxy.R`, `RRBS_analyse_nonCpG.R` → `RRBS_<contrast>_full.txt`, `RRBS_<contrast>_RNAseq.txt`; `RRBS_statisticsPaper.R`; `RRBS_plots.R`;
   6. TF analysis (RRBS): `ExtractSeq_TFsearching_RRBS.R` → FIMO → `Parse_FimoAll_RRBS.R` → `JASPAR_TF_Pairwise_RRBS.R age`;
   7. `Expression_methylation_cor.R`, `Sexually_dimorphic_genes.R`, `DosageEffect.R` and the TF analysis of the dosage and sex gene lists;
   8. RNA-seq coverage tracks: sorted BAM files of `fastq2counts.bash` (from `dir_mapping`) merged per age group (1 dph vs older) → `bamCoverage` → `young.bw`, `old.bw`, `old_scaled.bw`.

   FIMO and `Parse_Fimo*.R` use 40 cores; change `N=40` in `FullAnalysis.bash` and `mc.cores = 40` for smaller machines. Steps that run in parallel are awaited before the next step starts. The gene set analysis (Genomatix GePS, Metascape) was done outside this repository on the gene lists in `GeneLists/`.
2. **Sequencing statistics:** `Rscript Info_paper.R porthos` (RNA-seq mapping) or `Rscript Info_paper.R midas` (RRBS trimming), run on the servers with the mapping output.

### Main outputs

| Output | Content |
|---|---|
| `FullData_output.txt` | Differential expression results for all genes and contrasts |
| `GeneLists/ListGenes_*.txt` | Significant genes per contrast (`1dph`, `20dph`, `65dph`, `adult`, `age`), per expression pattern, dosage group and sex |
| `RRBS_<contrast>_full.txt`, `RRBS_<contrast>_RNAseq.txt` | Differential methylation per CpG for all CpGs and for the CpGs in DE genes (FDR recalculated); `_nonCpG` for CpH |
| `JASPAR_output/JASPAR_TF_*_comp_h.txt` | TF binding site enrichment per gene list |
| `dat_lme_1_prom_CpG.txt`, `dat_lme_1_genic_CpG.txt` (and `_nonCpG`) | Input of the methylation–expression mixed models |
| `*.tiff` | Figures (clustering, Circos, expression patterns, dosage compensation, methylation) |
| `*.BedGraph`, `young.bw`, `old.bw` | UCSC genome browser tracks (taeGut2) |

### Statistical design (short)

- **Differential expression:** edgeR GLM `~ age * sex` after TMM normalisation; genes with cpm > 1 in at least 3 samples are kept. Contrasts: sex effect at 1, 20 and 65 dph and adult, and the age effect (mean of 20 dph, 65 dph and adult vs 1 dph). Omnibus test over the five contrasts, then stagewise post-hoc testing (stageR, Holm, α = 0.05).
- **Differential methylation:** limma on M-values (`log2((β + 0.01)/(1 − β + 0.01))`), same design and contrasts. CpGs with at least 5 reads in all samples and at least 6 methylated reads in total are kept. Benjamini-Hochberg FDR; for the CpGs in DE genes the FDR is recalculated on that subset.
- **Expression patterns:** mean log-cpm at 20 dph, 65 dph and adult compared with 1 dph; higher/lower when the difference exceeds 1.65 SD, equal within 1 SD; 12 patterns.
- **TF binding sites:** FIMO matches (p < 1e-4) of the JASPAR 2018 non-redundant vertebrate motifs in promoters (−2,000 to +500 bp) or around DM CpGs (±50 bp); chi-squared test per motif against an outgroup, Benjamini-Hochberg; reported when FDR < 0.10 and |log2 FC| > 1.
- **Methylation vs expression:** linear mixed model (lmerTest) `log-cpm ~ M-value + (1 | sample) + (1 + M-value | gene)` for promoter and genic CpGs separately.
- **Sexually dimorphic genes and dosage compensation:** chi-squared tests (Bonferroni) for autosomal vs Z genes; dosage compensation from the 95% quantile of the autosomal M:F ratios; ANOVA with Tukey post hoc.

Promoters are defined as −2,000 to +500 bp around the transcription start site; genic CpGs are the remaining CpGs within a gene.

## Data

The RNA-seq and (ox)RRBS data are available in NCBI GEO under accession [GSE147974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147974).

The RNA-seq count table (`Plastoscine.txt`, `Plastoscine_annot.txt`) is part of this repository. Samples `JD-1` to `JD-24` are ordered by age (1 dph, 20 dph, 65 dph, adult) and, within age, three female then three male pools. The sample annotation of the RRBS files (pool number → sex and age) is written in the RRBS scripts. The labels of the bisulfite (`B16`) and oxidative bisulfite (`O16`) files of pool 16 are switched; the RRBS scripts correct this when reading the files.

The analysis also expects in the repository folder:
- the Bismark coverage files (via `environment.yml`);
- the cached annotation files `uniprot.Rda`, `Genes.Rda`, `Genes_RRBS.Rda` if the UniProt and Ensembl archive services are not available.

## License

This code is released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite:

Diddens J\*, Coussement L\*, Frankl-Vilches C, Majumdar G, Steyaert S, Ter Haar SM, Galle J, De Meester E, De Keulenaer S, Van Criekinge W, Cornil CA, Balthazart J, Van Der Linden A, De Meyer T, Vanden Berghe W. DNA methylation regulates transcription factor-specific neurodevelopmental but not sexually dimorphic gene expression dynamics in zebra finch telencephalon. Frontiers in Cell and Developmental Biology. 2021;9:583555. https://doi.org/10.3389/fcell.2021.583555.

\* shared first authorship
