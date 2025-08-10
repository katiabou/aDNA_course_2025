# aDNA_course_2025
Genotype imputation exercise for the 2025 aDNA course "Decoding the past: An introduction to ancient genomics"

### Introduction

This tutorial will show you how to run genotype imputation using GLIMPSE ([Rubinacci et al.,2021](https://doi.org/10.1038/s41588-020-00756-0)), an imputation and phasing method tailored for low coverage sequencing data. We will test how imputation works in an ancient dog from Siberia, dated to 2,000 years before present (ybp). This sample is relativaly high coverage for ancient standards (11.8x), so in order to test the imputation accuracy, we will downsample it to lower coverages (0.1x and 1x) and see how the imputed genotypes compare to the "true" high coverage ones.

We will run GLIMPSE from sequencing reads data (BAM files) to obtain imputed and phased genotypes.

For this tutorial, we will focus only on chromosome 22.

### Set paths to data folder
```
COURSE_PATH=/projects/course_1
DATA_PATH=${COURSE_PATH}/people/qcj125/data_exercise
```

### Load required modules
Run the following bash script:  
`. ${COURSE_PATH}/people/qcj125/activate.sh`


## Exercise 1: Explore input data
We will use the following files:  
1) Reference panel  
2) Genetic (recombination) map  
3) Reference fasta file  
4) Ancient dog bam files  

### Reference panel
`REF=${DATA_PATH}/reference_panel_vcf/ref-panel_chr22.vcf.gz`

Let's take a look:  
`bcftools view ${REF} | less -S`

#### Questions: 
1) How many sites are in the panel for chr22? (813342)
2) Are the sites phased?
3) How many samples are in the reference panel? (1519)
4) What species are present?

### Bam files
These include the initial high coverage bam file and the downsampled ones to 1x and 0.1x coverage:
```
${DATA_PATH}/bams/TRF.05.05_Merged.CanFam3.1Y.rmdup.bam
${DATA_PATH}/bams_down/TRF.05.05.chr22.0.1x.bam
${DATA_PATH}/bams_down/TRF.05.05.chr22.1x.bam
```

### Genetic map (recombination map)
`MAP=${DATA_PATH}/gen_map/chr22_average_canFam3.1_modified.txt`

Let's take a look:  
`less ${MAP}`

The three columns are genomic position, chromosome, and distance in cM.

#### Question: 
Why do you think we need a genetic map for imputation?

### Reference genome (fasta files)
`REFGEN=${DATA_PATH}/reference_fasta/CanFam31_chr22.fasta`

**Note!!!** It is important to make sure the reference panel samples and target samples have been mapped with the same reference genome, and that you use the correct one for imputation.


## Exercise 2: Imputing an ancient dog genome with GLIMPSE
We will now run different steps as part of the GLIMPSE imputation pipeline and impute the 1x ancient Siberian dog sample.

#Make output directories (ONLY RUN THESE ONCE!)
```
mkdir -p output/reference_panel
mkdir -p output/GLs_target_bams
mkdir -p output/chunks
mkdir -p output/GLIMPSE_imputed
mkdir -p output/GLIMPSE_ligated
```

### Step 1: Extract variable positions from the reference panel  
The first step is to extract the list of variable sites from the reference panel in a tsv format, to be used downstream for computing genotype likelihoods in the ancient samples (see below)

The -G option drops individual genotype information (not needed at this step)
The -m 2 -M 2 -v snps retains biallelic snps

```
REF=${DATA_PATH}/reference_panel_vcf/ref-panel_chr22.vcf.gz
VCF=output/reference_panel/chr22_ref_panel_sites.vcf.gz

bcftools view \
-G -m 2 -M 2 -v snps \
${REF} \
-Oz -o ${VCF}

bcftools index -f ${VCF}
```

Get the chrom, pos, ref and alt fields for each site and put it into a new tsv file:
```
VCF=output/reference_panel/chr22_ref_panel_sites.vcf.gz
TSV=output/reference_panel/chr22_ref_panel_sites.tsv.gz

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' \
${VCF} | \
bgzip -c > ${TSV}

tabix -s1 -b2 -e2 ${TSV}
```

Take a look at the tsv file:    
`less -S output/reference_panel/chr22_ref_panel_sites.tsv.gz`



### Step 2: Compute genotype likelihoods at each position

# Let's start with the 1x ancient Siberian dog
SAMPLE=TRF.05.05.chr22.0.1x

The bcftools mpileup command generates a VCF containing genotype likelihoods for a bam file at specified sites.
The output of this step is a VCF file format containing genotype likelihoods at each site (based on those presentin the reference panel).

```
BAM={DATA_PATH}/bams_down/${SAMPLE}.bam
VCF=output/reference_panel/chr22_ref_panel_sites.vcf.gz
TSV=output/reference_panel/chr22_ref_panel_sites.tsv.gz
REFGEN={DATA_PATH}/reference_fasta/CanFam31_chr22.fasta
GL=output/GLs_target_bams/${SAMPLE}.vcf.gz

bcftools mpileup \
-f ${REFGEN} \
-I -E -a 'FORMAT/DP' \
-T ${VCF} \
-r chr22 ${BAM} -Ou | \
bcftools call \
-Aim -C alleles -T ${TSV} \
-Oz -o ${GL} \
--threads 6

bcftools index -f ${GL}
```

### Step 3: Splitting each chromosome into chunks  
Glimpse has a tool (GLIMPSE_chunk) which splits the chromosome into "chunks" to run the imputation and phasing on each.

```
VCF=output/reference_panel/chr22_ref_panel_sites.vcf.gz
CHUNK=output/chunks/chr22_chunks.txt

GLIMPSE_chunk_static \
--input ${VCF} \
--region chr22 \
--window-size 2000000 --buffer-size 200000 \
--output ${CHUNK}
```

### Step 4: Impute! For each of the chunks estimated from above, impute using the genotype likelihoods 

```
GL=output/GLs_target_bams/${SAMPLE}.vcf.gz
REF_PANEL={DATA_PATH}/reference_panel_vcf/ref-panel_chr22.vcf.gz
MAP={DATA_PATH}/gen_map/chr22_average_canFam3.1_modified.txt
CHUNK=output/chunks/chr22_chunks.txt
OUT=output/GLIMPSE_imputed/${SAMPLE}.${ID}.bcf
PREFIX=output/GLIMPSE_imputed/${SAMPLE}

while IFS="" read -r LINE || [ -n "${LINE}" ];
do
  printf -v ID "%02d" $(echo ${LINE} | cut -d" " -f1)
  IRG=$(echo ${LINE} | cut -d" " -f3)
  ORG=$(echo ${LINE} | cut -d" " -f4)
  OUT=${PREFIX}.${ID}.bcf
  GLIMPSE_phase_static --input ${GL} \
  --reference ${REF_PANEL} \
  --map ${MAP} \
  --input-region ${IRG} \
  --output-region ${ORG} --output ${OUT} \
  --thread 2
  bcftools index -f ${OUT}
done < ${CHUNK}
```

Take a small break while this is running :)

### Step 5: Ligate the chunks together
After imputing, we can now merge all of the imputed sites for each chunk together, using the GLIMPSE_ligate tool:

```
#make list with file names to merge:
LIST=output/GLIMPSE_ligated/list.chr22.txt
ls output/GLIMPSE_imputed/${SAMPLE}.*.bcf > ${LIST}

LIGATED=output/GLIMPSE_ligated/${SAMPLE}.ligated.bcf

GLIMPSE_ligate_static \
--input ${LIST} \
--output ${LIGATED}

bcftools index -f ${LIGATED}
```
### Step 6: Final step is to phase the whole chromosome using the GLIMPSE_sample tool. 
We will use the mode which outputs the most likely haplotypes given the values in the FORMAT/HS field. (What is the HS field?)

```
LIGATED=output/GLIMPSE_ligated/${SAMPLE}.ligated.bcf
PHASED=output/GLIMPSE_ligated/${SAMPLE}.phased.bcf

GLIMPSE_sample_static \
--input ${LIGATED} \
--solve \
--output ${PHASED}

bcftools index -f ${PHASED}
```

#### Congrats, you have now imputed and phased an ancient dog sample!

