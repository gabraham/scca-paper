#!/bin/bash

set -e

# See ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/00README.txt
# for information about the genotyping QC procedures

# Download the genotype data:
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/relationships_w_pops_121708.txt

# Download the gene expression data
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-264/E-MTAB-264.processed.1.zip
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-198/E-MTAB-198.processed.1.zip

Rscript splitpops.R

unzip E-MTAB-264.processed.1.zip
unzip E-MTAB-198.processed.1.zip

# Rename CEU for consistency (ignore CEU RNA-seq data)
cp normalized_array_data.109.txt CEU_p3_expression.txt

bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2

# Convert genotypes to binary format
plink \
   --file hapmap3_r2_b36_fwd.consensus.qc.poly \
   --make-bed \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly

# Get the sample IDs for the gene expression data
grep Normalization *_p3_expression.txt \
   | cut -f 2- -d$'\t' | tr '\t' '\n' \
   > expression_samples.txt

wc -l expression_samples.txt

awk 'NR == FNR {a[$1]; next} $2 in a {print $1, $2}' \
   expression_samples.txt hapmap3_r2_b36_fwd.consensus.qc.poly.ped \
   > common_samples.txt

MERGE=merge_list.txt
rm -f $MERGE

# Split out genotypes by population
for f in *_p3_expression.txt
do
   POP=${f/_p3_expression.txt/}
   
   COM=common_samples_${f/_ap3_expression/}
   # Keep only samples with gene expression data
   awk 'NR == FNR {a[$2]=$1; next} $2 in a {print a[$2], $2}' \
      common_samples.txt \
      genotyped_samples_${POP}.txt > $COM

   wc -l $COM

   # Filter SNPs based on all individuals, not just
   # the ones with gene expression data
   plink \
      --bfile hapmap3_r2_b36_fwd.consensus.qc.poly \
      --make-bed \
      --maf 0.05 \
      --hwe 1e-6 \
      --mind 0.05 \
      --geno 0.1 \
      --autosome \
      --filter-founders \
      --keep genotyped_samples_${POP}.txt \
      --out hapmap3_r2_b36_fwd.consensus.qc.poly_${POP}_filtered

   # Keep only the individuals with the gene expression data
   plink \
      --bfile hapmap3_r2_b36_fwd.consensus.qc.poly_${POP}_filtered \
      --make-bed \
      --keep $COM \
      --autosome \
      --out hapmap3_r2_b36_fwd.consensus.qc.poly_${POP}_filtered_common

   echo hapmap3_r2_b36_fwd.consensus.qc.poly_${POP}_filtered_common.{bed,bim,fam} \
      >> $MERGE
done

plink \
   --merge-list merge_list.txt \
   --make-bed \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_filtered_common

plink \
   --bfile hapmap3_r2_b36_fwd.consensus.qc.poly_filtered_common \
   --chr 1 \
   --make-bed \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_filtered_common_chr1


