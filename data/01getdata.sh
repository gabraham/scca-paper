#!/bin/bash

set -e

# Download the genotype data:
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2

# Download the gene expression data from
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-264/E-MTAB-264.processed.1.zip
unzip E-MTAB-264.processed.1.zip

bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2


grep Normalization *_p3_expression.txt | cut -f 2- -d$'\t' | tr '\t' '\n' \
   > expression_samples.txt

wc -l expression_samples.txt

awk 'NR == FNR {a[$1]; next} $2 in a {print $1, $2}' \
   expression_samples.txt hapmap3_r2_b36_fwd.consensus.qc.poly.ped \
   > common_samples.txt

plink \
   --file hapmap3_r2_b36_fwd.consensus.qc.poly \
   --make-bed \
   --maf 0.05 \
   --hwe 5e-6 \
   --mind 0.05 \
   --geno 0.05 \
   --filter-founders \
   --keep common_samples.txt \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered

plink \
   --bfile hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered \
   --chr 1 \
   --make-bed \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_chr1


