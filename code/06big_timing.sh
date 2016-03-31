#!/bin/bash

set -e

ROOT=../data/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.qc.poly.filtered.wexpr

PLINK=plink
FPCA=flashpca

type $PLINK >/dev/null 2>&1 || {
   echo >&2 "plink is not available; set PLINK variable in 06big_timing.sh script or fix path"
   exit 1
}
type $FPCA >/dev/null 2>&1 || { 
   echo >&2 "flashpca is not available; set FPCA variable in 06big_timing.sh script or fix path"
   exit 1
}


EXPR1=expression_standardised.txt
EXPR2=expression_standardised_10000.txt

nreps=30

# Subset of genes
awk '{for(i=1;i<10002;i++){printf "%s ", $i} print $10002}' $EXPR1 > $EXPR2

CHR=$(seq 2 22)

for expr in $EXPR1 $EXPR2
do
   ################################################################################
   # Analyses of chr 1-22
   for chr in $CHR
   do
      OUT=$(basename $ROOT)_chr1to${chr}
   
      eval $PLINK \
         --bfile $ROOT \
         --chr 1-${chr} \
         --make-bed \
         --out ${OUT}
   done
   
   for rep in $(seq $nreps)
   do
      for chr in $CHR
      do
         OUT=$(basename $ROOT)_chr1to${chr}
   
         (time eval $FPCA \
            --scca \
            --bfile ${OUT} \
            --stand sd \
            --mem low \
            --pheno $expr \
            --lambda1 1e-3 \
            --lambda2 1e-3 \
            --ndim 1 \
            --v) 2>&1 | tee scca_${expr}_chr1to${chr}_rep${rep}.log
      done
   done
   
   ################################################################################
   # Analyses of subsets of chr 1
   for f in snps_*.txt
   do
      OUT=$(basename $ROOT)_chr1
      eval $PLINK \
         --bfile $ROOT \
         --extract $f \
         --make-bed \
         --out ${OUT}_${f/.txt/}
   done
   
   for rep in $(seq $nreps)
   do
      for f in snps_*.txt
      do
            OUT=$(basename $ROOT)_chr1_${f/.txt/}
      
            (time eval $FPCA \
               --scca \
               --bfile ${OUT} \
               --stand sd \
               --mem low \
               --pheno $expr \
               --lambda1 1e-3 \
               --lambda2 1e-3 \
               --ndim 1 \
               --v) 2>&1 | tee scca_${expr}_chr1_${f/.txt/}_rep${rep}.log
      done
   done
done

