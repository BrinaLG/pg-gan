#!/bin/bash

module load bcftools

# From github
#SUFFIX=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# for new Release
SUFFIX=.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

INPUT_CHR = /Users/blopezgf/Documents/Research/Project1/Data/1kG_20130502/
#INPUT_INDFILE = /Users/blopezgf/Documents/Research/Project1/Code/pg-gan/prep_data/IND_list/

CURR_dir = ${pwd}

# for each set of population(s)
for POP in CHB
do

  # for each chromosome
  for CHROM in `seq 1 22`
  do
      echo "bcftools view -S ${POP}_gan.txt --min-ac 1:minor -m2 -M2 -v snps -Oz -o ${POP}.chr${CHROM}${SUFFIX} ALL.chr${CHROM}${SUFFIX}"
      bcftools view -S ${POP}_gan.txt --min-ac 1:minor -m2 -M2 -v snps -Oz -o ${POP}.chr${CHROM}${SUFFIX} ALL.chr${CHROM}${SUFFIX}
  done

  # then merge into one vcf
  echo "bcftools concat -f ${POP}_filelist.txt -Oz -o ${POP}${SUFFIX}"
  bcftools concat -f ${POP}_filelist.txt -Oz -o ${POP}${SUFFIX}
done
