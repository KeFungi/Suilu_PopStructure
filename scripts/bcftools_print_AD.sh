#!/bin/bash
i=$1
bcftools view -m2 -M2 -H -s ${i} bwa_haplotypecaller_finalvcf/s275.vcf.gz | \
awk '{print $10}' | awk -F: '{print $1 OFS $2 OFS $3}' | awk -F, '{print $1 OFS $2}' > AF/${i}.AD.txt
