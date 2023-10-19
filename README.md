# metadata
Genome list (S. brevipes, S. brunnescens, S. luteus) : metadata/genome_list.txt  
Genome list (S. luteus only): metadata/s208_list.txt  
Family (Clade) assignments for plink: metadata/pop.fam.txt  
Individual to exclude according to Sumpplementary Appendix A: /metadata/clone.fam.txt, /metadata/hap.fam.txt, /metadata/mixed.fam.txt, /metadata/self.fam.txt  
Raw sequences: raw_fastq/  
Reference genome: ge/Suilu4.fasta  


# subsample deep reads
```bash
seql=$(cat metadata/genome_list.txt) #all samples
for seq in $seql
do 
  R1path=raw_fastq/${seq}_R1.fastq.gz
  R2path=raw_fastq/${seq}_R2.fastq.gz
  R1outpath=raw_fastq/${seq}_R1p.fastq.gz
  R2outpath=raw_fastq/${seq}_R2p.fastq.gz
  #subset reads
  seqtk sample -2 -s 100 $R1path 10000000 | gzip -c > $R1outpath
  seqtk sample -2 -s 100 $R2path 10000000 | gzip -c > $R2outpath
done
```

# trim adaptors
```bash
mkdir cutadapt_reports

for seq in $seql
do
    R1inpath=raw_fastq/${seq}_R1p.fastq.gz
    R2inpath=raw_fastq/${seq}_R2p.fastq.gz
    R1outpath=raw_fastq/${seq}_R1_trim.fastq.gz
    R2outpath=raw_fastq/${seq}_R2_trim.fastq.gz
    cutadapt -a  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 25 \
        -o $R1outpath -p $R2outpath \
        -m 80 \
        $R1inpath $R2inpath > cutadapt_reports/${seq}.log
done #min quality 25; min length 80bp
```

# mapping
## BWA
```bash
mkdir bwa_mapping
bwa index ge/Suilu4.fasta

seql=$(cat metadata/genome_list.txt)
for seq in $seql
do 
  R1path=raw_fastq/${seq}_R1_trim.fastq.gz
  R2path=raw_fastq/${seq}_R2_trim.fastq.gz
  outpath=bwa_mapping/${seq}.sam
  logpath=bwa_mapping/${seq}.log
  bwa mem ge/Suilu4.fasta $R1path $R2path > $outpath
done

seql=$(cat metadata/genome_list.txt)
for seq in $seql
do
  samfile=bwa_mapping/${seq}.sam
  bamfile=bwa_mapping/${seq}_sorted.bam
  logfile=bwa_mapping/${seq}_convert.log
  samtools view -F 4 -h $samfile | samtools sort -O BAM -o $bamfile -
done
```

## check depth
```bash
mkdir bwa_depth

seql=$(cat metadata/genome_list.txt)
for seq in $seql
do 
  bamfile=bwa_mapping/${seq}_sorted.bam
  outpath=bwa_depth/${seq}_depth.txt
  samtools depth $bamfile > $outpath
done

for seq in $seql
do 
  infile=bwa_depth/${seq}_depth.txt
  outpath=bwa_depth/${seq}_depth_avg.txt
  awk_arg=\'"{ sum += \$3; n++ } END { if (n > 0) print \"$seq\t\" sum / n}"\'
  awk $awk_arg $infile > $outpath
done

rm bwa_depth/all_depth_avg.txt
cat bwa_depth/*_depth_avg.txt > bwa_depth/all_depth_avg.txt
```

# call variants
## find repeat regions
```bash
cat ge/Suilu4.fasta | perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (sort(keys(%fa))){ print "$s\n$fa{$s}\n" }}' | perl -lne 'if(/^>(\S+)/){ $n=$1} else { while(/([a-z]+)/g){ printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_)) } }' > ge/Suilu4.mask.bed #get masked region shown in lowercase in the reference genome

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ge/Suilu4.mask.bed #print total length masked
```

## thin deep genomes
```bash
mkdir bwa_haplotypecaller_gvcf bwa_haplotypecaller_log
samtools faidx ge/Suilu4.fasta
picard CreateSequenceDictionary R=ge/Suilu4.fasta O=ge/Suilu4.dict

seql=$(cat metadata/genome_list.txt)
for seq in $seql
do
  tcov=30
  inbam="bwa_mapping/${seq}_sorted.bam"
  outbam="bwa_mapping/${seq}_thin.bam"
  cov=$(grep ${seq}$'\t' bwa_depth/all_depth_avg.txt | cut -d$'\t' -f 2)
  frac=$(bc <<< "scale=6; ${tcov}/${cov}")
  if (( $(echo "$frac < 1" | bc -l) )); then
    samtools view -h -s $frac -b $inbam > $outbam
    else samtools view -h -b $inbam > $outbam
  fi
done
```

## GATK haplotypecaller
```bash
#fix info for each sample
seql=$(cat metadata/genome_list.txt)
for seq in $seql
  do picard AddOrReplaceReadGroups I=bwa_mapping/${seq}_thin.bam  O=bwa_mapping/${seq}_info.bam RGLB=lib1 RGPL=illumina RGPU=unit1  RGSM=${seq}
done

for seq in $seql
do
  picard BuildBamIndex INPUT=bwa_mapping/${seq}_info.bam
done

#PL for each samples
for seq in $seql
do
  inbam=bwa_mapping/${seq}_info.bam
  outvcf=bwa_haplotypecaller_gvcf/${seq}.g.vcf.gz
  gatk --java-options '-Xmx8g' HaplotypeCaller \
  -R ge/Suilu4.fasta \
  --native-pair-hmm-threads 1 \
  -I $inbam -O $outvcf -ploidy 2 -ERC GVCF
done

#combine samples to whole .vcf
vflag275=$(cat metadata/genome_list.txt | awk '{print "--variant bwa_haplotypecaller_gvcf/" $1 ".g.vcf.gz "}' | xargs echo)
vflag208=$(cat metadata/s208_list.txt | awk '{print "--variant bwa_haplotypecaller_gvcf/" $1 ".g.vcf.gz "}' | xargs echo)
contigs=$(cat ge/Suilu4.fasta.fai | awk '{print $1}')

for chr in $contigs
do
  gatk CombineGVCFs -R ge/Suilu4.fasta $vflag275 -XL ge/Suilu4.mask.bed -L $chr -O bwa_haplotypecaller_gvcf/s275_${chr}.g.vcf.gz
  gatk CombineGVCFs -R ge/Suilu4.fasta $vflag208 -XL ge/Suilu4.mask.bed -L $chr -O bwa_haplotypecaller_gvcf/s208_${chr}.g.vcf.gz
done

#paralleling by contigs
contigs=$(cat ge/Suilu4.fasta.fai | awk '{print $1}')
for chr in $contigs
do
  gatk GenotypeGVCFs -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s275_${chr}.g.vcf.gz -O bwa_haplotypecaller_gvcf/s275_${chr}.vcf.gz
  gatk GenotypeGVCFs -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s208_${chr}.g.vcf.gz -O bwa_haplotypecaller_gvcf/s208_${chr}.vcf.gz
done
```

## filter SNP
```bash
# hardfilter
filter_exp='--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'
contigs=$(cat ge/Suilu4.fasta.fai | awk '{print $1}')
for chr in $contigs
do
  gatk VariantFiltration -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s275_${chr}.vcf.gz -O bwa_haplotypecaller_gvcf/s275_${chr}_hardfilter.vcf.gz --filter-name hard_filters $filter_exp
done

for chr in $contigs
do
  gatk VariantFiltration -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s208_${chr}.vcf.gz -O bwa_haplotypecaller_gvcf/s208_${chr}_hardfilter.vcf.gz --filter-name hard_filters $filter_exp
done

#SNPs only
contigs=$(cat ge/Suilu4.fasta.fai | awk '{print $1}')
for chr in $contigs
do
  gatk SelectVariants -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s275_${chr}_hardfilter.vcf.gz -O bwa_haplotypecaller_gvcf/s275_${chr}_filtered.vcf.gz --exclude-filtered --select-type-to-include SNP
  gatk SelectVariants -R ge/Suilu4.fasta -V bwa_haplotypecaller_gvcf/s208_${chr}_hardfilter.vcf.gz -O bwa_haplotypecaller_gvcf/s208_${chr}_filtered.vcf.gz --exclude-filtered --select-type-to-include SNP

for chr in $contigs
do
  python vcf_minitools/strict_SNP.py -i bwa_haplotypecaller_gvcf/s275_${chr}_filtered.vcf.gz -o bwa_haplotypecaller_gvcf/s275_${chr}_SNPfiltered.vcf.gz
  python vcf_minitools/strict_SNP.py -i bwa_haplotypecaller_gvcf/s208_${chr}_filtered.vcf.gz -o bwa_haplotypecaller_gvcf/s208_${chr}_SNPfiltered.vcf.gz
done

#combine vcf files
mkdir bwa_haplotypecaller_finalvcf
ls bwa_haplotypecaller_gvcf/s275_scaffold_*_SNPfiltered.vcf.gz > bwa_haplotypecaller_gvcf/s275_SNPfiltered.list.txt
ls bwa_haplotypecaller_gvcf/s208_scaffold_*_SNPfiltered.vcf.gz > bwa_haplotypecaller_gvcf/s208_SNPfiltered.list.txt

gatk MergeVcfs -I bwa_haplotypecaller_gvcf/s275_SNPfiltered.list.txt -O bwa_haplotypecaller_finalvcf/s275_SNPfiltered.vcf.gz
gatk MergeVcfs -I bwa_haplotypecaller_gvcf/s208_SNPfiltered.list.txt -O bwa_haplotypecaller_finalvcf/s208_SNPfiltered.vcf.gz

vcftools --depth --gzvcf bwa_haplotypecaller_finalvcf/s275_SNPfiltered.vcf.gz --out bwa_haplotypecaller_finalvcf/s275
vcftools --depth --gzvcf bwa_haplotypecaller_finalvcf/s208_SNPfiltered.vcf.gz --out bwa_haplotypecaller_finalvcf/s208
```

# test AD
```bash
mkdir AF
IID=$(cat metadata/genome_list.txt)

for i in $IID
do
  scripts/bcftools_print_AD.sh $i
done

Rscript scripts/AD_test.R
```

# subset individuals
```bash
#calculate heterozygote percent
plink --allow-extra-chr --vcf bwa_haplotypecaller_finalvcf/s275_SNPfiltered.vcf.gz --out plink/s275 --make-bed --set-missing-var-ids @:# --update-ids metadata/pop.fam.txt --double-id
plink --allow-extra-chr --bfile plink/s275 --out plink/s275 --geno 0.1 --het

# subset s225: (no mixed, no haploid, no self)
plink --allow-extra-chr --bfile plink/s275 --out plink/s224 --remove <(cat metadata/mixed.fam.txt metadata/hap.fam.txt metadata/self.fam.txt) --mac 1 --make-bed

# check clone
plink --allow-extra-chr --bfile plink/s224 --out plink/s224_con --geno 0.1 --distance square 1-ibs

# subset s218 (no mixed, no haploid, no self, no clone)
plink --allow-extra-chr --bfile plink/s224 --out plink/s214 --remove metadata/clone.fam.txt --mac 1 --make-bed

# subset s214 (no clone, no mixed, no haploid, no self): phylogenetic tree
plink --allow-extra-chr --bfile plink/s214 --out vcf_fasta/s214 --recode vcf-iid

# subset s213 (no clone, no mixed, no haploid, no self; S. luteus and S. brunnesens): Reynold's and Fst
plink --allow-extra-chr --bfile plink/s214 --out plink/s213 --family --remove-cluster-names Sbr --mac 1 --make-bed

## LD prune
plink --allow-extra-chr --bfile plink/s213 --out plink/s213_admixture --biallelic-only --mac 2 --indep-pairwise 50 10 0.1 --geno 0.1 
plink --allow-extra-chr --bfile plink/s213 --out plink/s213_admixture --extract plink/s213_admixture.prune.in --make-bed # biallelic
plink --allow-extra-chr --bfile plink/s213_admixture --out plink/s213_admixture --recode A

# subset s208 (no clone, no mixed, no haploid; S. luteus only): pca, admixture
plink --allow-extra-chr --vcf bwa_haplotypecaller_finalvcf/s208_SNPfiltered.vcf.gz --out plink/s208 --make-bed --set-missing-var-ids @:# --update-ids metadata/pop.fam.txt --double-id
plink --allow-extra-chr --bfile plink/s208 --out plink/s208 --recode A
plink --allow-extra-chr --bfile plink/s208 --out plink/s208 --recode

## LD prune
plink --allow-extra-chr --bfile plink/s208 --out plink/s208_admixture --biallelic-only --mac 2 --indep-pairwise 50 10 0.1 --geno 0.1 
plink --allow-extra-chr --bfile plink/s208 --out plink/s208_admixture --extract plink/s208_admixture.prune.in --make-bed # biallelic

plink --allow-extra-chr --bfile plink/s208_admixture --out plink/s208_admixture --recode A
plink --allow-extra-chr --bfile plink/s208_admixture --out plink/s208_admixture --recode

## avoid missing site: genetic diversity (Nucleotide diversity, Watterson's theta, Tajima's D)
plink --allow-extra-chr --bfile plink/s208 --out plink/s208_con --mac 1 --make-bed --geno 0.1
plink --allow-extra-chr --bfile plink/s208_con --out plink/s208_con --recode A

# s189 (no clone, no mixed, no haploid; Slu st): splittree
plink --allow-extra-chr --bfile plink/s208 --out plink/s189 --family --mac 1 --remove-cluster-names Slu1 Slu2 --make-bed
plink --allow-extra-chr --bfile plink/s189 --out plink/s189_admixture --biallelic-only --mac 2 --indep-pairwise 50 10 0.1 --geno 0.1
plink --allow-extra-chr --bfile plink/s189 --out plink/s189_admixture --extract plink/s189_admixture.prune.in --make-bed 
plink --allow-extra-chr --bfile plink/s189 --out plink/s186  --recode A

```

# Make tree
```bash
mkdir vcf_fasta
python scripts/vcf2phylip.py -i vcf_fasta/s214.vcf  -o vcf_fasta/s214.phy --remove-raxml-invar --skip-check --recode-vcf vcf_fasta/s214_phy.vcf.gz
mpiexec -bynode -np 20 raxmlHPC-HYBRID-SSE3 -T 8 -n s214-CAT -f a -s s214.phy -m GTRCAT -p 12345 -x 12345 -#1000 -V
```

# Kinship 
```bash
mkdir kinship
cp plink/s208.bed plink/s208.fam plink/s208.bim kinship
sed -i 's/scaffold_//g' kinship/s208.bim
king -b kinship/s208.bed --prefix ./kinship/s208 --kinship --sexchr 67
```

# Admixture
```bash
mkdir admixture_s208
cp plink/s208_admixture.bed plink/s208_admixture.bim plink/s208_admixture.fam plink/s208_admixture.map admixture_s208

cd admixture_s20
sed -i 's/scaffold_//g' s208_admixture.bim

mkdir seed1
cd seed1
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed2
cd seed2
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 12345 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed3
cd seed3
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 54321 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed4
cd seed4
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 0 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed5
cd seed5
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 1 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed6
cd seed6
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 2 --cv ../s208_admixture.bed $K"
done
cd ../


mkdir seed7
cd seed7
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 3 --cv ../s208_admixture.bed $K"
done
cd ../


mkdir seed8
cd seed8
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 4 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed9
cd seed9
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 5 --cv ../s208_admixture.bed $K"
done
cd ../

mkdir seed10
cd seed10
for K in {1..8}
do
sbatch -J admixture -c 4 --mem=16g -o s208_K${K}.log --wrap="admixture -s 6 --cv ../s208_admixture.bed $K"
done
cd ../
```

# G statistic, Fst, Reynold's distance, nucleotide diversity, Watterson's theta, Tajima's D
```bash
Rscript scripts/within_gstats.R 10000 #subset 10k SNPs
Rscript scripts/Fst_hier.R
Rscript scripts/fst.R
Rscript scripts/tajimaD.R 10000 #win=10000
Rscript scripts/reynold.R
```

# SplitsTree
```bash
plink --allow-extra-chr --bfile plink/s189_admixture --out plink/s189_admixture --recode vcf-iid bgz
python ./vcf2phylip/vcf2phylip.py -i plink/s189_admixture.vcf.gz -o splitstree/s189_admixture -p -n
```
