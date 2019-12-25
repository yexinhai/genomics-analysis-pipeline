## calling SNP 
Ye xinhai,2019.12.24

### 1.Create index for reference genome：
```
samtools faidx genome.fa
bwa index genome.fa
```

### 2.Mapping
```
bwa mem -t 28 -M genome.fa R1.fastq  R2.fastq > xxx.sam
```

Also, we can use bowtie2
```
bowtie2-bulid genome.fa genome

bowtie2 -p 28 -x genome -1 R1.fastq  -2 R2.fastq -S xxx.sam
```

### 3.Sam file conversion and bam file sorting
```
samtools view -b -S xxx.sam > xxx.bam

samtools sort xxx.bam -o xxx.sorted.bam

Mark Duplications using picard (optional?!,not sure)

samtools index xxx.sorted.bam
```

### Prepare data for SNP calling, can see details in https://www.plob.org/article/7009.html (not only for GATK)


### 4.SNP calling (we use freebayes here!)
Our default sample is diploid, use --ploidy 1, if your sample is haploid.
```
freebayes -f genome.fa -b xxx.sorted.bam >freebayes.vcf

```

### 5.VCF file filtering (ref: http://ddocent.com/filtering/)
```
vcffilter -f "QUAL > 20" freebayes.vcf > freebayes_qual20.vcf

#filtering out loci with an allele balance below 0.25 and above 0.75
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" freebayes_qual20.vcf > freebayes_qual20_AB.vcf

#filtering out sites that have reads from both strands. (optional, useful in GWAS or RNA-seq)
vcffilter -f "SAF / SAR > 100 whether or not their is a discrepancy in the properly paired status of for reads supporting reference or alternate alleles& SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s freebayes_qual20_AB.vcf > freebayes_qual20_AB_SBS.vcf

#ratio of mapping qualities between reference and alternate alleles
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" freebayes_qual20_AB_SBS.vcf > freebayes_qual20_AB_SBS_MQ.vcf

#whether or not their is a discrepancy in the properly paired status of for reads supporting reference or alternate alleles
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s freebayes_qual20_AB_SBS_MQ.vcf > freebayes_qual20_AB_SBS_MQ_PAIRED.vcf

#removing any locus that has a quality score below 1/4 of the depth
vcffilter -f "QUAL / DP > 0.25" freebayes_qual20_AB_SBS_MQ_PAIRED.vcf > freebayes_qual20_AB_SBS_MQ_PAIRED_DP.vcf

vcfallelicprimitives freebayes_qual20_AB_SBS_MQ_PAIRED_DP.vcf --keep-info --keep-geno > freebayes_qual20_AB_SBS_MQ_PAIRED_DP_prim.vcf

vcftools --vcf freebayes_qual20_AB_SBS_MQ_PAIRED_DP_prim.vcf  --remove-indels --recode --recode-INFO-all --out freebayes_qual20_AB_SBS_MQ_PAIRED_DP_prim_SNP.vcf
```