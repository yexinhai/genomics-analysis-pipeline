# EVidenceModeler (EVM) genome annotation pipeline

by Xinhai Ye, yexinhai@zju.edu.cn, 2019.08.01

NOTE: For our new developed genome annotation pipeline, please see https://github.com/meiyang12/Genome-annotation-pipeline

## Repeat annotation
### RepeatModeler
```
ln -s ../01_data/01_Genome/genome.fa .
mkdir 01_repeatModeler-denovo-repeat.lib && cd 01_repeatModeler-denovo-repeat.lib
BuildDatabase -name GDB -engine ncbi ../genome.fa &>BuildDatabase_run.log

RepeatModeler -engine ncbi -pa 28 -database GDB &>RepeatModeler_run.log
cd ../
```
### RepeatMasker
```
mkdir 02_delete-denovo-lib-result 03_delete-repeatmasker-lib-result 04_delete-repeamasker-noint-result
RepeatMasker -lib 01_repeatModeler-denovo-repeat.lib/RM_*/consensi.fa.classified -pa 28 -dir 02_delete-denovo-lib-result genome.fa &>RepeatMasker_run.log1
RepeatMasker -pa 28 -dir 03_delete-repeatmasker-lib-result 02_delete-denovo-lib-result/genome.fa.masked &>RepeatMasker_run.log2
RepeatMasker -pa 28 -dir 04_delete-repeamasker-noint-result -noint 03_delete-repeatmasker-lib-result/genome.fa.masked.masked &>RepeatMasker_run.log3
```

## De novo gene prediction (use masked genome!!!)

### augustus
```
augustus --species=nasonia --gff3=on --uniqueGeneId=true --noInFrameStop=true --strand=both genome.fa > augustus.out
perl ./releated_scripts/Convert_Augustus.pl augustus.out > augustus.gff
```
### snap
```
export ZOE=/gpfs/bioinformatics/software/snap/Zoe
snap -gff /gpfs/bioinformatics/software/snap/HMM/Nasonia.hmm genoma.fa > snap.out
perl ./releated_scripts/Convert_SNAP nap.out > snap.gff
```

## Homology-based gene prediction

### genomethreader
```
gth -genomic genome.fa -protein invertebrate.all.fasta -gff3out -intermediate -o genomethreader.out
perl ./releated_scripts/Convert_gth.pl genomethreader.out > genomethreader.gff
```

### exonerate
```
exonerate --model protein2genome --percent 50 -t genome.fa -q invertebrate.all.fasta --showtargetgff yes --showalignment no --score 100 --minintron 20 --maxintron 20000 >exonerate.out
perl ./releated_scripts/Convert_exonerate.pl exonerate.out > exonerate.gff
```

### split a big genome into small piece files, then run augustus, genomethreader, and exonerate
```
perl split_fasta_by_multiple_methods-yxh.pl Pven.genome.fasta -m 1 -p 28 -d ./out

find . -name "Pven.genome.*" -exec bash -c 'echo bsub -n 1 -oo gth.log -eo gth \"gth -genomic $0 -protein invertebrate.all.fasta -gff3out -intermediate -o $0.out\"' {} \; >gth.cmd

less gth.cmd

chmod a+x gth.cmd

./gth.cmd

find . -name "Pven.genome.*" -exec bash -c 'echo bsub -n 1 -oo exonerate.log -eo exonerate \"exonerate --model protein2genome --percent 50 -t $0 -q invertebrate.all.fasta --showtargetgff yes --showalignment no --score 100 --minintron 20 --maxintron 20000 \>$0.out\"' {} \; >exonerate.cmd

less exonerate.cmd

chmod a+x exonerate.cmd

./exonerate.cmd

```

## Transcripts based gene prediction (use masked genome!!!)
Hisat2-Stringtie methods to get the merged.gtf, then use transdecoder to get gff3 format file.

### An example for Hisat2-Stringtie pipeline
```
hisat2-build -p 28 genome.fa genome &>bowtie2-build_run.log
cut -d \",\" -f 1 $RNA_Seq |sed 's/$/.gtf/' >mergelist.txt #prepare a list for stringtie-merge

#for LSF, if use cluster server
bsub -J PvenCA-Y1_RRAS22634-V.hisat2 -n 28 -oo PvenCA-Y1_RRAS22634-V.hisat2.log -eo PvenCA-Y1_RRAS22634-V.hisat2.log 'source /gpfs/bioinformatics/software/scriptByShowky/genomeAnnotation.path; hisat2 -p 28 -x genome --new-summary -S PvenCA-Y1_RRAS22634-V.sam -1 /gpfs/home/yexinhai/project/16_Pven_genome/01_data/02_rna_seq/PvenCA-Y1_RRAS22634-V_1.clean.fq -2 /gpfs/home/yexinhai/project/16_Pven_genome/01_data/02_rna_seq/PvenCA-Y1_RRAS22634-V_2.clean.fq'
bsub -J PvenCA-Y1_RRAS22634-V.samtools -n 28 -oo PvenCA-Y1_RRAS22634-V.samtools.log -eo PvenCA-Y1_RRAS22634-V.samtools.log -w 'done("PvenCA-Y1_RRAS22634-V.hisat2")' 'source /gpfs/bioinformatics/software/scriptByShowky/genomeAnnotation.path; samtools sort -@ 8 -o PvenCA-Y1_RRAS22634-V.bam PvenCA-Y1_RRAS22634-V.sam'
bsub -J PvenCA-Y1_RRAS22634-V.stringtie -n 28 -oo PvenCA-Y1_RRAS22634-V.stringtie.log -eo PvenCA-Y1_RRAS22634-V.stringtie.log -w 'done("PvenCA-Y1_RRAS22634-V.samtools")' 'source /gpfs/bioinformatics/software/scriptByShowky/genomeAnnotation.path; stringtie -p 28 -o PvenCA-Y1_RRAS22634-V.gtf PvenCA-Y1_RRAS22634-V.bam'

#for normal
hisat2 -p 28 -x genome --new-summary -S PvenCA-Y1_RRAS22634-V.sam -1 /gpfs/home/yexinhai/project/16_Pven_genome/01_data/02_rna_seq/PvenCA-Y1_RRAS22634-V_1.clean.fq -2 /gpfs/home/yexinhai/project/16_Pven_genome/01_data/02_rna_seq/PvenCA-Y1_RRAS22634-V_2.clean.fq
samtools sort -@ 8 -o PvenCA-Y1_RRAS22634-V.bam PvenCA-Y1_RRAS22634-V.sam
stringtie -p 28 -o PvenCA-Y1_RRAS22634-V.gtf PvenCA-Y1_RRAS22634-V.bam

stringtie --merge -p 28 -o stringtie.merged.gtf mergelist.txt &>stringtie.merged.log

```
### Transdecoder pipeline
```
#prepare genome and transcripts.
ln -s ../01_data/01_Genome/genome.fa .
ln -s ../04_stringtie/stringtie.merged.gtf transcripts.gtf

/gpfs/bioinformatics//software/TransDecoder-2.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl transcripts.gtf genome.fa >transcripts.fasta
/gpfs/bioinformatics//software/TransDecoder-2.0.1/util/cufflinks_gtf_to_alignment_gff3.pl transcripts.gtf >transcripts.gff3
/gpfs/bioinformatics//software/TransDecoder-2.0.1/TransDecoder.LongOrfs -t transcripts.fasta

blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep  -db /gpfs/bioinformatics//database/uniprot_20180118/uniprot_sprot.20180118.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 11 -out blastp.outfmt6 &>blastp_run.log
hmmscan --cpu 11 --domtblout pfam.domtblout /gpfs/bioinformatics//database/pfam.20180118/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep &>hmmscan_run.log

/gpfs/bioinformatics//software/TransDecoder-2.0.1/TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
/gpfs/bioinformatics//software/TransDecoder-2.0.1/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
```

## Running EVM
### Preparing gff3 file
```
cat augustus.gff snap.gff > gene_predictions.gff
perl ./releated_scripts/denovo_change_2_gff3.pl gene_predictions.gff >gene_predictions.gff3

cat genomethreader.gff exonerate.gff > protein_alignments.gff
perl ./releated_scripts/homolog_change_2_gff3.pl protein_alignments.gff >protein_alignments.gff3
```

### Use EVM script to check the gff3 file
```
perl gff3_gene_prediction_file_validator.pl your.gff3
```

### Preparing a weights file (for example)
```
ABINITIO_PREDICTION     AUGUSTUS        2
ABINITIO_PREDICTION     SNAP    1
PROTEIN nucleotide_to_protein_match     5
TRANSCRIPT      transdecoder    10
```

### Partitioning the inputs
```
perl /gpfs/bioinformatics/software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gff/gene_predictions.gff3 --protein_alignments gff/protein_alignments.gff3 --transcript_alignments gff/transcripts.fasta.transdecoder.genome.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
```

### Generating the EVM command set
```
perl /gpfs/bioinformatics/software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights /gpfs/home/yexinhai/project/16_Pven_genome/02_annotation/08_EVM/weights.txt --gene_predictions gff/gene_predictions.gff3 --protein_alignments gff/protein_alignments.gff3 --transcript_alignments gff/transcripts.fasta.transdecoder.genome.gff3 --output_file_name evm.out  --partitions partitions_list.out > commands.list
```

### Run EMV !!!
```
perl /gpfs/bioinformatics/software/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
```
or just
```
bash commands.list
```

### Combining the partitions
```
perl /gpfs/bioinformatics/software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
```

### Convert to gff3 format
```
perl /gpfs/bioinformatics/software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome genome.fasta
```

### Combine gff3 files
```
find ./ -mindepth 2 -maxdepth 2 -type f -name "evm.out.gff3" -exec cat {} \; >evm.gff3
```

## Get your genome annotation
### Rename gff3
```
python ./releated_scripts/gffrename.py evm.gff3 Pven > Pven.gff3
```

### Get CDS.fasta and PEP.fasta
```
perl ./releated_scripts/getGene.pl --posformat gff Pven.gff3 Pven.genome.fasta > Pven.cds.fasta

perl ./releated_scripts/cds2aa.pl Pven.cds.fasta > Pven.pep.fasta
```

### Get annotation information for each protein
```
blastp -query Pven.pep.fasta -db /gpfs/bioinformatics/database/uniprot_20180118/uniprot_sprot.20180118.fasta -num_threads 28 -evalue 1e-3 -outfmt 6 -max_target_seqs 1 -out Pven.fasta.blastpout6 &>blastp_run.log

perl ./releated_scripts/add_swissprot_annotation.pl /gpfs/bioinformatics/database/uniprot_20180118/uniprot_sprot.20180118.fasta Pven.pep.fasta Pven.fasta.blastpout6  >Pven.pep.anno.fasta
```

