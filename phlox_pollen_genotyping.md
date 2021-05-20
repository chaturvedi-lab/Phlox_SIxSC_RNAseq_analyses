
**STEP 1: Calling variants using samtools mpileup**
I used samtools mpileup to call variants for to identify genotypes for pollen genes identification. For doing this, we need the sorted bam files created after alignment of the data to the reference genome. I used the following bash script for variant calling. I have specified the flags I used for samtools mpileup which can be adjusted as and when required. 

```bash
#!/bin/bash
#SBATCH -p shared # Partition to submit to
#SBATCH -n 16     # Number of cores
#SBATCH -N 1
#SBATCH -t 3-0:00  # Runtime in days-hours:minutes
#SBATCH --mem 100000
#SBATCH -J varcall
#SBATCH -o varcall_%A.out
#SBATCH -e varcall_%A.err
#SBATCH --mail-type=ALL
notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu

    echo ------------------------------------------------------
    echo SLURM: job identifier is $SLURM_JOBID
    echo SLURM: job name is $SLURM_JOB_NAME
    echo ------------------------------------------------------

#SAMTOOLs mpileup version 1.5 options used:
#C = adjust mapping quality; recommended:50, disable:0 [0]
#d = max per-file depth; avoids excessive memory usage [250]
#f = faidx indexed reference sequence file
#q = skip alignments with mapQ smaller than INT [0]
#Q = skip bases with baseQ/BAQ smaller than INT [13]
#g = generate genotype likelihoods in BCF format
#t = --output-tags LIST  optional tags to output:DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []

#BCFTOOLs call version 1.6 options used
#v = output variant sites only
#c/m = he original calling method (conflicts with -m) or alternative model for multiallelic and rare-variant calling (conflicts with -c)
#p = variant if P(ref|D)<FLOAT with -c [0.5]
#P =  --prior <float-o mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
#O =  output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' (here it is 'v')
#o = write output to a file [standard output]

module load samtools
module load bcftools

cd /n/holyscratch01/hopkins_lab/Chaturvedi/rnaseq_analyses/pollen_genotype/bamfiles/

samtools mpileup -C 50 -d 250 -f /n/holyscratch01/hopkins_lab/Chaturvedi/rnaseq_analyses/pollen_genotype/genome_assembly/assembly.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b phloxBam.txt -o variantsPhlox.bcf

bcftools call -v -c -p 0.01 -P 0.001 -O v -o variantsPhlox.vcf variantsPhlox.bcf
```

**STEP 2: Creating genotypes file using vcftools**
After I created the variants (vcf) file, I wanted to extract the genotypes called for each individual for each variant in the file. To get the genotype file, I used vcftools as follows:

```bash
ml vcftools
vcftools --vcf variantsPhlox.vcf --extract-FORMAT-info GT --out genotypes_file
```

**STEP 3: Comparing genotypes between individuals using R**
Once I got the genotype file, it was easy to compare the genotypes between individuals. I did this in R and identified the contigs of interest and wrote them out in a separate file. Here is the code I used to do this. This is example code which can be modified based on sample IDs. Additionally, I have tested this with a small set of loci but if there are a lot of variants, it would make sense to run this code on the cluster using a SLURM job submission script like above and just calling "R CMD mycode.R". mycode.R would be a possible R script with all the code below.


```R
#read in the genotype file
gtdat<-read.table("genotypes.txt", header=T)

#subset individuals from a single mom to make genotype comparisons; here it is all individuals with A in their sample ID. #AIN is the sample with immature pistil and no pollen treatment. I will compare this with AIS which is the sample with immature pistil and self pollen treatment. AIN should have pistil genes and AIS should have some pollen genes. These will differ between the two treatments. Therefore, I will identify all the loci where AIN and AIS have different genotypes.

ais<-cbind(gtdat[,c(1:2)],gtdat$sorted_Henry_AIS_S3_hisat2.bam)
ain<-cbind(gtdat[,c(1:2)],gtdat$sorted_Henry_AIN_S5_hisat2.bam)

#simulating missing genotypes to make sure we dont count those in our final comparison. This will already exist in the vcf file. This is just to show the example code.
#I am adding missing to AIS sample and not to AIN
miss<-rep("./.",4)
miss_con<-rep("contig_100",4)
miss_pos<-c(seq(1:4))
ais_miss<-cbind(miss_con, miss_pos, miss)
colnames(ais_miss)<-colnames(ais)

#final dataframe for ais
ais_m<-rbind(ais,ais_miss)

#for AIN we will simulate some genotypes which are not missing
ain_n<-ais_miss
colnames(ain_n)<-colnames(ain)
ain_n[,3]<-c("0/1","1/1","0/0","0/1")

#final dataframe for ais
ain_nm<-rbind(ain, ain_n)

#identify the missing genotypes in both samples and then drop those variants
miss_gt<-which(ais_m[,3] == "./." | ain_nm[,3] == "./.")

#remove those variants/rows from the ain and ais samples
ais_m<-ais_m[-miss_gt,]
ain_nm<-ain_nm[-miss_gt,]


#get the set of genotypes which differ between the two
diff_gt_ais_ain<-ais_m[which(ais_m[,3] != ain_nm[,3]),]

#write out the chromosome, position from this list to 
write.table(diff_gt_ais_ain[,c(1:2)], "pollen_genotype_vars.txt", col.names=T, row.names=F, sep =" ", quote = F)

```

**STEP 4: Associating variants with genes**
The variants from the list above have a chromosome ID and a position on the chromosome where the variant maps. However, for the actual gene expression analyses, ideally we should have a transcript ID with a start and stop position on the genome. Therefore, we need to find the gene on which the variant lies to conclude that it is a pollen gene.

I used the following script associate the variants with the genes: create_snp_annotations.py. We can run this script on the start or stop position of the transcripts to make sure that the variant is in the correct range of the positions.

I used only the scaffold and start or scaffold and stop position from the scafs_transcriptids.txt file to run this annotation.

python create_snp_annotations.py --map transcripts_start --ann genome_annotation.txt --out out_transcripts_start
python create_snp_annotations.py --map transcripts_stop --ann genome_annotation.txt --out out_transcripts_stop












