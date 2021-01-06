Henry RNA sequence - Analyses

This is a pipeline to run RNA sequence analyses. I start with filtering of reads, adapter trimming, filtering of transcripts based on annotation and then preparing a gene expression count matrix for differential expression analyses.

**STEP 1 Quality check, k-mer removal and quality check**

1. Ran FastQC for quality checking of the reads

```bash
#!/usr/bin/env/bash

#SBATCH -p bigmem    # Partition to submit to
#SBATCH -n 1         # Number of cores
#SBATCH -t 0-8:00    # Runtime in days-hours:minutes
#SBATCH --mem 200000   # Memory in MB
#SBATCH -J FastQC    # job name
#SBATCH -o FastQC.%A.out # File to which standard out will be written
#SBATCH -e FastQC.%A.err # File to which standard err will be written
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<schaturvedi@fas.harvard.edu>  # Email to which notifications will be sent



module purge
module load fastqc/0.11.5-fasrc01

for forward in /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/raw_sequences/*R1_001.fastq; do
        #prefix="${forward%R1_001.fastq}";
        reverse="${forward%R1_001.fastq}R2_001.fastq";
        #echo "$prefix $forward $reverse 12";
        fastqc -1 $forward -2 $reverse -t 5 -o /n/holyscratch01/hopkins_lab/Chaturvedi/fastqc/;
done

```

2. Ran rcorrector to identify erroneous k-mers in transcripts

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared                   # may consider running on a bigmem node for large dataset
#SBATCH -e rcorrector_%A.err            # File to which STDERR will be written
#SBATCH -o rcorrector_%A.out           # File to which STDOUT will be written
#SBATCH -J rcorrector_%A               # Job name
#SBATCH --mem=150000                 # Memory requested
#SBATCH --time=7-00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to


module purge
module load Rcorrector/20180919-fasrc01

for prefix in $(ls /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/raw_sequences/*.fastq | sed -r 's/_R[12]_001[.]fastq//' | uniq)

do

perl /n/helmod/apps/centos7/Core/Rcorrector/20180919-fasrc01/bin/run_rcorrector.pl -t 12 -1 "${prefix}_R1_001.fastq" -2 "${prefix}_R2_001.fastq" -od /n/holyscratch01/hopkins_lab/Chaturvedi/rcorrector/

done
```

After running Rcorrector I ended up with files where reads were marked "cor" or "unfixable". I then used a python script (FilterUncorrectabledPEfastq.py) to remove these "unfixable" reads from the read pairs.  This created log files for each pair of reads (51 files) and files for each read with prefix "rmunfixable" (102 files).
I then wanted to get a stats for total reads and removed reads. The log files give stats for total reads for both reads and total removed reads for both reads. I did this in bash and created two files:

Here is the python script:

```python
"""
author: adam h freedman
afreedman405 at gmail.com
data: Fri Aug 26 10:55:18 EDT 2016
This script takes as an input Rcorrector error corrected Illumina paired-reads
in fastq format and:
1. Removes any reads that Rcorrector indentifes as containing an error,
but can't be corrected, typically low complexity sequences. For these,
the header contains 'unfixable'.
2. Strips the ' cor' from headers of reads that Rcorrector fixed, to avoid
issues created by certain header formats for downstream tools.
3. Write a log with counts of (a) read pairs that were removed because one end
was unfixable, (b) corrected left and right reads, (c) total number of
read pairs containing at least one corrected read.
Currently, this script only handles paired-end data, and handle either unzipped
or gzipped files on the fly, so long as the gzipped files end with 'gz'.
"""

import sys
import gzip
from itertools import izip,izip_longest
import argparse
from os.path import basename

def get_input_streams(r1file,r2file):
    if r1file[-2:]=='gz':
        r1handle=gzip.open(r1file,'rb')
        r2handle=gzip.open(r2file,'rb')
    else:
        r1handle=open(r1file,'r')
        r2handle=open(r2file,'r')

    return r1handle,r2handle


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector fastq outputs")
    parser.add_argument('-1','--left_reads',dest='leftreads',type=str,help='R1 fastq file')
    parser.add_argument('-2','--right_reads',dest='rightreads',type=str,help='R2 fastq file')
    parser.add_argument('-s','--sample_id',dest='id',type=str,help='sample name to write to log file')
    opts = parser.parse_args()

    r1out=open('unfixrm_%s' % basename(opts.leftreads).replace('.gz',''),'w')
    r2out=open('unfixrm_%s' % basename(opts.rightreads).replace('.gz','') ,'w')

    r1_cor_count=0
    r2_cor_count=0
    pair_cor_count=0
    unfix_r1_count=0
    unfix_r2_count=0
    unfix_both_count=0

    r1_stream,r2_stream=get_input_streams(opts.leftreads,opts.rightreads)

    with r1_stream as f1, r2_stream as f2:
        R1=grouper(f1,4)
        R2=grouper(f2,4)
        counter=0
        for entry in R1:
            counter+=1
            if counter%100000==0:
                print "%s reads processed" % counter

            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
            head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]

            if 'unfixable' in head1 and 'unfixable' not in head2:
                unfix_r1_count+=1
            elif 'unfixable' in head2 and 'unfixable' not in head1:
                unfix_r2_count+=1
            elif 'unfixable' in head1 and 'unfixable' in head2:
                unfix_both_count+=1
            else:
                if 'cor' in head1:
                    r1_cor_count+=1
                if 'cor' in head2:
                    r2_cor_count+=1
                if 'cor' in head1 or 'cor' in head2:
                    pair_cor_count+=1

                head1=head1.split('l:')[0][:-1]
                head2=head2.split('l:')[0][:-1]
                r1out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))
                r2out.write('%s\n' % '\n'.join([head2,seq2,placeholder2,qual2]))

    total_unfixable = unfix_r1_count+unfix_r2_count+unfix_both_count
    total_retained = counter - total_unfixable

    unfix_log=open('rmunfixable_%s.log' % opts.id,'w')
    unfix_log.write('total PE reads:%s\nremoved PE reads:%s\nretained PE reads:%s\nR1 corrected:%s\nR2 corrected:%s\npairs corrected:%s\nR1 unfixable:%s\nR2 unfixable:%s\nboth reads unfixable:%s\n' % (counter,total_unfixable,total_retained,r1_cor_count,r2_cor_count,pair_cor_count,unfix_r1_count,unfix_r2_count,unfix_both_count))

    r1out.close()
    r2out.close()
    unfix_log.close()

```

```bash
grep ^total rmunfixable_Henry_*.log | cut -d':' -f3 > totalreads.txt
grep ^both rmunfixable_Henry_*.log | cut -d':' -f3 > removedreads.txt
```

3. Ran TrimGalore! for adapter trimming

```bash

#!/bin/bash
#SBATCH -J trimgalore
#SBATCH -n 16                     # Use 1 cores for the job
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 6-00:00                 # Runtime in D-HH:MM
#SBATCH -p shared         # Partition to submit to
#SBATCH --mem=100000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o trimgalore_PE.%A.out  # File to which STDOUT will be written
#SBATCH -e trimgalore_PE.%A.err  # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to



module purge
module load TrimGalore/0.5.0-fasrc01
module load python/3.6.3-fasrc02
source activate cutadapt



# $1 = R1 reads

# $2 = R2 reads



for prefix in $(ls unfixrm_Henry_*.cor.fq | sed -r 's/_R[12]_001[.]cor.fq//' | uniq)

do

trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 "${prefix}_R1_001.cor.fq" "${prefix}_R2_001.cor.fq"

done

```

**STEP 2 Functional annotation of transcripts and dropping transcripts with lower hits from the analyses**

1. Transdecoder (read more here: https://github.com/TransDecoder/TransDecoder/wiki)
TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments to the genome using Tophat and Cufflinks.

TransDecoder identifies likely coding sequences based on the following criteria:

a. a minimum length open reading frame (ORF) is found in a transcript sequence

b. a log-likelihood score similar to what is computed by the GeneID software is > 0.

the above coding score is greatest when the ORF is scored in the 1st reading frame as compared to scores in the other 2 forward reading frames.

c. if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF, the longer one is reported. However, a single transcript can report multiple ORFs (allowing for operons, chimeras, etc).

d. a PSSM is built/trained/used to refine the start codon prediction.

e. optional the putative peptide has a match to a Pfam domain above the noise cutoff score.

Here is the bash script I used to run Transdecoder on the cluster:

```bash

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared                   # may consider running on a bigmem node for large dataset
#SBATCH -e transdec_%A.err            # File to which STDERR will be written
#SBATCH -o transdec_%A.out           # File to which STDOUT will be written
#SBATCH -J transdec_%A               # Job name
#SBATCH --mem=150000                 # Memory requested
#SBATCH --time=7-00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to



module load TransDecoder/5.3.0-fasrc01

#Running TransDecoder is a two-step process. First run the TransDecoder step that identifies all long ORFs.

$TRANSDECODER_HOME/TransDecoder.LongOrfs -t /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trinity/Trinity_assembly/Trinity.fasta

#Now, run the step that predicts which ORFs are likely to be coding.

$TRANSDECODER_HOME/TransDecoder.Predict -t /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trinity/Trinity_assembly/Trinity.fasta

```

You'll now find a number of output files containing 'transdecoder' in their name:

```bash
ls -1 |grep transdecoder
```

```
Trinity.fasta.transdecoder.bed
Trinity.fasta.transdecoder.cds
Trinity.fasta.transdecoder.gff3
Trinity.fasta.transdecoder.mRNA
Trinity.fasta.transdecoder.pep    
Trinity.fasta.transdecoder_dir/
```

The file we care about the most here is the 'Trinity.fasta.transdecoder.pep' file, which contains the protein sequences corresponding to the predicted coding regions within the transcripts.

2. Sequence homology searches using transdecoder output with blastx and blastp

Earlier, I ran blastx against our mini SWISSPROT phlox database to identify likely full-length transcripts. I am running blastx again to capture likely homolog information, and I'll lower our E-value threshold to 1e-5 to be less stringent than earlier.

First create a database for the pep file created by transdecoder

```
makeblastdb -in ./Trinity.fasta.transdecoder.pep -dbtype prot
```

Then run blastx using the transdecoder output.

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared                   # may consider running on a bigmem node for large dataset
#SBATCH -e blastx_%A.err            # File to which STDERR will be written
#SBATCH -o blastx_%A.out           # File to which STDOUT will be written
#SBATCH -J blastx_%A               # Job name
#SBATCH --mem=150000                 # Memory requested
#SBATCH --time=7-00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to


module load blast/2.6.0+-fasrc01

blastx -db ../transdecoder/Trinity.fasta.transdecoder.pep \
         -query /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trinity/Trinity_assembly/Trinity.fasta \ -num_threads 12 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > td_swissprot.blastx.outfmt6
```

Now, let's look for sequence homologies by just searching our predicted protein sequences rather than using the entire transcript as a target. To do this I ran blastp. For this I first downloaded the uniprot_sprot fasta file as follows:

```bash
mkdir uniprot_sprot

URL=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
URL=ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz

wget -b $URL
wget -nv -P data/ $URL
gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta

#make database of the fasta for next step
makeblastdb -in ./uniprot_sprot.fasta -dbtype prot
```

Then I ran blastp using the following bash script:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared                   # may consider running on a bigmem node for large dataset
#SBATCH -e blastp_%A.err            # File to which STDERR will be written
#SBATCH -o blastp_%A.out           # File to which STDOUT will be written
#SBATCH -J blastp_%A               # Job name
#SBATCH --mem=150000                 # Memory requested
#SBATCH --time=7-00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to


module load blast/2.6.0+-fasrc01

blastp -query /n/holyscratch01/hopkins_lab/Chaturvedi/transdecoder/Trinity.fasta.transdecoder.pep -db /n/holyscratch01/hopkins_lab/Chaturvedi/transdecoder/blastp/uniprot_sprot/uniprot_sprot.fasta -num_threads 10 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > td_uniprot.blastp.outfmt6

```
Output format:

Column headers:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

 1.	 qseqid	 query (e.g., unknown gene) sequence id
 2.	 sseqid	 subject (e.g., reference genome) sequence id
 3.	 pident	 percentage of identical matches
 4.	 length	 alignment length (sequence overlap)
 5.	 mismatch	 number of mismatches
 6.	 gapopen	 number of gap openings
 7.	 qstart	 start of alignment in query
 8.	 qend	 end of alignment in query
 9.	 sstart	 start of alignment in subject
 10.	 send	 end of alignment in subject
 11.	 evalue	 expect value
 12.	 bitscore	 bit score

3. Running Grouper (tutorial is here: https://github.com/COMBINE-lab/grouper/blob/master/README.md)

Before running grouper, we have to run SALMON. From the tutorial: Run Salmon on each sample in your experiment, passing it the --dumpEq option. This will tell Salmon to dump a representation of the fragment equivalence classes that it computed during quasi-mapping of each sample. If you wish to use orphan read information for joining contigs in Grouper, use the --writeOrphanLinks option as well, which will dump orphan read pair information to a file. Apart from these additional option, Salmon should be run normally (i.e. passing in whatever other options are appropriate for your samples).

SALMON tutorial:

First, build an index of the transcriptome.  The index is a structure that salmon uses to quasi-map RNA-seq reads during quantification. The index need only be constructed once per transcriptome, and it can then be reused to quantify many experiments. We use the index command of salmon to build our index:

```bash
salmon index -t /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trinity/Trinity_assembly/Trinity.fasta -i phlox_salmon_index
```
Now that we have our index built and all of our data downloaded, we’re ready to quantify our samples. Since we’ll be running the same command on each sample, the simplest way to automate this process is, again, a simple shell script (quant_tut_samples.sh):

```bash

for prefix in $(ls /n/holyscratch01/hopkins_lab/Chaturvedi/trim_galore/trimmed_reads/unfixrm_Henry_*001.cor_val_*.fq | sed -r 's/_R[12]_001[.]cor_val_[12].fq//' | uniq)

do

salmon quant -i phlox_salmon_index -l A -1 "${prefix}_R1_001.cor_val_1.fq" -2 "${prefix}_R2_001.cor_val_2.fq" -p 8 --validateMappings -o ${prefix}_quant

done

mv /n/holyscratch01/hopkins_lab/Chaturvedi/trim_galore/trimmed_reads/*quant /n/holyscratch01/hopkins_lab/Chaturvedi/grouper/salmon_quant/

```

On Jan 6th 2020
Adam suggested to rerun salmon with the following comments: So .. perhaps do two salmon runs, both with --dumpEq, and one with and the other without orphan reads. To integrate orphan reads with salmon you need to provide salmon the --writeOrphanLinks argument.

Here are the scripts for the new runs:

```bash
#!/bin/bash
#SBATCH -J salmon
#SBATCH -n 16                     # Use 1 cores for the job
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 6-00:00                 # Runtime in D-HH:MM
#SBATCH -p shared         # Partition to submit to
#SBATCH --mem=100000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o salmon_ind.%A.out  # File to which STDOUT will be written
#SBATCH -e salmon_ind.%A.err  # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to



module purge
module load salmon/0.12.0-fasrc01

for prefix in $(ls /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trimmed_filtered_seqs/unfixrm_Henry_*001.cor_val_*.fq | sed -r 's/_R[12]_001[.]cor_val_[12].fq//' | uniq)

do

salmon quant -i phlox_salmon_index -l A -1 "${prefix}_R1_001.cor_val_1.fq" -2 "${prefix}_R2_001.cor_val_2.fq" -p 8 --validateMappings --dumpEq -o ${prefix}_quant

done

mv ls /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trimmed_filtered_seqs/*quant /n/holyscratch01/hopkins_lab/Chaturvedi/salmon/salmon_quant_run2/
```

```bash
#!/bin/bash
#SBATCH -J salmon3
#SBATCH -n 16                     # Use 1 cores for the job
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 6-00:00                 # Runtime in D-HH:MM
#SBATCH -p shared         # Partition to submit to
#SBATCH --mem=100000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o salmon_ind3.%A.out  # File to which STDOUT will be written
#SBATCH -e salmon_ind3.%A.err  # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=schaturvedi@fas.harvard.edu # Email to send notifications to



module purge
module load salmon/0.12.0-fasrc01

for prefix in $(ls /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trimmed_filtered_seqs/unfixrm_Henry_*001.cor_val_*.fq | sed -r 's/_R[12]_001[.]cor_val_[12].fq//' | uniq)

do

salmon quant -i phlox_salmon_index -l A -1 "${prefix}_R1_001.cor_val_1.fq" -2 "${prefix}_R2_001.cor_val_2.fq" -p 8 --validateMappings --dumpEq --writeOrphanLinks -o ${prefix}_quant

done

mv ls /n/holystore01/LABS/hopkins_lab/Lab/schaturvedi/Henry_rnaseq/trimmed_filtered_seqs/*quant /n/holyscratch01/hopkins_lab/Chaturvedi/salmon/salmon_quant_run3/
```

SAM IS HERE
