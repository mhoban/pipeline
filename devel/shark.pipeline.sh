#!/bin/bash

#Notes/Questions
#cutadapt: why does the R2 read have the -m 30 otions set, but the R1 read doesn't?
#sample stats: How did you generate these? (esp "reads starting with cut site")

#Dependencies managed by homebrew:
#fastqc
#fastx_toolki

#Dependencies installed manually (unless otherwise noted, symlinked from $basedir/bin/ to $basedir/src/<project>)
#cutadapt (installed via python pip)
#PEAR
#rainbow
#Trim Galore


#General setup variables
########################################################################################
#data directories:
basedir=/Users/deadbilly/code/pipeline    #base directory for pipeline junk
datadir=$basedir/data/shark               #base directory for sequence data / config
configdir=$datadir/config                 #various config files live in here
fastqdir=$datadir/sample.fastq/           #directory where raw fastq files are stored
                                          #fastq filenames are expected to be in the 
                                          #format *_SAMPLEID_*_R1_*.fasta.gz/*_SAMPLEID_*_R2_*.fasta.gz
#working directories:
outputdir=$datadir/output
working_dir=$fastqdir


#executable directories:
sysutil=/usr/local/bin                    #utility binary directory (systemwide)
locutil=$basedir/bin                      #utility binary directory (locally installed)
########################################################################################
#Pipeline-specific config
#
#These variables will tell fastx toolkit to filter entire reads that
#do not have a quality score of at least 20 in 90% of bases
min_quality=20
quality_percent=90
#
# What files are we currently dealing with
working_files=*.fastq.gz
#Have our reads been paired?
paired=false
########################################################################################
#setup stuff:
mkdir -p $outputdir       #create working output directory
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $PIPEDIR/funx

function genQC() {
  local qc=$1
  local qcdir=$2
  if [[ ! -d $qcdir ]]; then
    mkdir -p $qcdir || return 1
  fi

  shift; shift

  if [[ $qc == true ]]; then
    # Call FastQC
    $sysutil/fastqc -o $qcdir $@
  fi

  # Generate HTML index
  echo "<html><body>" > $qcdir/index.html
  for html in $qcdir/*_fastqc.html; do
    echo "<a href=\"$html\">$(basename $html)</a><br>" >> $qcdir/index.html
  done
  echo "</html></body>" >> $qcdir/index.html

  # Prompt to open browser
  if ask "Launch browser to view QC reports? [y/N]"; then
    open $qcdir/index.html
    pause
  fi
}


############  Prepare REFERENCE GENOME #################
#import data from server
#wget --user=biocore --password=1234qwerty 'http://courge.ics.hawaii.edu/~mahdi/11_22_2013/Project_Jon_Whitney.tar.gz'


########### Preparing Raw Sequences ############

#Generate quality stats with fastqc
#####################################
if ask "Generate initial FastQC reports? [y/N]"
then
  genQC true $outputdir/QC $working_dir/$working_files
fi

#Trim Illumina adapters with Trim Galore
########################################
#Alternatively, pairing the reads with PEAR seems to do a lot of quality filtering for you
#We call a modified version of the trim glore script which uses the -b search algorithm
#(searches both 3' and 5' ends for adapter sequence)
if ask "Keep paired reads separate (trim galore, rainbow, etc.)? [Y/n]" 1
then
  echo "Using Trim Galore to trim Illumina adapters..."
  working_files=*_R1_*.fastq.gz

  # Forward and reverse TruSeq adapters
  fwd_adapt=GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
  rvs_adapt=GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
  for fq in $working_dir/$working_files
  do
    #Setup filenames
    fq_f=fq             #forward read file
    fq_r=${fq_f/R1/R2}  #reverse read file

    # pass to Trim Galore
    $locutil/trim_galore \
      --phred33 --paired \
      --length 70 \
      --retain_unpaired \
      --dont_gzip \
      -a $fwd_adapt \
      -a2 $rvs_adapt \
      --stringency 10 \
      -e 0.1 -q 10 -r1 75 -r2 75 \
      --output_dir $outputdir \
      --fastqc --fastqc_args="-o $outputdir/QC" \
      --clip_R1 4 --clip_R2 4 --three_prime_clip_R1 1 --three_prime_clip_R2 1 \
      $fq_f $fq_r
  done
  genQC false $outputdir/QC
  #Set working fileset to adapter-trimmed fastq files
  working_files=*_trim.fastq.gz
fi
#--------------------------------------------------
# ##################################################################
# #   Paired-end merging with PEAR (ignore this for the time being)
# ##################################################################
# elif ask "Attempt to merge paired reads with PEAR? [y/N]"
# then
#   echo "Merging paired-end reads..."
# 
#   mkdir -p $outputdir
#   working_files=*_R1_*.fastq.gz
#   for fq in $working_dir/$working_files
#   do
#     #PEAR options
#     #-n 33: minimum assembly length 33
#     #-t 33: minimum trip length 33
#     #-q 10: minimum quality score 10
#     #-j 2: set number of threads (processors)
#     #-u 0: disallow uncalled bases (throw out N's)
#     #-m 550: maximum length of reads
#     #-y 1g: use 1G of memory
#     $locutil/pear -f $fq \
#       -r ${fq/R1/R2} \
#       -o $outputdir/$(basename ${fq/_R1_/_merged_} .fastq.gz) \
#       -n 33 \
#       -t 33 \
#       -q 10 \
#       -j 2 \
#       -u 0 \
#       -m 550 \
#       -y 1g
#   done
#   #Set working fileset to assembled read files
#   working_files=*_merged_*.assembled.fastq
#   #Remember that our reads have been merged
#   paired=true
# fi
#-------------------------------------------------- 

#Switch working directory to location of proccessed files
working_dir=$outputdir

#--------------------------------------------------
#TODO: figure out if this is necessary
# echo "Quality-filtering reads:"
# echo "Discarding whole reads with <90% of bases having quality score >= 20..."
# pause
# # Filtering reads based on % quality (filters out poor quality reads overall)
# # -q 20 -p 90 means it filters whole reads that do not have a Qscore of 20 in 90% of bases
# # uses $min_quality and $quality_percent variables
# for fq in $working_dir/$working_files
# do
#   barefile=${fq##*/} #similar to $(basename fq)
#   barefile=${barefile%#.*} #does what $(basename fq .ext) doe for any extension
#   #Make sure we are transparent to whether the fastq file is gzipped or not
#   #if we used PEAR, they're not zipped, otherwise they probably are
#   #the product of this section will be gzipped fastq files with "_Qual20" tagged on the end
#   is_zipped $fq && catcmd=gzcat || catcmd=cat
#   $catcmd $fq | $sysutil/fastq_quality_filter \
#     -q $min_quality \
#     -p $quality_percent \
#     -z -v \
#     -o $(barefile)_Qual20.fastq.gz
# done
#-------------------------------------------------- 




#--------------------------------------------------
# ############# Filtering Unpaired Sequences #######################
# #TODO: I think we only need to filter unpaired reads if we've fastq_quality_filter'd first
#        otherwise, trim galore handles it for us
# #This is unneccessary if we've paired with PEAR
# if ! $paired
# then
# #If we have unpaired reads, we might have had to do this before calling cutadapt
# #usage: Perl filterNonPairedReads.pl <output_prefix> <read_1.fastq.gz> <read_2.fastq.gz>
# #requires reads to be compressed in .gz format
# #perl /home/jw2/scripts/filterNonPairedReads.pl MEL_filtered_cQTF MEL_R1_cleanQualTrimFilt.fastq.gz MEL_R2_cleanQualTrimFilt.fastq.gz
#   echo "Filtering out unpaired reads from forward and reverse fastq files..."
#   echo "Nothing actually happens here right now"
#   pause
# fi
#-------------------------------------------------- 

######################################################################################################################## 
################ RAINBOW: create de novo assemblies
######################################################################################################################## 
#make de novo assembly with rainbow 
#-m is number of mismatches allowed between two reads being compared for clustering (default = 4)
#Rainbow Assembly with MEL_filtered50 reads, allow 6 mismatches (-m 6) and select best+read1
/home/jw2/programs/rainbow_2.0.3/rainbow cluster -m 6 -1 MEL_filtered_cQTF.pair1.fastq  -2 MEL_filtered_cQTF.pair2.fastq > MEL_rbcluster.out 2> MEL_rblog
/home/jw2/programs/rainbow_2.0.3/rainbow div -i MEL_rbcluster.out -o MEL_rbdiv.out
/home/jw2/programs/rainbow_2.0.3/rainbow merge -o MEL_rbasm.out -a -i MEL_rbdiv.out
 
/home/jw2/programs/rainbow_2.0.3/select_best_rbcontig.pl MEL_rbasm.out > MEL_best_contig.fasta
/home/jw2/programs/rainbow_2.0.3/select_all_rbcontig.pl  MEL_rbasm.out > MEL_all_contig.fasta
/home/jw2/programs/rainbow_2.0.3/select_best_rbcontig_plus_read1.pl MEL_rbasm.out MEL_rbdiv.out > MEL_best+read1_contig.fasta


#Count # of contigs in each output
grep -c "^>" MEL_best_contig.fasta
grep -c "^>" MEL_best+read1_contig.fasta
grep -c "NNNNNNNNNN" MEL_best+read1_contig.fasta


#Best contigs:         922,782  (contigs where R1,R2 overlap)
#Best+read 1 contigs:  902,927  (contigs where R1,R2 DONT overlap, 2 contigs joined by 10 Ns)

#histogram of contig sizes (http://wiki.bioinformatics.ucdavis.edu/index.php/Count_fasta.pl)
perl /home/jw2/scripts/count_fasta.pl MEL_best_contig.fasta

# Total length of sequence: 146,353,824 bp
# Total number of sequences:  922,782

# 0:99            3
# 100:199   857,277
# 200:299    64,818
# 300:399       389
# 400:499       278
# 500:599        16
# 600:699         1

perl /home/jw2/scripts/count_fasta.pl MEL_best+read1_contig.fasta

# Total length of sequence: 254,119,690 bp
# Total number of sequences:  922,782

# 0:99            1
# 100:199     1,679
# 200:299   872,621
# 300:399    28,911
# 400:499    19,207
# 500:599       347
# 600:699        13
# 700:799         3


#Here is an example of the same contig in Best and Best+Read1 (E1_L296) 

# GATCGGAGCTGAGTCTAAATGCGGACGTGTAGGCAGAAGGATTGGAGCACTGTATGACCATTCAGAGCATTGTGACAGTGTGATGTTTTGCCTCGACTAGCCGTGAATTGTCATTCTGTGGTCATGTTAGCAATAAATACATATGCAATGACAAAAGTCCCAGTAACACAAGATAACGAGCCACATTTTTGAAGCTGAAGAAAACCATCTGCAAAGGATGGTTAGATATATCTTTTGTGTATGATGTAAAACTATATACAATAAAGAAATATATTGATTTGTTGCATTTACATGTT

# GATCGGAGCTGAGTCTAAATGCGGACGTGTAGGCAGAAGGATTGGAGCACTGTATGACCATTCAGAGCATTGTGACAGTGTGATGTTTTGCCTCGACTAGCCGTGAATTGTCATTCTGTGGTCATGTTAGCAATAAATACATATGCAATGACAAAAGTCCCAGTAACACAAGATAACGAGCCACATTTTTGAAGCTGAAGAAAACCATCTGCAAAGGATGGTTAGATATATCTTTTGTGTATGATGTAAAACTATATACAATAAAGAAATATATTGATTTGTTGCATTTACATGTTNNNNNNNNNNGTACAGCTTTCAACACAGCAGCAGTAAACAGCTTTTCTGAAAAAGTGTGTGGGCAGCAGAGCATCTTGCCACATACGTGTGAATGTGTCCTTGAGTGACTCTATGCTTGATC



###################################################################################################### 
# REMOVE mtDNA contigs from Rainbow assembly to create Nuclear contigs only

cd ~/ref/mito-blast

#Blast Parcatus assembly to Parcatus_Mtgenome
formatdb -i Parcatus_mt-genome.fasta -p F 
blastall -p blastn -i run3/MEL_best+read1_contig.fasta -d Parcatus_mt-genome.fasta -e 1e-3 -F F -o MEL_mito.out -b 5 -v 5 -a 6 -m 8

#use awk to filter output file
#e.g., extract alignments > 50 bases 
awk '$4>=50' MEL_mito.out > align50.out
#RESULT = 176 unique contigs align (>100 bases and 98-100% similarity) to Parcatus-mt-genome

#now extract/print queryID from that output
awk '{print $1}' align50.out > mt_contigs.ids

#extract IDfile from all contigs
grep -P "^>" run3/MEL_best+read1_contig.fasta | awk -F ">" '{print $2}'  > allcontigs.ids

#use COMM to print lines unique to ALL file  --> thus nuclear
comm -13  <(sort mt_contigs.ids | uniq) <(sort allcontigs.ids | uniq) > nuclear_contigs.ids
#RESULT = 932,801 IDs unique to allcontigs (i.e., nuclear) adds up 176 mito + 922,293 nuclear = 932,977 ALL

#create new Nuclear_contigs.fasta using nuclear_contigs.ids from MEL_best+read1_contigs
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' nuclear_contigs.ids run3/MEL_best+read1_contig.fasta > MEL_nuclear_contigs.fasta

#Double-check nuclear contigs for MTreads by blasting again
blastall -p blastn -i MEL_nuclear_contigs.fasta -d Parcatus_mt-genome.fasta -e 1e-3 -F F -o temp.out -b 5 -v 5 -a 6 -m 8
#RESULT: no alignments >50 bases = SUCCESSFULLY pulled out all mt-DNA contigs

perl /home/jw2/scripts/count_fasta.pl MEL_nuclear_contigs.fasta

# Total length of sequence: 227201953 bp
# Total number of sequences:  932801

# 0:99          139
# 100:199   219,668
# 200:299   582,955
# 300:399   126,860
# 400:499     2,932
# 500:599       198
# 600:699        33
# 700:799        14
# 800:899         1
# 900:999         0
# 1000:1099       1

#NOW append Parcatus-mt-genome.fasta to MEL_nuclear+contigs.fasta --> Parcatus_ref_cons95_May14.fasta
cp ~/rainbow/MEL_nuclear_contigs.fasta  ~/usearch
cd ~/usearch

########################################################################################################################################   Usearch ######################

#Usearch to cluster rainbow contigs and remove duplicate / reverse complement contigs

#1) dereplicate contigs (filter out duplicate contigs w/ 100% identity including reverse complement contigs)
~/programs/usearch -derep_fulllength MEL_nuclear_contigs.fasta -output uniques.fasta -sizeout -strand both -minseqlength 70

#Cluster contigs using usearch
#2) prep-step: sort by contig length from longest to shortest (to speed up algorithm and use less memory)
~/programs/usearch -sortbylength uniques.fasta -output uniques_sorted.fasta -minseqlength 70

#3) cluster sequences using cluster_smallmem (only one that allows searching for reverse complements)
~/programs/usearch -cluster_smallmem uniques_sorted.fasta -id 0.95 -centroids centroids.fasta -strand both -uc clusters.uc -consout consensus95.fasta -maxseqlength 1000 -dbstep 8
    
    #NOTE: first run failed because limited to only 4.4gb memory. Ran subsequent steps with -dbstep 8 (which matches the word length so it covers all bases with non-overlapping searches). default word length is 8. Increasing -wordlength 16 -dbstep 16 should reduce memory (if needed). The -maxseqlength 10000 is meant to exclude the Mt-genome contig (17k) to reduce memory.

    #rerun using shared usearch and larger non-overlapping words of 16
    #/programs/usearch -cluster_smallmem uniques_sorted.fasta -id 0.95 -centroids centroids.fasta -strand both -uc clusters2.uc -consout consensus95-2.fasta -maxseqlength 1000 -wordlength 16 -dbstep 16

#RESULTS
# 932,802   Initial sequences (from rainbow)  min 71, avg 244, max 17336nt
# -24,659   100% identical sequences (2.7%)
# 908,143   Unique sequences (97.3%), min 71, avg 245, max 17336nt
#-240,500   Redundant contigs
# 667,643   Clusters (73.5%) based on 95% similarity

#usearch will output both centroids and consensus .fasta files. Consensus is a consensus of all aligned contigs, whereas centroids is the contig used as the centroid of the cluster. The usearch manual states that consensus contigs will probably reduce error (taking majority of multiple alignments) as opposed to centroids, which represent a single contig. 

# #consensus95           input (rainbow out) 
#   0:099         648        139
# 100:199     123,717    219,668
# 200:299     426,749    582,955
# 300:399     113,779    126,860
# 400:499       2,583      2,932
# 500:599         142        198
# 600:699          19         33
# 700:799           6         14
# 800:899           0          1
# 900:999           0          0
# 1000:1099         0          1
# 17300(mt-geno)    0          1

# #centroids95
#   0:099          46
# 100:199     122,210
# 200:299     425,976
# 300:399     116,325
# 400:499       2,855
# 500:599         183
# 600:699          32
# 700:799          14
# 800:899           1
# 900:999           0
# 1000:1099         1


#re-run usearch just changing -id (similarity) to see effect of varying threshold
#~/programs/usearch -cluster_smallmem uniques_sorted.fasta -id 0.98 -centroids centroids98.fasta -strand both -uc clusters98.uc -consout consensus98.fasta -maxseqlength 10000 -dbstep 8

#RESULTS
# 932,802   Initial sequences (from rainbow)  min 71, avg 244, max 17336nt
# 908,143   Unique sequences (97.3%), min 71, avg 245, max 17336nt
# 741,295   98% threshold   (18% contigs 98-99.9% similar)
# 708,661   97%     "       (3.6% contigs 97-98% similar)
# 688,685   96%     "       (2.2% contigs 96-97% similar)
# 667,643   95%     "       (2.3% contigs 95-96% similar)

#NOTE - there was an issue with the consensus files - there were some sequences in the fasta that only had a few bases and others that were very truncated. I think that it takes a consensus only of the overlap of all reads (so the size of the consensus is the length of the smallest aligned sequence - is this right?) Regardless, I think the centroid sequences are better - as it pulls the centroid that is most dissimilar from other centroids (and they happen to be longer)


##########################################
#Fixing Headers from consensus and centroid output
#fasta headers come out all funky from usearch, like this:
#>centroid=E541063_L749;size=1;;seqs=10;

#use sed to remove the "centroid=" prefix | then remove everything after the contigID
sed 's/>centroid=/>/' consensus95.fasta | sed 's/\;.*$//' > cons95_clean.fasta
#NOTE:(semicolon is a special character in sed so it needs to be preceded with \ backslash to tell it too look for literal ";"). So in words, we're saying look for pattern s/  /  where ; and any character (.) one or more times (*) will be substituted ($) with (nothing) 
#try to pipe it

#clean up CENTROID.fasta
#centroid headers are like this:
#>E388599_L789;size=1;

#so just need to remove everything after the contigID
sed 's/\;.*$//' centroids95.fasta > cen95_clean.fasta


#remove partial mt-DNA sequence
#1) extract IDs from file
grep -o -E "^>\w+" cen95_clean.fasta | tr -d ">" > cen95.ids
#2) remove Mt-DNA ID
grep -v "^gi" cen95.ids > cen95.ids2
#3) use perl to extract IDs from fasta
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' cen95.ids2 cen95_clean.fasta > cen95_clean2.fasta

mv cen95_clean2.fasta centroid95_clean.fasta
#####################################

#export centroids95.fasta to reference folder
cp ~/usearch/centroid95_clean.fasta ~/ref

#blast consensus contigs with Mt-genome to double check for any mt-contigs
cd ~/ref/mito-blast 
blastall -p blastn -i ~/usearch/centroid95_clean.fasta -d ~/ref/mito-blast/Parcatus_mt-genome.fasta -e 1e-3 -F F -o Mt-check.out -b 5 -v 5 -a 6 -m 8
awk '$4>=50' Mt-check.out | less
#results not a single match above alignment50 --> no mt-contigs in reference

cd ~/ref


#concatenate consensus contigs with Mt-genome to create reference
cat centroid95_clean.fasta Parcatus_mt-genome.fasta > Parcatus_ref_cent95_May14.fasta

#check that it contains 667,643 contigs and one big mt-genome at end
perl /home/jw2/scripts/count_fasta.pl Parcatus_ref_cent95_May14.fasta
#yes it does --> good to proceed.

#also check that all headers in fasta are clean (without extra crap)
#GOOD TO GO




#######################################################################################

######################################## 
################ Alignments of sample to reference 
########################################

#first prep contigs from de novo assembly (and clustering) for alignment by indexing it (in other words your turning your contigs into your reference genome.
#generate BWA index 
bwa index -a is Parcatus_ref_cent95_May14.fasta
#generate the fasta file index 
samtools faidx Parcatus_ref_cent95_May14.fasta
#generate the sequence dictionary
java -jar /home/jw2/programs/picard-tools-1.108/CreateSequenceDictionary.jar R=Parcatus_ref_cent95_May14.fasta O=Parcatus_ref_cent95_May14.dict

cd ~/alignments

#Use BWA to align each read set (read 1 and 2 separate) to the same reference genome
~/programs/bwa-0.7.8/bwa aln ../ref/Parcatus_ref_cent95_May14.fasta   -t 16 MEL_trimmed_paired_1.fastq  -f MEL_1.sai
~/programs/bwa-0.7.8/bwa aln ../ref/Parcatus_ref_cent95_May14.fasta   -t 16 MEL_trimmed_paired_2.fastq  -f MEL_2.sai
~/programs/bwa-0.7.8/bwa aln ../ref/Parcatus_ref_cent95_May14.fasta   -t 16 PWS_trimmed_paired_1.fastq  -f PWS_1.sai
~/programs/bwa-0.7.8/bwa aln ../ref/Parcatus_ref_cent95_May14.fasta   -t 16 PWS_trimmed_paired_2.fastq  -f PWS_2.sai

#use sampe to align the paired ends together into a single .sam file
~/programs/bwa-0.7.8/bwa sampe ../ref/Parcatus_ref_cent95_May14.fasta  MEL_1.sai   MEL_2.sai   MEL_trimmed_paired_1.fastq   MEL_trimmed_paired_2.fastq -f MEL.sam

~/programs/bwa-0.7.8/bwa sampe ../ref/Parcatus_ref_cent95_May14.fasta   PWS_1.sai   PWS_2.sai   PWS_trimmed_paired_1.fastq   PWS_trimmed_paired_2.fastq  -f PWS.sam

#then samtools converts .sam into .bam (binary file thats highly compressed)
samtools import ../ref/Parcatus_ref_cent95_May14.fasta   MEL.sam   MEL.bam
samtools import ../ref/Parcatus_ref_cent95_May14.fasta   PWS.sam   PWS.bam

#use samtools to sort .bam files by the contig ID in the reference file
samtools sort MEL.bam  MEL_sorted
samtools sort PWS.bam  PWS_sorted

#MAPPING STATS - MEL
samtools flagstat MEL_sorted.bam > MEL_sorted.bam.flagstat

# 31937030 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 23991625 + 0 mapped (75.12%:-nan%)
# 31937030 + 0 paired in sequencing
# 15968515 + 0 read1
# 15968515 + 0 read2
# 22880516 + 0 properly paired (71.64%:-nan%)
# 23629058 + 0 with itself and mate mapped
# 362567 + 0 singletons (1.14%:-nan%)
# 689600 + 0 with mate mapped to a different chr
# 567178 + 0 with mate mapped to a different chr (mapQ>=5)



#EXTRACT ONLY PROPERLY PAIRED READS

#Extract properly paired (-f 0x02) AND exclude un-mapped reads (-F 0x04) (PM = paired mapped)
samtools view -f 0x02 -F 0x04 -b MEL_sorted.bam > MEL_sorted_PM.bam

samtools flagstat MEL_sorted_PM.bam > MEL_sorted_PM.bam.flagstat

# 22827292 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 22827292 + 0 mapped (100.00%:-nan%)
# 22827292 + 0 paired in sequencing
# 11407568 + 0 read1
# 11419724 + 0 read2
# 22827292 + 0 properly paired (100.00%:-nan%)
# 22827292 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00%:-nan%)
# 19696 + 0 with mate mapped to a different chr
# 19237 + 0 with mate mapped to a different chr (mapQ>=5)

#########################################

#MAPPING STATS - PWS
samtools flagstat PWS_sorted.bam > PWS_sorted.bam.flagstat

# 30510322 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 22533656 + 0 mapped (73.86%:-nan%)
# 30510322 + 0 paired in sequencing
# 15255161 + 0 read1
# 15255161 + 0 read2
# 21426014 + 0 properly paired (70.23%:-nan%)
# 22123025 + 0 with itself and mate mapped
# 410631 + 0 singletons (1.35%:-nan%)
# 665325 + 0 with mate mapped to a different chr
# 542813 + 0 with mate mapped to a different chr (mapQ>=5)

###############################

#Extract properly paired (-f 0x02) AND exclude un-mapped reads (-F 0x04) (PM = paired mapped)
samtools view -f 0x02 -F 0x04 -b PWS_sorted.bam > PWS_sorted_PM.bam

samtools flagstat PWS_sorted_PM.bam > PWS_sorted_PM.bam.flagstat

# 21373995 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 21373995 + 0 mapped (100.00%:-nan%)
# 21373995 + 0 paired in sequencing
# 10688257 + 0 read1
# 10685738 + 0 read2
# 21373995 + 0 properly paired (100.00%:-nan%)
# 21373995 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00%:-nan%)
# 21225 + 0 with mate mapped to a different chr
# 20781 + 0 with mate mapped to a different chr (mapQ>=5)

################################

# use picard tools to add read groups to each .bam file (gives an individual or pool)
java -XX:MaxPermSize=64g -jar /home/jw2/programs/picard-tools-1.108/AddOrReplaceReadGroups.jar \
   I=MEL_sorted_PM.bam  O=MEL_sorted_PM_WRG.bam  RGLB=MEL RGPL=illumina RGPU=MEL  RGSM=MEL_pool30 RGCN=GenomicsCore  VALIDATION_STRINGENCY=SILENT 
 
java -XX:MaxPermSize=64g -jar /home/jw2/programs/picard-tools-1.108/AddOrReplaceReadGroups.jar  \
   I=PWS_sorted_PM.bam O=PWS_sorted_PM_WRG.bam  RGLB=PWS RGPL=illumina RGPU=PWS RGSM=PWS_pool30 RGCN=GenomicsCore VALIDATION_STRINGENCY=SILENT 

#index bams with ReadGroup    
#samtools index MEL_sorted_PM_WRG.bam
#samtools index PWS_sorted_PM_WRG.bam

#Use MergeSamFiles (Picard-tools) to merge BAM files
java -XX:MaxPermSize=64g -jar /home/jw2/programs/picard-tools-1.108/MergeSamFiles.jar \
   INPUT=MEL_sorted_PM_WRG.bam INPUT=PWS_sorted_PM_WRG.bam OUTPUT=BOTH_sorted_PM_WRG.bam

#Validation that both RGs are in merged file
~/programs/freebayes/bamtools/bin/bamtools header -in BOTH_sorted_PM_WRG.bam | grep -P "@RG" > RG_validate.log 

#@RG     ID:1    CN:GenomicsCore LB:MEL  PL:illumina     PU:MEL  SM:MEL_pool30
#@RG     ID:1.1  CN:GenomicsCore LB:PWS  PL:illumina     PU:PWS  SM:PWS_pool30

#index merged bam file
samtools index BOTH_sorted_PM_WRG.bam

cp BOTH_sorted_PM_WRG.bam ~/freebayes
cp BOTH_sorted_PM_WRG.bam.bai ~/freebayes
cd ~/freebayes 

#########################################################################################
########### Variant Calling with FreeBayes ##############################################

#options
# -C 2   #variable in at least 2 reads  (from a single sample)
# -F .2  #variable in at least 20% of reads (from a single sample)
# -E    #--max-complex-gap   #usually set to half of read length (default: 3)
# -p    #ploidy of samples (default = 2)
# -J    #pooled samples. discrete model
# -K    #pooled samples. continous model
# -l   #use all variant calls in VCF file to produce genotype likelihoods
# -m    #minimum mapping quality (default 0)
# -q    #exclude alleles if their base quality is < Q (default: 0)
# -j    #use mapping quality of alleles when calculating likelihoods
# --genotype-qualities    #report GQ in VCF
# --populations Pop.file  #designate a pop for each sample file
# --failed-alleles Failed-Alleles.out 


#FreeBayes
#we can use 3 variant calling methods for pooled samples: J (discrete), K (continuous), and Naive.  In short, -K (pooled-continous) is the best method that generates frequency-based calls for all variants. J-mode calculates genotypes from pooled data and looks for variants between pools (this is over-simplified and produces very few SNPS - highly underestimated), and Naive takes way too long and calls a SNP for each site and counts #s of each variant (it will likely overestimate SNPs)


#K-mode Continuous-pooling (no-ploidy set and ignores genotype outcome). Generates frequency-based calls for all variants passing input thresholds
~/programs/bin/freebayes -f ../ref/Parcatus_ref_cent95_May14.fasta  --bam BOTH_sorted_PM_WRG.bam \
    -m 1 -q 3 -E 3 -j -F 0.05 -C 2 --pooled-continuous --min-coverage 10  --vcf BOTH_poolK_cent95_May7.vcf 

#J-mode Discrete-Pooling (modeled using discrete genotypes across pools)
#~/programs/bin/freebayes -f ../ref/Parcatus_ref_cent95_May14.fasta --bam BOTH_sorted_PM_WRG.bam \
    #-m 1 -q 3 -E 3 -p 60 -j --pooled-discrete --min-coverage 10 --vcf BOTH_poolJ_rawVar.vcf 

#naive variant calling: simply annotate observation counts of SNPs and indels
 #~/programs/bin/freebayes -f ../ref/Parcatus_ref_cent95_May14.fasta  --bam BOTH_sorted_PM_WRG.bam \
  #  --haplotype-length 0  --min-alternate-count 1  --min-alternate-fraction 0  --pooled-continuous \
 #   --report-monomorphic  --haplotype-length 0  --vcf BOTH_naive_rawVar.vcf




############### VcfLib Tools to work with VCF files ##################

#Filter vcf results
~/programs/vcflib/bin/vcffilter -f "QUAL > 20" BOTH_poolK_cent95_May7.vcf > BOTH_poolK_Fq20.vcf 

### PoolK vcf files are too big for excel - so cut down to SNPS 
~/programs/vcflib/bin/vcffilter -f "TYPE = snp" BOTH_poolK_Fq20.vcf > BOTH_poolK_Fq20_snps.vcf 

#Filter to include only SNPs with total Depth > 20
~/programs/vcflib/bin/vcffilter -f "DP > 20" BOTH_poolK_Fq20_snps.vcf  > BOTH_poolK_Fq20_snps_20x.vcf 

#Annotate variants with the distance to the nearest variant
~/programs/vcflib/bin/vcfdistance < BOTH_poolK_Fq20_snps_20x.vcf   > BOTH_poolK_Fq20_snps_20x_v2.vcf  

#Export vcf file as tab-delimited file
~/programs/vcflib/bin/vcf2tsv BOTH_poolK_Fq20_snps_20x_v2.vcf  > BOTH_poolK_Fq20_snps_20x_v2.tsv

#gzip output files to scp
gzip -c -1 -v BOTH_poolK_Fq20_snps_20x_v2.vcf > BOTH_poolK_Fq20_snps_20x_v2.vcf.gz
gzip -c -1 -v BOTH_poolK_Fq20_snps_20x_v2.tsv > BOTH_poolK_Fq20_snps_20x_v2.tsv.gz


#-----------------------------------------------------------------------------

cd ~/alignments 


############## SNP and INDEL Realigner with GATK   ################

#Determine Depth of Coverage 
#java -Xmx16g -jar /home/jw2/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
#   -R ../ref/Parcatus_ref_cent95_May14.fasta \
#   -T Coverage \
#   -o DepthCov_MEL_filt50_cd95_ref \
#   -I MEL_sorted_WRG.bam

#Realign Targets - MEL + PWS together
java -XX:MaxPermSize=64g -jar /home/jw2/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -R ../ref/Parcatus_ref_cent95_May14.fasta \
 -T RealignerTargetCreator  \
 -I BOTH_sorted_PM_WRG.bam -o both_realigner.intervals -nt 16

#Indel Realigner - MEL + PWS together
java -XX:MaxPermSize=64g -jar /home/jw2/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T IndelRealigner  -R ../ref/Parcatus_ref_cent95_May14.fasta \
     -I BOTH_sorted_PM_WRG.bam  \
     -targetIntervals both_realigner.intervals \
     -o BOTH_sorted_PM_WRG_realigned.bam 



########### SNP calling ########################################################################################################################################

#need to manually edit freebayes output and convert to .sync file

#import .sync file from local (created in excel)

cd ~/freebayes 
mkdir run5_freebayes
cp BOTH_poolK_SimpleSNPS_syncmake_EXPORT_2.prn run5_freebayes
cd run5_freebayes
mv BOTH_poolK_SimpleSNPS_syncmake_EXPORT_2.prn BOTH_poolK_snps_sync_raw

#When importing file from excel need to remove ^M carriage returns
sed 's/\cM/\n/g' BOTH_poolK_snps_sync_raw > BOTH_poolK_snps.sync

#1) #Calculate allele frequency differences
perl /home/jw2/programs/popoolation2/snp-frequency-diff.pl --input BOTH_poolK_snps.sync --output-prefix BOTH_poolK_snps --min-count 6 --min-coverage 20 --max-coverage 1000

#This script creates two output files having two different extensions:
#_rc: this file contains the major and minor alleles for every SNP in a concise format
#_pwc: this file contains the differences in allele frequencies for every pairwise comparision of the populations present in the synchronized file


######## Fst-values: measure differentiation between populations ###########

#Calculate Fst for every SNP
perl /home/jw2/programs/popoolation2/fst-sliding.pl --input BOTH_poolK_snps.sync --output BOTH_poolK_snps.fst --suppress-noninformative --min-count 6 --min-coverage 20 --max-coverage 1000 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 30
#the option --suppress-noninformative suppresses output for windows not containing a SNP. When applied to windows of size one, this option suppresses output for bases that are no SNP.

######## Fishers Exact Test #################################################

#Fisher's Exact Test: estimate the significance of allele frequency differences
#The Fishers exact test can be used to test whether any differences in allele frequencies
#are statistically significant. At low coverages the absolute changes of allele frequencies 
#or the Fst values may be strongly influenced by sampling effects, therefore the 
#Fishers exact test may be used to identify significant changes in allele frequency.
#NOTE: Fisher-test.pl requires Text-NSP to be installed (see below for install instructions)

#Run Fisher's test to calculate p-values for allele frequency differences
perl /home/jw2/programs/popoolation2/fisher-test.pl --input BOTH_poolK_snps.sync --output BOTH_poolK_snps_ALL.fet --min-count 6 --min-coverage 20 --max-coverage 1000

######## Combine Fst w/ FET outputs ########################################## 

#Both outputs have same file format and can be combined.

#1. use awk to filter FET to only include SNP sites
awk '$4>0' BOTH_poolK_snps_ALL.fet | awk '$3>0' > BOTH_poolK_snps.fet

#2. use awk to bind column 6 (-log10p-value) from fet to all columns of fst when the first 2 columns (CHROM and POS) match both fst and fet
awk 'NR==FNR{a[$1,$2]=$6;next} ($1,$2) in a{print $0, a[$1,$2]}' BOTH_poolK_snps.fet BOTH_poolK_snps.fst > BOTH_poolK_snps_join.txt 

#3. use sed to remove the "1:2=" prefix before the fst and -logP columns
#NOTE: fst is column 6 and -log10P is column 7
sed 's/1:2=//g' BOTH_poolK_snps_join.txt > BOTH_poolK_snps_join_clean.txt




#########################################################################
####### BLAST contigs that contain SNPs of interest######################

#extract contig ids that meet whatever threshold
#e.g., extract SNP ids that have Fst > .1 and are significant p < 0.001 (-log10P > 3)
awk '$6>0.1' BOTH_poolK_snps_join_clean.txt | awk '$7>3' | awk '{print $1}' | sort | uniq > sigContigs.ids 

#use perl script to extract contigs from fasta file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' sigContigs.ids ../reference/Parcatus_ref_cent95_May14.fasta > sigContigs.fasta




######################################################################################################
#BELOW HERE LIES JUNK (some useful notes)
######################################################################################################
######################################################################################################
######Commented-out stuff from Jon's pipeline moved down here###########
######################################################################################################
#--------------------------------------------------
# #post-adapter cleaning concatenating
#TODO: why this?
# #concatenate R2 untrimmed reads + trimmed reads
# cat MEL_R2_untrim.fastq MEL_R2_trimL30.fastq > MEL_R2_clean.fastq 
# cat PWS_R2_untrim.fastq PWS_R2_trimL30.fastq > PWS_R2_clean.fastq 
#-------------------------------------------------- 

#unzip tar.gz
#tar -zxvf Project_Jon_Whitney
#gunzip *.fastq.gz   #unzip all sequences with .gz

#concatenate files (if each read set was parsed into multiple files, e.g., if more than 4 million reads on GAIIx)
#cat MEL_CTTGTA_L005_R1_00*.fastq > MEL_R1.fastq
#cat MEL_CTTGTA_L005_R2_00*.fastq > MEL_R2.fastq
#cat PW5_GCCAAT_L005_R1_00*.fastq > PWS_R1.fastq
#cat PW5_GCCAAT_L005_R2_00*.fastq > PWS_R2.fastq
#/home/mahdi/programs/FastQC/fastqc  MEL_R1.fastq
#/home/mahdi/programs/FastQC/fastqc  PWS_R1.fastq
#/home/mahdi/programs/FastQC/fastqc  MEL_R2.fastq
#/home/mahdi/programs/FastQC/fastqc  PWS_R2.fastq
#cutadapt stuff:
#Read 2 - MEL
#--------------------------------------------------
# cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
# -e 0.1 -O 10 -m 30 --untrimmed-output=MEL_R2_untrim.fastq MEL_R2.fastq \
# -o MEL_R2_trimL30.fastq
# #Read 1 - PWS
# cutadapt -g ^AGATCGGAAGAGCACACGTCTGAACTCCAGTCACgccaatATCTCGTATGCCGTCTTCTGCTTG \
# -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACgccaatATCTCGTATGCCGTCTTCTGCTTG \
# -e 0.10 -O 10 --untrimmed-output=PWS_R1_untrim.fastq PWS_R1.fastq -o PWS_R1_trim.fastq
# #Read 2 - PWS
# cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
# -e 0.1 -O 10 -m 30 --untrimmed-output=PWS_R2_untrim.fastq PWS_R2.fastq \
# -o PWS_R2_trimL30.fastq
#-------------------------------------------------- 

#The way I was initially tryin to use cutadapt
  #--------------------------------------------------
  #i5=AATGATACGGCGACCACCGAGATCTACAC_BARCODE_ACACTCTTTCCCTACACGACGCTCTTCCGATCT
  #i7=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC_BARCODE_TCGTATGCCGTCTTCTGCTTG
  # #This will (primitively)read in barcodes from config file
  # #config file must be in proper format
  # #sample forward reverse
  # #lines can be commented with " (double quote)
  # 
  # #B0rk if cannot read barcodes file
  # #TODO: user can enter barcodes file location
  # if [ ! -f $barcodes ]
  # then
  #   echo "Barcodes file not found"
  #   exit 1
  # fi
  # 
  # #Read barcode file, line by line
  # while read sample fwd rvs
  # do
  #   #putting the ^# pattern right in code effs up syntax
  #   #highlighting in vi, so we store it in this variable
  #   comment="^#" 
  #   if [[ ! $sample =~ $comment ]]
  #   then
  #     #get forward and reverse adapters, with appropriate indices
  #     #we leave the indices lowercase, for fun and legibility
  #     forward=${i5/_BARCODE_/$(echo $fwd | tr A-Z a-z)}
  #     reverse=${i7/_BARCODE_/$(echo $rvs | tr A-Z a-z)}
  # 
  #     #setup filenames
  #     fwd_fastq=$(echo $fastqdir/*_${sample}_*_R1_*.fastq.gz)
  #     fwd_fastq_trim=$(basename $fwd_fastq .fastq.gz)_trim.fastq.gz
  #     fwd_fastq_untrim=$(basename $fwd_fastq .fastq.gz)_untrim.fastq.gz
  #     rvs_fastq=$(echo $fastqdir/*_${sample}_*_R2_*.fastq.gz)
  #     rvs_fastq_trim=$(basename $rvs_fastq .fastq.gz)_trim.fastq.gz
  #     rvs_fastq_untrim=$(basename $rvs_fastq .fastq.gz)_untrim.fastq.gz
  #     
  #     #run cutadapt
  #     #options:
  #     #-O 10: minimum overlap length
  #     #-a forward adapter
  #     #-A reverse adapter
  #     #-o, -p trimmed forward/reverse output file
  #     #--untrimmed(-paired?)-output untrimmed output files
  #     $sysutil/cutadapt \
  #       -O 10 \
  #       -a $forward \
  #       -A $reverse \
  #       -o $outputdir/$fwd_fastq_trim \
  #       -p $outputdir/$rvs_fastq_trim \
  #       --untrimmed-output=$outputdir/$fwd_fastq_untrim \
  #       --untrimmed-paired-output=$outputdir/$rvs_fastq_untrim \
  #       $fwd_fastq $rvs_fastq   #input fastq files
  #   fi
  # done < $barcodes
  #-------------------------------------------------- 
#Sample Stats on Hawkfish Reads
###### Initial number of sequences in ALL_R1:            34,195,059
###### Num seqs after filtering in    ALL_R1_QualFilt20: 28,237,739 (82.5%)
###### Initial number of sequences in All_R2:            34,195,059
###### Num seqs after filtering in    ALL_R2_QualFilt20: 16,481,952 (48.2%) 
#NOTE: last 50bases of read 2 are low quality, so first need to quality trim then quality filter
 
###### Breakdown by individuals:
### MEL R1 Before:    17,432,597
###        After-q20: 17,415,442
### MEL R2 Before:    17,432,597
###        After-q20:  8,590,968
### PWS R1 Before:    16,762,462
###        After-q20: 13,770,990
### PWS R2 Before:    16,762,462
###        After-q20:  7,890,984

####### Reads starting with cut site (avg ~ 75% start with cut site)
### MEL R1 Before:    17,432,597
### start w/ GATC:    12,582,727 (72%)
### MEL R2 Before:    17,432,597
### start w/ GATC:    13,678,349 (78%)
### PWS R1 Before:    16,762,462
### start w/ GATC:    12,464,395 (74%)
### PWS R2 Before:    16,762,462
### start w/ GATC:    13,406,845 (80%)
#--------------------------------------------------
# ## Alternative to filtering Orphaned Reads  using FLEXBAR
# flexbar -r read1.fastq -p read2.fastq -f i1.8
# 
# perl flexbar -r MEL_R1_untrim.fastq -p MEL_R2_clean.fastq -f i1.8 -n 12 -q 20 -m 30 -u 3 -j -s
# 
# -n   # threads to employ
# -r   #read1
# -p   #read2
# -f   #format (sanger, illumina, etc)
# -q   #trim 3'end until specific q score is reached
# -m   #min read length to remain after trimming
# -u   #max number of Ns allowed (default = 0)
# -j   #generate histogram of read lengths
# -s   #write orphaned reads to separate file
#-------------------------------------------------- 



#--------------------------------------------------
# ### Alternative to filtering using TrimGalore
# #written in perl so all commands must be preceded by perl
# perl ~/programs/bin/trim_galore  -v
# 
# --phred33  #sanger/Illumina 1.9+ (default)
# -q        #quality filter (default = 20)
# -s        #stringency for overlap (-O setting in cutadapt)
# -e        #error rate (default: 0.1)
# --paired  #maintains paired-reads
# --length  #keep sequences above set length
# --retain_unpaired  #keeps orphaned reads
# --dont_gzip
# -r1     #length cutoff to retain read1
# -r2     #length cutoff to retain read2
# 
# 
# perl ~/programs/bin/trim_galore --phred33 --paired --length 70 --retain_unpaired --dont_gzip \
#     -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG \
#     -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT   \
#     --stringency 10 -e 0.1 -q 10 -r1 75 -r2 75 \
#     MEL_1.fq MEL_2.fq
#-------------------------------------------------- 

#OUTPUT
#Results are not as stringent as doing it doing it in cutadapt 


#Run again for just mate-pairing function

#--------------------------------------------------
# perl ~/programs/bin/trim_galore --phred33 --paired --length 30 --retain_unpaired --dont_gzip \
#     -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG \
#     -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT   \
#     --stringency 58 -e 0.1 -r1 75 -r2 50 -q 10 \
#     MEL_cQTF_1.fq MEL_cQTF_2.fq
#-------------------------------------------------- 

#OUTPUT
#Input R1: 
#Input R2: 
#paired:   
#TODO: this perl one-liner gives sorted length distribution of fastq files
#Get length distribution of fastq
#cat <input_fastq> | perl -ne '$s=<>;<>;<>;chomp($s);print length($s)."\n";' | sort | uniq -c
###### Alternative Filtering Method with Read 2s
#wondering if extra long reads on read 2s are reason why 2x more read2s are filtered
#testing to see if we truncate all R2 at 110 if we get around this issue
#/home/jw2/programs/fastx/fastx_trimmer -l 110 -v -i MEL_R2_clean.fastq -o MEL_R2_cleanTrim110.fastq -Q 33
#/home/jw2/programs/fastx/fastq_quality_filter -q 20 -p 90 -v -i MEL_R2_cleanTrim110.fastq -o MEL_R2_cleanTrim110QualFilt.fastq -Q 33  

### MEL R2 Before:      17,431,143
### Qualtrim + q20,p90: 15,761,208  (-10% total)
### Trim110 + q20,p90:  14,326,930  (-17% total)

#Process R1 differently to keep intact (same read length - for rainbow assembly)
#/home/jw2/programs/fastx/fastq_quality_filter -q 20 -p 90 -v -i  MEL_R1_untrim.fastq -o MEL_R1_cleanQualFilt.fastq -Q 33

### MEL R1 Before:          17,423,468
### Qualtrim + q20,p90:     16,325,778  (-6% total)
### q20,p90 w/out qualtrim: 14,464,798  (-16% total)

#### RESULT: --> filtering by quality q20,p90 without trimming poor quality 3' ends results in an additional
#loss of 10% of reads (because we are throwing out reads with lower quality 3' ends).  Better to trim them


### Alternative - try Mahdi's Combine PE 
#/home/jw2/scripts/fastqCombinePairedEnd.py

#USAGE: python fastqCombinePairedEnd.py input1 input2 separator

# input1 = LEFT  fastq file (R1)
# input2 = RIGHT fastq file (R2)
# separator = character that separates the name of the read from the part that
#     describes if it goes on the left or right, usually with characters '1' or
#     '2'.  The separator is often a space, but could be another character. A
#     space is used by default.

/home/jw2/scripts/fastqCombinePairedEnd.py MEL-cQTF_1.fastq MEL-cQTF_2.fastq _
/home/jw2/scripts/fastqCombinePairedEnd.py PWS-cQTF_1.fastq PWS-cQTF_2.fastq _
