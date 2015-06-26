#!/bin/bash

#--------------------------------------------------
# chrom=$1
# pos=$2
# rg=$3
#-------------------------------------------------- 


bases=$(mktemp --tmpdir basesXXX)
echo -e "A\nC\nG\nT" > "$bases"



baseCount() {
  GREP_OPTIONS=
  ## this perl code interprets the SAM file CIGAR string and outputs the correct sequence
perlcode=$(cat <<'PERL'
chomp;
$p = 0;
($pos,$start,$cigar,$seq)=split /\t/;
next if ($start > $pos);
$thispos = $pos-$start+1;
while ($cigar =~ m/(\d+)([MIDNSHP=X])/g) {
  last if ($p >= $thispos);
  if ($2 eq "M" || $2 eq "=" || $2 eq "X") { $p += $1; }
  elsif ($2 eq 'D') { $thispos -= $1; }
}
$p >= $thispos ? print substr($seq,$thispos-1,1),"\n" : print "\n";
PERL
)
  ## specific to my whole deal
  ## it's an associative array mapping SAM read groups to population name
  declare -A samplemap
  samplemap["GROUP-1"]="red sea"
  samplemap["GROUP-3"]="indian"
  samplemap["GROUP-2"]="se pacific"
  samplemap["GROUP-5"]="ne pacific"
  samplemap["GROUP-4"]="n central pacific"
  samplemap["GROUP-6"]="w central pacific"
  samplemap["GROUP-7"]="atlantic"
  samplemap["GROUP-8"]="taiwan"
  ## Only use SNPs
  if [[ ${#5} == 1 ]]
  then
    echo -en "${samplemap[$1]}\t$2\t$3\t$5\t"
    ## samtools spits out the alignment for this SNP
    samtools view -r $1 alignments.merged.bam $2:$3 | 
      awk 'BEGIN {FS="\t"; OFS="\t";} {if ($4<'$3') print "'$3'",$4,$6,$10 }' | 
      perl -n -e "$perlcode" | 
      grep -h "^[ACGT]$" - "$4" | 
      sort | uniq -c | awk '{print $1-1}' | tr "\n" "\t"
    echo
  fi
}
export -f baseCount

#Output header line
echo -e "pop\tcontig\tpos\tref\tA\tC\tG\tT"

## Loop through VCF file for each contrig, position, and reference allele
## tail -n+25 skips 25 header lines (may be different per VCF file)
tail -n+25 all_indiv_concat.vcf | cut -f1,2,4 |
while read contig pos ref
do
  ## IDEA to speed this up: load the alignment once for each contig, rather than each damn SNP
  ## This does basically nothing to speed the whole thing up, but looks clever
  parallel baseCount {} $contig $pos $bases $ref ::: GROUP-{1..8}
done


rm "$bases"