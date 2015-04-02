#!/usr/bin/perl
use strict;

my $outName = $ARGV[0];
my $pair1 = $ARGV[1];
my $pair2 = $ARGV[2];
my $minlength = $ARGV[3] || 25;
my $pair1FH;
my $pair2FH;
my $pair1outFH;
my $pair2outFH;
my $pairSingleFH;

print "pair1: $pair1\n";
print "pair2: $pair2\n";
print;

open ($pair1FH, '-|', 'zcat', "$pair1") || die "Cannot Open pair 1 File";
print "------\n";
print $pair1,"\n";
print "------\n";

<STDIN>;

open ($pair2FH, '-|', 'zcat', "$pair2") || die "Cannot Open pair 2 File";
print "------\n";
print $pair2,"\n";
print "------\n";

<STDIN>;

print "I am here\n";
<STDIN>;


my @pairedFiles = ($pair1FH,$pair2FH);

my %readCount;
my $totalRecords = 0;
for my $pairFH (@pairedFiles) {
  while (my $head = <$pairFH>) {
    <$pairFH>;
    <$pairFH>;
    <$pairFH>;
    $totalRecords++;
    my $key = substr($head, 0,  index($head, " "));
    my $value = $readCount{$key};
    if(defined $value) {
      $value++;
      $readCount{$key} = $value;
    }
    else {
      $readCount{$key} = 1;
    }
  }
  close($pairFH);
}

open ($pair1FH, '-|', 'zcat', "$pair1") || die "Cannot Open pair 1 File";
open ($pair2FH, '-|', 'zcat', "$pair2") || die "Cannot Open pair 2 File";
@pairedFiles = ($pair1FH,$pair2FH);
#open ($pair1outFH, "| gzip >$outName.pair1.fastq.gz") || die "Cannot Open pair out 1 File";
#open ($pair2outFH, "| gzip >$outName.pair2.fastq.gz") || die "Cannot Open pair out 2 File";
#open ($pairSingleFH, "| gzip >$outName.single.fastq.gz") || die "Cannot Open single File";
open ($pair1outFH, ">$outName.pair1.fastq") || die "Cannot Open pair out 1 File";
open ($pair2outFH, ">$outName.pair2.fastq") || die "Cannot Open pair out 2 File";
open ($pairSingleFH, ">$outName.single.fastq") || die "Cannot Open single File";

my @outputFiles = ($pair1outFH,$pair2outFH);

my $idx=0;
my $records=0;
my $pairs=0;
my $single=0;
my $discarded=0;
my %ignore;

print "Copying reads\n";
for my $pairFH (@pairedFiles) {
  my $outFH = $outputFiles[$idx];
  $idx++;

  while(my $head = <$pairFH>) {
    my $seq = <$pairFH>;
    my $head2 = <$pairFH>;
    my $qual = <$pairFH>;

    my $key = substr($head, 0, index($head, " "));

    my $value = $readCount{$key};
    my $fh;
    $records++;
    if($value == 1) {
      $fh = $pairSingleFH;
      $single++;
    }
    elsif($value > 1) {
      $fh = $outFH;
      $pairs++;
    }
    else {
      die "No value found for: $key\n";
    }

    if ($ignore{$key} != 1 && length($seq) >= $minlength) {
      print $fh $head;
      print $fh $seq;
      print $fh $head2;
      print $fh $qual;
    } else {
      $ignore{$key} = 1;
      $discarded++;
    }
#    print "\rTotal: $totalRecords, Records: $records, Pairs: $pairs, Single: $single                      ";
  }
  close($outFH);
  close($pairFH);
}
close($pairSingleFH);
print "Total: $totalRecords, Records: $records, Pairs: $pairs, Single: $single\n";
print "Discarded: $discarded reads (although possibly half this many)"
