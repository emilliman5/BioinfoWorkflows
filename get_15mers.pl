#!/usr/local/bin/perl

use strict;

open (FILE," < chrM.fa");
open(OUT, ">> mtDNA_15mers.fastq");
my $seq_raw='';
while(<FILE>){
    chomp;
    unless($_=~/>/){
        $seq_raw .=$_;
    }
}

my $seq=$seq_raw.substr($seq_raw,0,18);

for(my $i=0;$i<length($seq)-18;$i++){
    my $start=$i+1;
    my $stop=$i+20;
    print OUT "\@mtDNA:".$start.":".$stop."\n".substr($seq,$i,20)."\n"."\+mtDNA:".$start.":".$stop."\n"."FFFFFFFFFFFFFFFFFFFF\n";
    }

$seq=reverse($seq_raw).substr(reverse($seq_raw),0,18);
$seq=~tr/ATCG/TAGC/;

for(my $i=0;$i<length($seq)-15;$i++){
    my $start=$i+1;
    my $stop=$i+21;
    print OUT "\@mtDNA_Rev_and_Comp:".$start.":".$stop."\n".substr($seq,$i,20)."\n"."\+mtDNA_Rev_and_Comp:".$start.":".$stop."\n"."FFFFFFFFFFFFFFFFFFFF\n";
    }