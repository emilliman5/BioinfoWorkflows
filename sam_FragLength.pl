#!/usr/bin/perl


use strict;
use File::Basename;

my @suffixlist=(".sam", ".bam", ".bam.bai", ".fastq");

open (SAM, "<$ARGV[0]") || die "Could not open $ARGV[0]\n";

my $file=basename($ARGV[0], @suffixlist);
my $frag_file=$file."_fragLength.sam";

open (TXT, ">>$frag_file");

while (<SAM>){

	chomp;
	
	my @tmp=split(/\t/, $_);
		
		if ($tmp[8]>=0){
			print TXT $_."\n";
		}	
	}