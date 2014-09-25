#!/usr/bin/perl


use strict;
use File::Basename;

my @suffixlist=(".sam", ".bam", ".bam.bai", ".fastq");

open (SAM, "<$ARGV[0]") || die "Could not open $ARGV[0]\n";

my $file=basename($ARGV[0], @suffixlist);

my $nfr_file=$file."_NFR.sam";
my $nuc_file=$file."_Nuc.sam";

open (NFR, ">>$nfr_file");
open (NUC, ">> $nuc_file");

while (<SAM>){

	chomp;
	
	if( $_=~/^\@/){
		print NFR $_."\n";
		print NUC $_."\n";
		}
	else {
	
		my @tmp=split(/\t/, $_);
		
		if (abs($tmp[8]) <= 130){
	
			print NFR $_."\n";
			}
			
		else {
			
			print NUC $_."\n";
			}
		}
	}