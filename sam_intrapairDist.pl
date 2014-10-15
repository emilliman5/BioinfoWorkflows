#!/usr/bin/perl


use strict;
use File::Basename;

my @suffixlist=(".sam", ".bam", ".bam.bai", ".fastq");
my $frag=$ARGV[0];
my $nuc=180;

open (SAM, "<$ARGV[1]") || die "Could not open $ARGV[0]\n";

my $file=basename($ARGV[1], @suffixlist);

# my $nfr_file=$file."_NFR.sam";
# my $nuc_file=$file."_Nuc.sam";

open (NFR, ">>$ARGV[2]");
open (NUC, ">>$ARGV[3]");
open (OUT, ">>$ARGV[4]");
open (DIST, ">>$ARGV[5]");

while (<SAM>){

	chomp;
	
	if( $_=~/^\@/){
		print NFR $_."\n";
		print NUC $_."\n";
		}
	else {
	
		my @tmp=split(/\t/, $_);
		
		if ($tmp[8] >=0 ){
			print DIST $tmp[8]."\n";
		}
		
		if (abs($tmp[8]) <= $frag){
	
			print NFR $_."\n";
			}
			
		elsif (abs($tmp[8]) >= $nuc) {
			
			print NUC $_."\n";
			}
		
		else {
		
			print OUT $_."\n";
			}
		}
	}