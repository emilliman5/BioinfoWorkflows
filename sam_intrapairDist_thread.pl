#!/usr/bin/perl


use strict;
use File::Basename;
use threads;

my $number_of_threads=4;

my @suffixlist=(".sam", ".bam", ".bam.bai", ".fastq");
my $frag=$ARGV[0];
open (SAM, "<$ARGV[1]") || die "Could not open $ARGV[0]\n";

my $file=basename($ARGV[0], @suffixlist);

# my $nfr_file=$file."_NFR.sam";
# my $nuc_file=$file."_Nuc.sam";

open (NFR, ">>$ARGV[2]");
open (NUC, ">> $ARGV[3]");
open (DIST, ">>$ARGV[4]");

while (<SAM>){

	chomp;	
	if( $_=~/^\@/){
		print NFR $_."\n";
		print NUC $_."\n";
		}
	else {
	
		
	}
	
sub initThreads{
	my  @initThreads;
		foreach(my $i=1; $i<=$number_of_threads; $i++){
			push(@initThreads, $i);
			}
		return @initThreads;
		}
	
sub splitSam{
	
	my $id =threads->tid();
	my $line=shift;
	
	my @tmp=split(/\t/, $line);
		
		if ($tmp[8] >=0 ){
			print DIST $tmp[8]."\n";
		}
		
		if (abs($tmp[8]) <= $frag){
			print NFR $line."\n";
			}
			
		else {			
			print NUC $line."\n";
			}
		threads->exit();
		}