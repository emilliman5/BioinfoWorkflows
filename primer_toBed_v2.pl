#!/usr/local/bin/perl

#script that will map oligonucleotide sequences (specifically the output from PCRtiler) and output bed file. Ideally used for small genome snapshots on IGV.

#08/14/2014 EJM

use strict;

open (MMTV, "<$ARGV[0]") || "Could not open: $ARGV[0]\n"; #input sequence

open(PRIMERS, "<$ARGV[1]")  || "Could not open: $ARGV[1]\n"; #primer sequences

open (BED, ">>$ARGV[2]"); #output file

my $genome_seq="";
my $seqID='';
print "$seqID\n";

while(<MMTV>){

	chomp;
	if(/^>/){
	
		$seqID=substr($_,1);
		}
	
	else{
		
		$genome_seq .=$_;
		
		}
	
	}
	
print BED "#ggfTags\n";

while(<PRIMERS>){
	
	next if (!/^\d/);
	chomp;
	my $left='';
	my $right='';
	
	my @line=split("\t",$_);
	
		if($genome_seq=~/($line[1])/){
			#print "Found a match\n";
			$left=$-[0];
			my $length=$left+length($line[1]);
			print BED $seqID."\t".$left."\t".$length."\tPrimer_".$line[0]."_F"."\t1\t+\n";
		}
		
		my $rc=reverse($line[2]);
		$rc=~tr/AGTC/TCAG/;
			
		if($genome_seq=~/($rc)/){
			$right=$+[0];
			my $length=$right-length($line[2]);
			print BED $seqID."\t".$length."\t".$right."\tPrimer_".$line[0]."_R"."\t1\t-\n";

			}
	
	}
	
		