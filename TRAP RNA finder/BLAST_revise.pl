#!/usr/bin/perl

use strict;
use Data::Dumper;

open (SPECIES, "<$ARGV[0]") || die "Could not open file:\t$ARGV[0]\n";
open (RESULTS, "<$ARGV[1]") || die "Could not open file:\t$ARGV[1]\n";

my ($file, $ext)=split(/\./, $ARGV[1]);
open (OUTFILE, ">>$file._with_species.txt");

my %species=();

my @species_file=<SPECIES>;
chomp(@species_file);

for (my $i=0; $i<=$#species_file; $i=$i+4){
	my @line=split(/,/, $species_file[$i]);
	my ($NC, $everything_else) = split(/\s/, $species_file[$i+2]);
	
	$species{$NC}=substr($line[0], 3);
	$species{$NC}=~s/^\s//;
	}

while(<RESULTS>){
	chomp;
	if ($_=~/(Query:\t)(NC_\d+\.\d+)\.(\d+)/i){
		print OUTFILE "\nQuery:\t".$species{$2}."\t$2.$3\n\n";
		}
		
	elsif ($_=~/^>/){
		my ($NC_base, $NC_revision, @tmp)=split(/\./,substr($_, 1));
		print OUTFILE ">".$NC_base.".".$NC_revision."\t".$species{$NC_base.".".$NC_revision}."\n";
		
		}
		
	elsif ($_=~/^(NC_\d+\.\d+\.\d+)/){
		my ($NC_base, $NC_revision, @tmp)=split(/\.|\t/, $_);
		print OUTFILE $species{$NC_base.".".$NC_revision}."\t".$_."\n";
		}
	else {
		print OUTFILE $_."\n";
		
		}
	}