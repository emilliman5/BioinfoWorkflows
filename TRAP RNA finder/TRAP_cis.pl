#!/usr/bin/perl

#TRAP cis-element finder script

use strict;

open (ANNONTATIONS, "<$ARGV[0]") || die "Error: could not find/open $ARGV[0]\n";

my @base_file = split( /\//, $ARGV[0]); 

my ($file, $ext) = split(/\./, $base_file[-1]);

open (CISFILE, ">>$file"."_cis.txt");

while (<ANNONTATIONS>){
	chomp;
	my @line =split(/\t/, $_);
	
	for my $index (10, 16){
	
		if ($line[$index]=~/tandem/i){
			
			if ($line[$index]=~/upstream/i){
							
				if ($line[$index-1] <= 150){
					print CISFILE join("\t", @line[0..4], @line[$index-5..$index],"\n");
					}
				}
			
			elsif($line[$index]=~/5' overlap/i){
					print CISFILE join("\t", @line[0..4], @line[$index-5..$index],"\n");
					}
				
			elsif ($line[$index]=~/within/i){
				if ($line[$index-1] <= 50){
					print CISFILE join("\t", @line[0..4], @line[$index-5..$index],"\n");
					}
				}
			
			}
		}
	}
		

