#!/usr/local/bin/perl

use strict;

open (INFILE, "<$ARGV[0]") || die "Could not open $ARGV[0]\n";
open (OUTFILE, ">>$ARGV[1]");

while (<INFILE>){
	chomp;
	
	my @line=split("\t", $_);
	my $exp=2**$line[4];
	print OUTFILE join("\t", @line[0..2], $exp, "+\n");
	}
