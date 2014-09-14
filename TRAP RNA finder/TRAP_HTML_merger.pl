#!/usr/local/bin/perl

use strict;
my $file;
my @files=();
opendir (DIR,"/Users/eric_milliman/Desktop/shraddha_project/results/html/");# || die "Could not open Directory\n";

my @files = grep /^NC_\d+/, readdir(DIR);
print join("\n", @files)."\n";

open (HTML,">>TRAP_results_all.html");

print HTML <<HTML1;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html>
	<head>
		<script src="sorttable.js"></script>
		<link rel="stylesheet" type="text/css" href="mapper.css" />
	</head>
	<body>
	<table class="results sortable">
		<tr class="header">
			<td>ID</td>
			<td>Score</td>
			<td>F/R</td>
			<td>Start</td>
			<td>End</td>
			<td>Repeats</td>
			<td>Sequence</td>
		</tr>
HTML1

foreach my $file(@files){
	open (INFILE, "</Users/eric_milliman/Desktop/shraddha_project/results/html/$file") || die "Could not open file: $file\n";
	my $copy=0;
	
	while (<INFILE>){
		chomp;
		if ($copy==1){
			unless ($_=~/<\/table>/){
				print HTML $_."\n";
				}	
			else { last;	}
			}
		elsif ($_=~/<td class="ID">NC_\d+/){
			$copy=1;
			print HTML "<tr>\n$_\n";
			}
		}
	close (INFILE);
	}
		
	print HTML <<HTML2;
	
		</table>
	</body>
	</html>	
HTML2
				
				