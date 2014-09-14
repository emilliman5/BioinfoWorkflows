#!/usr/bin/perl

use strict;
use Net::FTP;
use Cwd;

open (SOURCE,">>Source_genomes.txt") || die "Could not create file: $!\n";
my $current_path=getcwd();
print $current_path."\n";

if (@ARGV){
	print "..."; ##Help and usage for obtaining genome sequences and annotations via this FTP module...
	}

my $ftp = Net::FTP->new("ftp.ncbi.nlm.nih.gov", Debug => 0)
      or die "Cannot connect to ftp.ncbi.nlm.nih.gov: $@";
$ftp->binary();
$ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;
   
$ftp->cwd("/genomes/Bacteria")
      or die "Cannot change working directory ", $ftp->message;
   
my @list=$ftp->ls()
      or die "ls failed ", $ftp->message;
      
my @sorted_list=sort(@list);
    
for (my $i=0; $i<=$#sorted_list; $i++){
	print $i."\t".$sorted_list[$i]."\n";
	}
      
print "Choose bacteria to acquire (separate multiple indices by a space): ";
my @bacteria=split(/\s/,<STDIN>); 
chomp(@bacteria);

my $storage_path=$current_path."/genomes/";

my $directory_path=mkdir "$storage_path", 0777 unless -d "$storage_path";

foreach my $selection(@bacteria){	
	my $directory_path=mkdir "$storage_path"."$sorted_list[$selection]", 0777 unless -d "$storage_path"."$sorted_list[$selection]";
	print SOURCE $storage_path.$sorted_list[$selection]."/\n";
	my @files=$ftp->ls("$sorted_list[$selection]/");
	
	foreach my $file(@files){
		if($file=~/(fna)|(gff)|(ptt)/){
			my $touch=qx(touch $storage_path$file);
			my $get=$ftp->get("$file", "$storage_path"."/"."$file") 
			 or die "Cannot get file ", $ftp->message;
		}
	}
}

$ftp->quit;

