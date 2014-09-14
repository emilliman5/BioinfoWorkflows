#!/usr/bin/perl

##ToDo:		Search for Rho-dependent terminators High-C low G, and single strandedness... 60-100nt rut site.
##			Trp codon UGG anticodon ACC
##			Add non-complete genome hack (i.e., B. strearo contigs)

##			PSSM conversion -- Create frequency matrix then... score(j,i)=log m(j,i)/f(j) m(j,i) = frequency of character j at position i 
##			f(j) = frequency overall frequency of j.

use strict;
use Getopt::Long;
use TRAPAlgorithm3;
use BinarySearch;
use Data::Dumper;
use IO::File;

my ($paths, $help, $cut_off, $fuzznuc_patterns, $tbox, $tbox_cutoff, $CM_model) ;
my @matrix;
GetOptions ( 	"help|h" => \$help,
				"matrix|m=s" => \@matrix,			##PSSM files for promoter predictions
				"genome|g=s" => \$paths,			##genome sequence and annotation path file (from FTP script);
				"fuzznuc|f=s" => \$fuzznuc_patterns,
				"cm_model|r=s" => \$CM_model);
			
my @gffVersion3Fields=qw(SeqID Source Feature Start End Score Strand Frame Attribute);
my @transTermFields=qw(ID Start Stop Conf HP Tail Distance);

if ($help || !$paths || !$fuzznuc_patterns){
	exit &show_help();
	}

open (BACTERIA, "<$paths") || die "Could not find or open $paths\n";

my @GenomesToProcess=<BACTERIA>;
chomp @GenomesToProcess;
my %Genomes=();

foreach my $genome(@GenomesToProcess){
	chomp;

	opendir(DIR, $genome) || die "Could not open directory $genome\n";
	my @files=readdir(DIR);
	
	foreach my $File(@files){
		next if ($File=~/^\./);
		my ($name, $ext)=split(/\./, $File);
		$Genomes{$name}{$ext}=$genome.$File;
		}
	}

foreach my $genome(keys(%Genomes)){
	my %FinalROIs=();
	my %GenomeSeq=();
	my %ROIAnnotations=();
	my %ROItranscript=();
	my @fuzznuc=();
	
	my $dir_path=mkdir "results", 0777 unless -d "results";
	my $directory_path=mkdir "results/$genome", 0777 unless -d "results/$genome/";
	my $directory_path_html=mkdir "results/html", 0777 unless -d "results/html/";
	print "\nStarting TRAP search on $genome...\n";
	
	open (my $seqfile,"<$Genomes{$genome}{'fna'}") || die "Could not open $Genomes{$genome}\n";
	open (my $genomeannotations,"<$Genomes{$genome}{'gff'}") || die "Could not open $Genomes{$genome}\n";
	open (my $html,">>results/html/$genome"."_parsed.html");
	open (my $fuzznuc_parse,">>results/$genome/$genome"."_fuzznuc_merged.gff");
	open (my $seq_out,">>results/$genome/$genome"."_seq_fuzznuc.fa");
	open (my $annotations,">>results/$genome/$genome"."_annotatations.txt");
	open (my $intergenic, ">>results/$genome/$genome"."_intergenic_annotations.txt");
	open (my $transcript, ">>results/$genome/$genome"."_in_silico_transcript.txt");
	open (my $intergenic_transcript, ">>results/$genome/$genome"."_intergenic_in_silico_transcript.txt");
	
	foreach my $file($annotations, $intergenic){
		print $file "Fuzznuc ID\tStrand\tStart\tStop\tScore\tFeature 1\tStart\tStop\tStrand\tDistance from Binding site\tOrientation\tFeature 2\tStart\tStop\tStrand\tDistance from Binding site\tOrientation\n";
		}
	foreach my $file($transcript, $intergenic_transcript){
		print $file "Fuzznuc ID\tStrand\tStart\tStop\tScore\tPromoter ID\tDistance\tPvalue\tTbox ID\tDistance\tUpstream Terminator ID\tStart\tStop\tConf\tHP score\tTail Score\tDistance From Binding site\tDownstream Terminator ID\tStart\tStop\tConf\tHP score\tTail Score\tDistance From Binding site\n";
		}
		
	my $header;
	
	while(<$seqfile>){
		chomp;
		if (/^>.*(NC_\d+)(\.\d)/){
			$header=$1.$2;
			$GenomeSeq{$header}='N';
			}
		elsif(/^>/){
			$_=~s/>//;
			$header=$_;
			$GenomeSeq{$header}='N';
			}
		else{
			$_=~tr/a-z/A-Z/;
			$GenomeSeq{$header}.=$_;
			}
		}
		
######################################################################################################
##						Genome GC content for Promoter search
######################################################################################################	
	my @sequence=split(//,$GenomeSeq{$header});
	my $GC_content = calculate_GC(@sequence) ;
	my $AT_content = (1-$GC_content);
	my $Patser_content="a:t $AT_content g:c $GC_content";
	
######################################################################################################
######################################################################################################
##						TRAP binding site aearch using Fuzznuc
######################################################################################################	

	&fuzznuc(\$Genomes{$genome}{'fna'}, \$fuzznuc_patterns, \@fuzznuc);
	print "Fuzznuc search finished return: $?\n";	
	
######################################################################################################
######################################################################################################
##						Terminator Search using TransTermHP
######################################################################################################
	my %Terminators=();
	print "Starting intrinsic terminator search with transtermHP...\n";
	my %termCoordinates=&TransTermHP( \$Genomes{$genome}{'fna'}, \$Genomes{$genome}{'ptt'}, \%Terminators);

######################################################################################################
######################################################################################################
##						Promoter Search using Patser v3a 2008
######################################################################################################
	my %Promoters=();
	my %promCoordinates;
	print "Starting promoter search with Patser-v3e...\n";
	
	%promCoordinates=&Patser( \$Genomes{$genome}{'fna'}, \@matrix, \%Promoters, $Patser_content);
	open (TEST, ">>Promter_dump.txt");
	print TEST Dumper(%Promoters)."\n";
######################################################################################################
######################################################################################################
##									Tbox Search using Infernal
##				To run infernal, a covariance model will need to be generated first,
##	using the cmbuild, cmalign, and cmcalibrate, the alignments can be retreveived from the Rfam database
##
######################################################################################################
	my %Tboxes=();
	print "Starting Tbox search with Infernal Covariance model...\n";
	my %TboxCoordinates=&Infernal( \$Genomes{$genome}{'fna'}, \$CM_model, \%Tboxes);

######################################################################################################

	my %preROIs=();
	my %ROIs=();
	my %GenomeFeatures=();
	my @Annotations=<$genomeannotations>;
	gff3Parser(\@fuzznuc, \%preROIs, \@gffVersion3Fields);
	gff3Parser(\@Annotations, \%GenomeFeatures, \@gffVersion3Fields);
	my %sortedGenomeFeatures=();
	
	foreach my $SeqID(keys(%preROIs)){
		foreach my $ROI(keys(%{$preROIs{$SeqID}})){
			%{$ROIs{$SeqID}{$preROIs{$SeqID}{$ROI}{"Strand"}}{$ROI}}=%{$preROIs{$SeqID}{$ROI}};
			}
		}
######################################################################################################
##						Time to coalesce the results from fuzznuc
######################################################################################################	

	my ($ROIBuild, $ROIBuildStart, $ROIBuildEnd, $ROIBuildStrand, $ROIBuildSeqID, $strand);
	foreach my $SeqID(keys(%ROIs)){
		my $reset=1;
		foreach $strand (keys(% {$ROIs{$SeqID}})){
		
			foreach my $ROI (sort { ${$ROIs{$SeqID}{$strand}}{$a}->{"Start"} <=> ${$ROIs{$SeqID}{$strand}}{$b}->{"Start"} } keys(% {$ROIs{$SeqID}{$strand}})) {
				
				if ($reset){
					$reset=0;
					$ROIBuild=$ROI;
					start_end_strand($ROIs{$SeqID}{$strand}{$ROIBuild}{'Start'}, $ROIs{$SeqID}{$strand}{$ROIBuild}{'End'}, \$ROIBuildStart, \$ROIBuildEnd);
					}
					
				else{
				
					my ($CurrentStart, $CurrentEnd);
					start_end_strand($ROIs{$SeqID}{$strand}{$ROI}{'Start'}, $ROIs{$SeqID}{$strand}{$ROI}{'End'}, \$CurrentStart, \$CurrentEnd);
					my $left_between=between($ROIBuildStart, $ROIBuildEnd, $CurrentStart);
					my $right_between=between($ROIBuildStart, $ROIBuildEnd, $CurrentEnd);
					
					if  ((($left_between || $right_between) || ( between($CurrentStart, $CurrentEnd, $ROIBuildStart) && 
							between($CurrentStart, $CurrentEnd, $ROIBuildEnd)))) { 
		
						if(	between($CurrentStart, $CurrentEnd, $ROIBuildStart) && 
												between($CurrentStart, $CurrentEnd, $ROIBuildEnd)){
							$ROIBuildStart=$CurrentStart;
							$ROIBuildEnd=$CurrentEnd;
							}	
						elsif($left_between && !$right_between){
							$ROIBuildEnd=$CurrentEnd;
							}
						elsif($right_between && !$left_between){
							$ROIBuildStart=$CurrentStart;
							}
						}
					else {
						$FinalROIs{$SeqID}{$ROIBuild}{"Start"}=$ROIBuildStart;		
						$FinalROIs{$SeqID}{$ROIBuild}{"End"}=$ROIBuildEnd;
						$FinalROIs{$SeqID}{$ROIBuild}{"Strand"}=$strand;
						##restarts the build process....
						$ROIBuild=$ROI;
						start_end_strand($ROIs{$SeqID}{$strand}{$ROIBuild}{'Start'}, $ROIs{$SeqID}{$strand}{$ROIBuild}{'End'}, \$ROIBuildStart, \$ROIBuildEnd);
						}
					}
				}
				$FinalROIs{$SeqID}{$ROIBuild}{"Start"}=$ROIBuildStart;		
				$FinalROIs{$SeqID}{$ROIBuild}{"End"}=$ROIBuildEnd;
				$FinalROIs{$SeqID}{$ROIBuild}{"Strand"}=$strand;	
			}
		}
	
######################################################################################################
##			Now time to extract sequence of fuzznuc search and score the results...
######################################################################################################	

	&html_head(\$html, $Genomes{$genome}{'fna'});
		
		foreach my $SeqID(keys(%FinalROIs)){
			my @sortedCoordinates;
			push @{ $sortedGenomeFeatures{$SeqID} },	sort { ${$GenomeFeatures{$SeqID}}{$a}->{"Start"} <=> ${$GenomeFeatures{$SeqID}}{$b}->{"Start"} } keys(%{$GenomeFeatures{$SeqID}});
			foreach my $feat(@{$sortedGenomeFeatures{$SeqID}}){
				push @sortedCoordinates, $GenomeFeatures{$SeqID}{$feat}{'Start'};
				}
				
			foreach my $roi(keys(%{$FinalROIs{$SeqID}})){
				
				my $midpoint=int(($FinalROIs{$SeqID}{$roi}{'End'}-$FinalROIs{$SeqID}{$roi}{'Start'})/2)+$FinalROIs{$SeqID}{$roi}{'Start'};
				my $left=$FinalROIs{$SeqID}{$roi}{"Start"};
				my $right=$FinalROIs{$SeqID}{$roi}{"End"};
				my $strand=$FinalROIs{$SeqID}{$roi}{"Strand"};
				my $match_length=abs($left-$right);
				my $seq_match=substr($GenomeSeq{$SeqID}, ($left-5), ($match_length+10));
				
				if ($strand eq "-"){
					my $tmp=$seq_match;
					$seq_match=&rev_comp($tmp);
					}	
					
				my ($successFlag,$returnedSequence,$score,$count,$startShift,$endShift) = TRAPalgorithm3::algorithmScore($seq_match);			
				
				if (!$successFlag){
					$seq_match=substr($GenomeSeq{$SeqID}, ($left-25), ($match_length+50));
					
					if ($strand eq "-"){
						my $tmp=$seq_match;
						$seq_match=&rev_comp($tmp);
						}
						
					($successFlag,$returnedSequence,$score,$count,$startShift,$endShift) = TRAPalgorithm3::algorithmScore($seq_match);
					}
					
				$FinalROIs{$SeqID}{$roi}{'Score'}=$score;

				print $html "\n\t\t<tr>\n\t\t\t<td class=\"ID\">$roi</td>\n\t\t\t<td class=\"score\">$score</td>\n\t\t\t<td class=\"direction\">$strand</td>\n\t\t\t<td class=\"start\">$left</td>\n\t\t\t<td class=\"end\">$right</td>\n\t\t\t<td class=\"repeats\">$count</td>\n\t\t\t<td class=\"sequence\">$returnedSequence</td>\n\t\t</tr>";
				print $fuzznuc_parse qq/$SeqID\tfuzznuc\tTRAP_pattern_match\t$left\t$right\t$score\t$strand\t.\tID=$roi; NAME=$roi\n/;	

######################################################################################################
##			Time to annontate nearby promoters, terminators and T-boxes...
######################################################################################################	

			my $terminator_index=BinarySearch::BinSearch($midpoint, \@{$termCoordinates{$strand}});
			my $term_index=int($terminator_index);
			my $promoter_index=BinarySearch::BinSearch($midpoint, \@{$promCoordinates{$strand}});
			my $prom_index=int($promoter_index);
			my $T_box_index=BinarySearch::BinSearch($midpoint, \@{$TboxCoordinates{$strand}});
			my $tbox_index=int($T_box_index);
			
			
			my $prom_index2= &upstream_element(\%ROItranscript, $roi, $prom_index, $strand, \%Promoters, $midpoint, "Promoter");
			my $tbox_index2= &upstream_element(\%ROItranscript, $roi, $tbox_index, $strand, \%Tboxes, $midpoint, "TBox");			
			my $term_index2= &upstream_element(\%ROItranscript, $roi, $term_index, $strand, \%Terminators, $midpoint, "Terminator");
			my $term_index3;
			
			if($term_index2 == "NULL"){
				$term_index3=$term_index;
				}
			elsif($strand eq "-"){
				$term_index3=$term_index2-1;
				}
			else{	$term_index3=$term_index2+1;	}
			
			
			if ($strand eq "-"){
				
				if (&smallest_difference($midpoint, $Terminators{$strand}{$term_index3}{'Start'}, $Terminators{$strand}{$term_index3}{'Stop'}) <= 500 && $Terminators{$strand}{$term_index3}{"Start"}-$midpoint <=0 ){
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}=$Terminators{$strand}{$term_index3};
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}{"Distance"}=abs($midpoint-$Terminators{$strand}{$term_index3}{'Start'});
					}
				
				elsif (&smallest_difference($midpoint, $Terminators{$strand}{$term_index3-1}{'Start'}, $Terminators{$strand}{$term_index3-1}{'Stop'}) <= 500 && $Terminators{$strand}{$term_index3-1}{"Start"}-$midpoint <=0 ){
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}=$Terminators{$strand}{$term_index3-1};
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}{"Distance"}=abs($midpoint-$Terminators{$strand}{$term_index3-1}{'Start'});
					}
				
				}
				
			elsif ($strand eq "+") {
				
				if (&smallest_difference($midpoint, $Terminators{$strand}{$term_index3}{'Start'}, $Terminators{$strand}{$term_index3}{'Stop'}) <= 500 && $Terminators{$strand}{$term_index3}{"Start"}-$midpoint >=0 ){
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}=$Terminators{$strand}{$term_index3};
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}{"Distance"}=abs($midpoint-$Terminators{$strand}{$term_index3}{'Start'});
					}
					
				elsif (&smallest_difference($midpoint, $Terminators{$strand}{$term_index3+1}{'Start'}, $Terminators{$strand}{$term_index3+1}{'Stop'}) <= 500 && $Terminators{$strand}{$term_index3+1}{"Start"}-$midpoint <=0 ){
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}=$Terminators{$strand}{$term_index3+1};
					$ROItranscript{$roi}{"Terminator"}{"Downstream"}{"Distance"}=abs($midpoint-$Terminators{$strand}{$term_index3+1}{'Start'});
					}
				
				}
				
			else { die "No Directionality given\n";	}
			
			my %Txn_elements= {	"Promoter" => $ROItranscript{$roi}{"Promoter"}{"Upstream"}{"Distance"},
								"Tbox" => $ROItranscript{$roi}{"TBox"}{"Upstream"}{"Distance"},
								"Term_up" => $ROItranscript{$roi}{"Termiator"}{"Upstream"}{"Distance"},
								};
			
			my @sorted= sort {$Txn_elements{$b} <=> $Txn_elements{$a}} keys %Txn_elements;
			my ($seq_dist, $upstream_dist, $downstream_dist, $left_coord);
			
			if ($sorted[0]>1){
				$upstream_dist=$sorted[0];
				}
			else {
				$upstream_dist=500;
				}
			
			if($ROItranscript{$roi}{"Termiator"}{"Downstream"}{"Distance"}>1){
				$downstream_dist=$ROItranscript{$roi}{"Termiator"}{"Downstream"}{"Distance"};
				}
			else {
				$downstream_dist=500;
				}
				
			$seq_dist=$upstream_dist+$downstream_dist;
			
			$left_coord=$left-$upstream_dist;
			
			my $sequence=substr($GenomeSeq{$SeqID}, $left_coord, $seq_dist);
			if ($strand eq "-"){
				my $tmp=$sequence;
				$sequence=&rev_comp($tmp);
				}
			print $transcript $roi."\t".$FinalROIs{$SeqID}{$roi}{'Strand'}."\t".$FinalROIs{$SeqID}{$roi}{'Start'}."\t".$FinalROIs{$SeqID}{$roi}{'End'}."\t".$FinalROIs{$SeqID}{$roi}{'Score'}."\t";				##	Fix this it does not R&C minus strands results
			print $seq_out qq/>$roi\tStrand:$FinalROIs{$SeqID}{$roi}{"Strand"}\t$left\t$right\tTRAPalgorithm3 Score: $score\n/."$sequence\n";																	##

			if (exists $ROItranscript{$roi}{"Promoter"}{"Upstream"}{"Distance"}){
				print $transcript $ROItranscript{$roi}{"Promoter"}{"Upstream"}{"Matrix"}."_".$ROItranscript{$roi}{"Promoter"}{"Upstream"}{"ID"}."\t".$ROItranscript{$roi}{"Promoter"}{"Upstream"}{"Distance"}."\t".$ROItranscript{$roi}{"Promoter"}{"Upstream"}{"Pvalue"}."\t";
				}
			else { print $transcript "\t\t\t";	}
			
			if (exists $ROItranscript{$roi}{"TBox"}{"Upstream"}{"Distance"}){
				print $transcript $ROItranscript{$roi}{"TBox"}{"Upstream"}{"ID"}."\t".$ROItranscript{$roi}{"TBox"}{"Upstream"}{"Distance"}."\t";
				}
			else { print $transcript "\t\t";	}
			
			
			foreach my $side("Upstream", "Downstream"){
				if ($ROItranscript{$roi}{"Terminator"}{$side}{"Distance"}<=500){
					foreach my $field(@transTermFields){
						print $transcript $ROItranscript{$roi}{"Terminator"}{$side}{$field}."\t";	
						}
					}
				else {
					print $transcript "\t\t\t\t\t\t\t";
					}
				}
				
			print $transcript "\n";
			
######################################################################################################
##			Now time to annotate the putative binding sites...
######################################################################################################	
						
			my ($dist1, $dist2)= (100000,100000);
			my ($current_match1, $current_match2);
			
			my $index=BinarySearch::BinSearch($midpoint, \@sortedCoordinates);
				
			foreach my $genome_feat (@{$sortedGenomeFeatures{$SeqID}}[int($index-5)..int($index+5)]){
				my $gstart=$GenomeFeatures{$SeqID}{$genome_feat}{"Start"};
				my $gend=$GenomeFeatures{$SeqID}{$genome_feat}{"End"};
				
					if ($genome_feat ne $current_match2 && (abs($midpoint-$gstart)<$dist1 || abs($midpoint-$gend)<$dist1)){
						$ROIAnnotations{$roi}{"1"}[0]=$genome_feat;
						$dist2=$dist1;
						$current_match2=$current_match1;
						$ROIAnnotations{$roi}{"2"}[0]=$current_match2;
						$ROIAnnotations{$roi}{"2"}[1]=$dist2;
						$current_match1=$genome_feat;
						$dist1=smallest_difference($midpoint, $gstart, $gend);
						$ROIAnnotations{$roi}{"1"}[1]=$dist1;
						}
					elsif ($genome_feat ne $current_match1 && (abs($midpoint-$gstart)<$dist2 || abs($midpoint-$gend)<$dist2)){
						$ROIAnnotations{$roi}{"2"}[0]=$genome_feat;
						$current_match2=$genome_feat;
						$dist2=smallest_difference($midpoint, $gstart, $gend);
						$ROIAnnotations{$roi}{"2"}[1]=$dist2;
						}
					}
							
				foreach my $feat (keys(%{$ROIAnnotations{$roi}})){
				
					my $genome_feat=$ROIAnnotations{$roi}{$feat}[0];
					my $gleft=$GenomeFeatures{$SeqID}{$genome_feat}{"Start"};
					my $gright=$GenomeFeatures{$SeqID}{$genome_feat}{"End"};
					my $gstrand=$GenomeFeatures{$SeqID}{$genome_feat}{"Strand"};
					my $left_between=between($gleft, $gright, $left);
					my $right_between=between($gleft, $gright, $right);
					
					if($left_between && $right_between){
						if($gstrand eq "+"){
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "within, tandem orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "within, divergent orientation";
								}
							}
						else{ ##$gstrand eq "-"
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "within, divergent orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "within, tandem orientation";
								}
							}
						}
						
					elsif($left_between){
						if($gstrand eq "+"){
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "3' overlap, tandem orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "3' overlap, convergent orientation";
								}
							}
						else{ ##$gstrand eq "-"
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "5' overlap, divergent orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "5' overlap, tandem orientation";
								}
							}
						}	
					elsif($right_between){
						if($gstrand eq "+"){
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "5' overlap, tandem orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "5' overlap, divergent orientation";
								}
							}
						else{ ##$gstrand eq "-"
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "3' overlap, convergent orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "3' overlap, tandem orientation";
								}
							}
						}
					elsif($left > $gright){
						if($gstrand eq "+"){
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "downstream, tandem orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "downstream, convergent orientation";
								}
							}
						else{ ##$gstrand eq "-"
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "upstream, divergent orientation";
								}
							else{	##$strand eq "-"
								push @{$ROIAnnotations{$roi}{$feat}}, "upstream, tandem orientation";
								}
							}
						}
						elsif($left < $gright){
							if($gstrand eq "+"){
								if($strand eq "+"){
									push @{$ROIAnnotations{$roi}{$feat}}, "upstream, tandem orientation";
									}
								else{
									push @{$ROIAnnotations{$roi}{$feat}}, "upstream, divergent orientation";
									}
								}
							else{ 		##$gstrand eq "-"
								if($strand eq "+"){
									push @{$ROIAnnotations{$roi}{$feat}}, "dowstream, convergent orientation";
									}
								else{	##$strand eq "-"
									push @{$ROIAnnotations{$roi}{$feat}}, "dowstream, tandem orientation";
									}
								}
							}
					elsif($right > $gleft){
						if($gstrand eq "+"){
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "downstream, tandem orientation";
								}
							else{
								push @{$ROIAnnotations{$roi}{$feat}}, "downstream, convergent orientation";
								}
							}
						else{		##$gstrand eq "-"
							if($strand eq "+"){
								push @{$ROIAnnotations{$roi}{$feat}}, "upstream, divergent orientation";
								}
							else{	##$strand eq "-"
								push @{$ROIAnnotations{$roi}{$feat}}, "upstream, tandem orientation";
								}
							}
						}
					}
				}
			}

		&html_end(\$html);
		
		foreach my $SeqID (keys(%FinalROIs)){
		
			foreach my $ROI(keys(%{$FinalROIs{$SeqID}})){
		
				print $annotations "$ROI\t$FinalROIs{$SeqID}{$ROI}{'Strand'}\t$FinalROIs{$SeqID}{$ROI}{'Start'}\t$FinalROIs{$SeqID}{$ROI}{'End'}\t$FinalROIs{$SeqID}{$ROI}{'Score'}\t";
				my $IntergenicString="$ROI\t$FinalROIs{$SeqID}{$ROI}{'Strand'}\t$FinalROIs{$SeqID}{$ROI}{'Start'}\t$FinalROIs{$SeqID}{$ROI}{'End'}\t$FinalROIs{$SeqID}{$ROI}{'Score'}\t";
				my $Intergenic=1;
				foreach my $feat (keys(%{$ROIAnnotations{$ROI}})){
					my $gene=$ROIAnnotations{$ROI}{$feat}[0];
					print $annotations "$gene\t$GenomeFeatures{$SeqID}{$gene}{'Start'}\t$GenomeFeatures{$SeqID}{$gene}{'End'}\t$GenomeFeatures{$SeqID}{$gene}{'Strand'}\t$ROIAnnotations{$ROI}{$feat}[1]\t$ROIAnnotations{$ROI}{$feat}[2]\t";
					
					if ($ROIAnnotations{$ROI}{$feat}[2]=~/(upstream|dowstream)/gi){
						$IntergenicString .="$gene\t$GenomeFeatures{$SeqID}{$gene}{'Start'}\t$GenomeFeatures{$SeqID}{$gene}{'End'}\t$GenomeFeatures{$SeqID}{$gene}{'Strand'}\t$ROIAnnotations{$ROI}{$feat}[1]\t$ROIAnnotations{$ROI}{$feat}[2]\t";
						}	
					else{
						$Intergenic=0;
						}
					}
				print $annotations "\n";
				
				if($Intergenic){
					print $intergenic $IntergenicString."\n";
					print $intergenic_transcript "$ROI\t$FinalROIs{$SeqID}{$ROI}{'Strand'}\t$FinalROIs{$SeqID}{$ROI}{'Start'}\t$FinalROIs{$SeqID}{$ROI}{'End'}\t$FinalROIs{$SeqID}{$ROI}{'Score'}\t";
					
						if (exists $ROItranscript{$ROI}{"Promoter"}{"Upstream"}{"Distance"}){
							print $intergenic_transcript $ROItranscript{$ROI}{"Promoter"}{"Upstream"}{"Matrix"}."_".$ROItranscript{$ROI}{"Promoter"}{"Upstream"}{"ID"}."\t".$ROItranscript{$ROI}{"Promoter"}{"Upstream"}{"Distance"}."\t".$ROItranscript{$ROI}{"Promoter"}{"Upstream"}{"Pvalue"}."\t";
							}
						else { print $intergenic_transcript "\t\t\t";	}
						
						if (exists $ROItranscript{$ROI}{"TBox"}{"Upstream"}{"Distance"}){
							print $intergenic_transcript $ROItranscript{$ROI}{"TBox"}{"Upstream"}{"ID"}."\t".$ROItranscript{$ROI}{"TBox"}{"Upstream"}{"Distance"}."\t";
							}
						else { print $intergenic_transcript "\t\t";	}
						
						
						foreach my $side("Upstream", "Downstream"){
							if ($ROItranscript{$ROI}{"Terminator"}{$side}{"Distance"}<=500){
								foreach my $field(@transTermFields){
									print $intergenic_transcript $ROItranscript{$ROI}{"Terminator"}{$side}{$field}."\t";	
									}
								}
							else {
								print $intergenic_transcript "\t\t\t\t\t\t\t";
								}
							}
							
						print $intergenic_transcript "\n";
				}
			}
		}
	closedir(DIR);
}

####################################################################################################
#											SubRoutines											   #
####################################################################################################
sub upstream_element{
	
	my ($roi_transcript_ref, $roi, $index, $strand, $element_hashref, $midpoint, $element)=@_;
	my %Txn_feature=%{$element_hashref};
	my $return_index="NULL";
	my $distance=500;
	
	if ($strand eq "-"){
		if (&smallest_difference($midpoint, $Txn_feature{$strand}{$index}{'Start'}, $Txn_feature{$strand}{$index}{'Stop'}) <= $distance && $midpoint - $Txn_feature{$strand}{$index}{"Stop"} <= 0 ){
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}=$Txn_feature{$strand}{$index};
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}{"Distance"}=abs($midpoint-$Txn_feature{$strand}{$index}{'Stop'});
			}
		elsif (&smallest_difference($midpoint, $Txn_feature{$strand}{$index+1}{'Start'}, $Txn_feature{$strand}{$index+1}{'Stop'}) <= $distance && $Txn_feature{$strand}{$index+1}{"Stop"}-$midpoint >=0 ){
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}=$Txn_feature{$strand}{$index+1};
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}{"Distance"}=abs($midpoint-$Txn_feature{$strand}{$index+1}{'Stop'});
			$return_index=$index+1;
			}
		}		

	elsif ($strand eq "+"){
		if (&smallest_difference($midpoint, $Txn_feature{$strand}{$index}{'Start'}, $Txn_feature{$strand}{$index}{'Stop'}) <= $distance && $midpoint - $Txn_feature{$strand}{$index}{"Stop"} >=0 ){
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}=$Txn_feature{$strand}{$index};
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}{"Distance"}=abs($midpoint-$Txn_feature{$strand}{$index}{'Stop'});
		}
		elsif (&smallest_difference($midpoint, $Txn_feature{$strand}{$index-1}{'Start'}, $Txn_feature{$strand}{$index-1}{'Stop'}) <= $distance && $Txn_feature{$strand}{$index-1}{"Stop"}-$midpoint >=0 ){
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}=$Txn_feature{$strand}{$index-1};
			${$roi_transcript_ref}{$roi}{$element}{"Upstream"}{"Distance"}=abs($midpoint-$Txn_feature{$strand}{$index-1}{'Stop'});
			$return_index=$index-1;
			}
	}
	
	else {
	
			die "No directionality given\n";
		}
		return $return_index;
				
}				

sub fuzznuc {

#####################################################################################################
#																									#
#																									#
#								Fuzznuc	as part of the EMBOSS Suite									#
#																									#		
#																									#
#####################################################################################################
	my $sequence=${$_[0]};
	my $pattern='@'.${$_[1]};
	$sequence=~/(NC_\d{6})\.(\w+)/;
	@{$_[2]}=qx(fuzznuc -sequence $sequence -pattern $pattern -complement Y -rformat2 gff -outfile stdout);
	return;
}

sub TransTermHP {
#####################################################################################################
#																									#
#											TransTermHP												#
#																									#
#		C. Kingsford, K. Ayanbule and S.L. Salzberg. Rapid, accurate, computational discovery of 	#
#	Rho-independent transcription terminators illuminates their relationship to DNA uptake. Genome	#
#										Biology 8:R22 (2007).										#		
#																									#
#####################################################################################################
	my $sequence=${$_[0]};
	my $ptt_file=${$_[1]};
	my $Terminator_ref=$_[2];
	my %coordinates;
	my ($index_p, $index_m)=(0,0);
	my $transterm=qx(terminator_pred/transterm_hp_v2.06/transterm $sequence $ptt_file -p terminator_pred/transterm_hp_v2.06/expterm.dat --all-context);
	
	my @terms=split(/\n/,$transterm);
	my @lineToProcess;
	foreach my $line(@terms){
		if ($line=~/term\s\d+/gi){
			@lineToProcess=split (/\s{1,}/, $line);
			my $ID=$lineToProcess[1].$lineToProcess[2];		
			my ($left, $right)=&start_end_strand($lineToProcess[3], $lineToProcess[5]);		##double check that the left ("start") coordinate is lower than the right ("stop") coordinate.
			if($lineToProcess[6] eq "-"){
				${$Terminator_ref}{$lineToProcess[6]}{$index_m}= {
					"ID"	=> $ID,
					"Start" => $lineToProcess[3],
					"Stop"	=> $lineToProcess[5],
					"HP"	=> $lineToProcess[9],
					"Conf"	=> $lineToProcess[8],
					"Tail"	=> $lineToProcess[10]
					};
				push @{$coordinates{$lineToProcess[6]}}, $left;
				$index_m++;
			}
			else {
				${$Terminator_ref}{$lineToProcess[6]}{$index_p}= {
					"ID"	=> $ID,
					"Start" => $lineToProcess[3],
					"Stop"	=> $lineToProcess[5],
					"HP"	=> $lineToProcess[9],
					"Conf"	=> $lineToProcess[8],
					"Tail"	=> $lineToProcess[10]
					};
				push @{$coordinates{$lineToProcess[6]}}, $left;
				$index_p++;
			}
		}

	}
	#print Dumper(%results)."\n";
	print "TranstermHP search finished return: $?\n";
	return %coordinates;
}

sub Patser {
#####################################################################################################
#											Patser													#
#									Part of the Consensus Package									#
#										From the Stormo Lab											#
#####################################################################################################
	my $sequence=${$_[0]};
	my @matrix_files=@{$_[1]};
	my $promoter_ref=$_[2];
	my $gc_content=$_[3];
	my %coordinates;
	my %Temp_promoter=();
	my ($index_p, $index_m, $index)=(0,0,0);
	open (TMP ,">Genomes.txt") || die "Could not open Genomes.txt\n";
	print TMP "$sequence\n";
	close (TMP);

	foreach my $matrix_file (@matrix_files){
		my $patser=qx(~/Desktop/shraddha_project/Consensus-2008/patser-v3e.2008/patser-v3e -m $matrix_file -f Genomes.txt -A $gc_content -c -M -1 -li);
		my $motif_size;
			my @file_name=split(/\\|\./,$matrix_file);
		my @Promoters=split(/\n/,$patser);
		my @lineToProcess;		
		
		foreach my $line(@Promoters){
			my $strand="+";
			if ($line=~/width[\s|\w|\W](\d+)/){
				$motif_size=$1;
				}
			elsif( $line=~/\/\w+\/\w+\//){
				@lineToProcess=split (/\s{1,}/, $line);
				my $start=$lineToProcess[2];
				my $p_value=$lineToProcess[$#lineToProcess];
				
				if ($lineToProcess[2]=~/C/){
					$strand="-";
					$start=substr($lineToProcess[2],0,length($lineToProcess[2])-1);
					}
					
				$Temp_promoter{$index}={
						"ID"	=> $index,
						"Matrix"=>$file_name[$#file_name-1],
						"Strand" => $strand,
						"Start" => $start,
						"Stop"	=> $start-$motif_size,
						"Pvalue" => $p_value
						};
						$index++;
					}
				}
			}
			foreach my $Key(sort { $Temp_promoter{$a}->{"Start"} <=> $Temp_promoter{$b}->{"Start"} } keys(%Temp_promoter)){
				if ($Temp_promoter{$Key}{"Strand"} eq "-"){
					${$promoter_ref}{$Temp_promoter{$Key}{"Strand"}}{$index_m}={
						"ID"	=> $Temp_promoter{$Key}{"ID"},
						"Matrix"=>$Temp_promoter{$Key}{"Matrix"},
						"Strand" => $Temp_promoter{$Key}{"Strand"},
						"Start" => $Temp_promoter{$Key}{"Start"},
						"Stop"	=> $Temp_promoter{$Key}{"Stop"},
						"Pvalue" => $Temp_promoter{$Key}{"Pvalue"}
						};	
					push @{$coordinates{$Temp_promoter{$Key}{"Strand"}}}, $Temp_promoter{$Key}{"Start"};
					$index_m++;
					}
				else {
					${$promoter_ref}{$Temp_promoter{$Key}{"Strand"}}{$index_p}={
						"ID"	=> $Temp_promoter{$Key}{"ID"},
						"Matrix"=>$Temp_promoter{$Key}{"Matrix"},
						"Strand" => $Temp_promoter{$Key}{"Strand"},
						"Start" => $Temp_promoter{$Key}{"Start"},
						"Stop"	=> $Temp_promoter{$Key}{"Stop"},
						"Pvalue" => $Temp_promoter{$Key}{"Pvalue"}
						};
					push @{$coordinates{$Temp_promoter{$Key}{"Strand"}}}, $Temp_promoter{$Key}{"Start"};
					$index_p++;
					}
				}
	print "Patser-v3e search finished return code: $?\n";
	return %coordinates;
}		

sub Infernal {
#####################################################################################################
#																									#
#											Infernal												#
#																									#
#		E. P. Nawrocki, D. L. Kolbe, and S. R. Eddy, Infernal 1.0: Inference of RNA alignments,		#
#								Bioinformatics 25:1335-1337 (2009)									#		
#																									#
#####################################################################################################
	my $sequence=${$_[0]};
	my $cm_file= ${$_[1]};
	my $Tbox_ref=$_[2];
	my %coordinates;
	my ($index_p, $index_m, $start, $stop, $ID)=(0,0,0,0,1);
	my $infernal=qx(cmsearch $cm_file $sequence);
	
	my @tbox=split(/\n/,$infernal);	
	my @lineToProcess;
	
	for (my $i=0; $i<=$#tbox; $i++){
		my $line=$tbox[$i];
		if ($line=~/Target/gi){
			my @Scores=split(/\s{1,}/, $tbox[$i+1]);
			if ($Scores[3]>=5){
				my $strand='';
				@lineToProcess=split (/\s{1,}/, $line);
			
				if ($lineToProcess[8] >= $lineToProcess[10]){
					$strand = "-";
					$start=$lineToProcess[10];
					$stop=$lineToProcess[8];
					}
				elsif($lineToProcess[8] <= $lineToProcess[10]){
					$start=$lineToProcess[8];
					$stop=$lineToProcess[10];
					$strand = "+";
					}
				else {
					die "Could not determine directionality in Tbox parsing \n";
					}
			
				if ($strand eq "-"){
					${$Tbox_ref}{$strand}{$index_m}= {
						"ID"	=> "TBox".$ID,
						"Start" => $start,
						"Stop"	=> $stop,
						};
					push @{$coordinates{$strand}}, $start;
					$index_m++;
					}		
				else {
					${$Tbox_ref}{$strand}{$index_p}= {
						"ID"	=> "TBox".$ID,
						"Start" => $start,
						"Stop"	=> $stop,
						};
					push @{$coordinates{$strand}}, $start;
					$index_p++;
					
				}
			$ID++;
			}
		}
	}
	#print Dumper(%results)."\n";
	print "Infernal T-box search finished return: $?\n";
	return %coordinates;
}

sub rev_comp{
	my $seq=shift;
	my $seq_rev=reverse($seq);
	$seq_rev=~tr/ATGC/TACG/;
	return $seq_rev;
	}

sub start_end_strand{

	my $start=shift;
	my $end=shift;
	my $left=shift;
	my $right=shift;
	
	if ($start>$end){
		$$left=$end;
		$$right=$start;
		}
	else{
		$$right=$end;
		$$left=$start;
		}
	}

sub between{

	my $left=shift;
	my $right=shift;
	my $valuetotest=shift;
	
	if($valuetotest==$left || $valuetotest==$right){
		return 1;
		}
	
	elsif($left<=$right){
		if($valuetotest>=$left && $valuetotest<=$right){
			return 1;
			}
		else{ return 0;	}
		}
		
	elsif($left>=$right){	
		if($valuetotest<=$left && $valuetotest>=$right){
			return 1;
			}
		else{ return 0;	}
		}
	else {	print "Numbers are not a valid range\n$left\t$right\n";
			die;
		}
	}
	
sub gff3Parser{
	my @InputFile=@{$_[0]};
	my $FeatureContainer=$_[1];
	my $gffFields=$_[2];
	chomp (@InputFile);
	
	foreach (@InputFile){
		unless (/^#/){
			my @Fields=split(/\t/,$_);
			my @AttributeTags=split(/;/,pop(@Fields));
			my %AttributeField=();
			my $ID;
			
			foreach my $AttributePair(@AttributeTags){
				my ($Tag,$Value)=split(/=/,$AttributePair,2);
				$Tag=~s/'|"//g;
				$Value=~s/'|"//g;
				$AttributeField{$Tag}=$Value;
				}
			
			if ($Fields[2] eq "gene" || $Fields[1] eq "fuzznuc"){
				if($AttributeField{'ID'}){
					$ID=$AttributeField{'ID'};
					}
				else {
					$ID=$AttributeField{'locus_tag'};
					}
					
				for (1..$#Fields){	
					${$FeatureContainer}{$Fields[0]}{$ID}{${$gffFields}[$_]}=$Fields[$_];
					}
				while (my ($k, $v)=each %AttributeField){
					${$FeatureContainer}{$Fields[0]}{$ID}{$k}=$v;
					}
				}
			}
		}
	}
	
sub smallest_difference {
	my $anchor_value=shift;
	my $value1=shift;
	my $value2=shift;
	
	if (abs($anchor_value-$value1)<abs($anchor_value-$value2)){
		return abs($anchor_value-$value1);
		}
	else{
		return abs($anchor_value-$value2);
		}
	}
	
sub closest{

	my $test_value=shift;
	my $offset=shift;
	my $distance=abs($test_value-$_[0]);
	my $index;
	for(my $i=1; $i<=$#_; $i++){
		my $test_distance=abs($test_value-$_[$i]);
		if ($distance>$test_distance){
			$distance=$test_distance;
			$index=$i;
			}
		}
	my $index_return=$offset+$index;
	return $index_return, $distance;
}		

sub html_head{
	my $file_out=${$_[0]};
	my $file=$_[1];
	my $algoName=TRAPalgorithm3::algorithmName();
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $timestamp=($year+1900)."-".++$mon."-$mday $hour:$min:$sec";
	my $repeatRegEx=TRAPalgorithm3::regexDump();

print $file_out <<HTML1;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html>
	<head>
		<script src="sorttable.js"></script>
		<link rel="stylesheet" type="text/css" href="mapper.css" />
	</head>
	<body>
	<table class="header">
		<tr>
			<td class="cat">Execution Time:</td>
			<td class="param">$timestamp</td>
		</tr>
		<tr>
			<td class="cat">Input file:</td>
			<td class="param">$ARGV[1]</td>
		</tr>
		<tr>
			<td class="cat">Output file:</td>
			<td class="param">$file</td>
		</tr>
		<tr>
			<td class="cat">Repeat RegEx:</td>
			<td class="param">$repeatRegEx</td>
		</tr>
		<tr>
			<td class="cat">Algorithm:</td>
			<td class="param">$algoName</td>
		</tr>
	</table>
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



}

sub html_end {
	my $file_out=${$_[0]};
	print $file_out <<HTML2;
	
		</table>
	</body>
	</html>	
HTML2
	
}

sub calculate_GC {
  my (@sequence) = @_;
  my ($n, $GC, $AT) = (0, 0, 0); # initialise counters
  my $seqread  ;
  
  foreach my $seqread(@sequence) {  
    $GC++ if $seqread =~ /[GCgc]/ ;
    $AT++ if $seqread =~ /[ATat]/ ;
  }
  $n = $GC + $AT ;
  my $GC_content = $GC / $n ;
  return $GC_content ;
}

sub show_help{

	print <<HELP;

This script is designed to search for "regular expression-like" patterns in 
nucleic acid sequence (specifically genomes, though with some tweaking cDNAs, 
ESTs and Contigs could be supported as well) using the EMBOSS program fuzznuc.  
These patterns are then annontated against the genomic features of the given 
organism, as well as for promoters and rho-independent terminators.  The 
fuzznuc results are also scored (currently for TRAP binding, but any scoring 
module could be used). Outputs include the final fuzznuc results in GFF3 format, 
seqeunces for each fuzznuc  result in fasta, HTML for scoring and binding site 
mark-up, and tab-delimited for annontations (all, and intergenic only).

This script does not require any extra modules besides those that are 
distributed with the script (TRAPAlgorithm3, BinarySearch, and Hit) or
are part of the Core Perl distribution (GetOpt::Long and File::IO).

To run this script, the EMBOSS suite needs to be installed (6.0 or later) 
as well as TransTermHP see (Kingsford et al Genome Biology 2007). Both of which
are freely available.

Before running this script, it is advisable, if possible, to first run
TRAP_source_FTP.pl to acquire the necessary files for analysis.  By doing this
you will ensure that the proper files are present and formatted correctly for use by
all programs (specifically TransTermHP has strict file requirements).

Usage: 
	
	
HELP

}