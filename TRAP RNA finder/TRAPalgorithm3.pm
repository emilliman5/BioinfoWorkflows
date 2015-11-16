package TRAPalgorithm3;

use strict;

#written by C. Jackson 2/2011 University at Buffalo

BEGIN {
	#Define scores for the repeat types in the scoring algorithm
	use constant {
		GSCORE => 12,
		CSCORE => 2,
		ASCORE => 3,
		TSCORE => 12,
		MINREPEAT => 6,
		MAXREPEAT => 11,
		MINSCORE => -200
	};

	#Define scores for the spacer size in the scoring algorithm
	use constant {
		ZEROSPACE => -5,
		ONESPACE => 2,
		TWOSPACE => 4,
		THREESPACE => 2,
		GAPOPEN => -5,
		GAPEXTEND => -1, 		#Score for spacers greater then 3N, per nucleotide over 3
		MAXSKIP => 5,			#Maximum number of spacers greater than 3N before throwing in the towel
		MAXSPACERSIZE => 20		#Defines the largest possible spacer
	};
	
	use constant REGEXP => "[GATC]AG";
	
	use constant NORMALIZE => 1;	#Set as 1 to normalize all scores by repeat count
									#I don't know why I wrote this or how it could help anything
	
	use constant {
		WEIGHTPOSITION => 0,		#Boolean: Set to 1 to activate weighting
		WEIGHTLENGTH => 2,			#Number of repeats on the outside to apply EXTERNALMOD to
		EXTERNALMOD => 1.5,
		INTERNALMOD => 1.0
	};
	
	use constant VERSION => "1.4";
}

#########################################################################
#																		#
#	algorithmName returns the constants used in the scoring algorithm	#
#		-Useful for keeping track of what fiddling's been done			#
#																		#
#########################################################################

sub algorithmName {
	return "TRAP @{[VERSION]} ([@{[MINSCORE]},[@{[GSCORE]},@{[CSCORE]},@{[ASCORE]},@{[TSCORE]},@{[MINREPEAT]},@{[MAXREPEAT]}],[@{[ZEROSPACE]},@{[ONESPACE]},@{[TWOSPACE]},@{[THREESPACE]},@{[GAPOPEN]},@{[GAPEXTEND]},@{[MAXSKIP]}, @{[MAXSPACERSIZE]}],@{[NORMALIZE]})";
}

#########################################################################
#																		#
#	regexDump returns the regexp used in finding repeats 				#
#		-Useful for keeping track of what fiddling's been done			#
#																		#
#########################################################################

sub regexDump {
	return "@{[REGEXP]}";
}

#########################################################################
#																		#
#	algorithmScore														#
#		Arguments:														#
#			string	sequence											#
#		Returns:														#
#			boolean	successFlag											#
#			string	HTMLsequence										#
#			int		score												#
#			int		repeats												#
#			int		start (of scored area)								#
#			int		end (of scored area)								#
#########################################################################

sub algorithmScore { #returns SuccessFlag, HTMLSequence, Score, Repeats, Start, End
	my $aSeq = $_[0];
	my $seqLen=length($aSeq);
	
	my @repeatLoc=();
	my @repeatN=();
	
	while($aSeq =~ /[GT]AG/gi) {
		if(substr($aSeq,$-[0]+2,3) eq "GAG" && (($-[0]-$repeatLoc[-1]-3)<2 || substr($aSeq,$-[0],1) eq "T")) {push(@repeatLoc,($-[0]+2)); next;}
		if(($-[0]-$repeatLoc[-1])<3 && $#repeatLoc>-1) {next;}
		push (@repeatLoc,$-[0]);	
	}
		
	SCAN: while($aSeq =~ /[AC]AG/gi) {
		my $currentLoc=$-[0];
		foreach(@repeatLoc) {
			if(abs($currentLoc-$_)<4) {next SCAN;}
		}
		push(@repeatLoc,$currentLoc);
	}
	
	@repeatLoc=sort {$a <=> $b} @repeatLoc;
	
	foreach(@repeatLoc) {
		push(@repeatN,substr($aSeq,$_,1));
	}
	
	my $repeatCount=$#repeatLoc+1;
	if(&MINREPEAT()>$repeatCount) {return (0,"",0,0,0,0);} #Kick back fail flag if there's not enough repeats
	
	my $score=&MINSCORE();
	my @highScore=(); #(Start, End)
	
	for(my $start=0;$start<=($repeatCount-&MINREPEAT()-1);$start++) { #This walks through the possible starting positions for the scoring frame
		my @maxEnd=($repeatCount-1,$start+&MAXREPEAT()-1);
		@maxEnd=sort {$a <=> $b} @maxEnd;
		my $lastLoc=0;
		FRAME: for (my $end=($start+&MINREPEAT()-1);$end<=$maxEnd[0];$end++) { #This walks through the possible ending positions for this scoring frame anchored to a starting position
			my $currentScore=0;
			my $currentRepeatCount=$end-$start;
			my @currentSliceL=@repeatLoc[$start..$end];
			my @currentSliceN=@repeatN[$start..$end];
			my $last=shift(@currentSliceL);
			my $currentSpacerCount=0;
			my $currentRepeatPosition=0;
			foreach(@currentSliceL) {
				my $spacer=abs($_-$last-3);
				if($spacer>&MAXSPACERSIZE()) {next FRAME;}
				if(0>$spacer) 	{die "REGEX ERROR";}
				elsif(0==$spacer)	{$currentScore+=ZEROSPACE; $last=$_;}
				elsif(1==$spacer)	{$currentScore+=ONESPACE; $last=$_;}
				elsif(2==$spacer)	{$currentScore+=TWOSPACE; $last=$_;}
				elsif(3==$spacer)	{$currentScore+=THREESPACE; $last=$_;}
				elsif(3<$spacer) 	{$currentScore+=(&GAPOPEN()+&GAPEXTEND()*($spacer-3)); $currentSpacerCount++;  $last=$_;}
			}
			if(&MAXSKIP()<$currentSpacerCount) {next FRAME;}
			foreach(@currentSliceN) {
				my $thisRepeatScore=0;
				if($_ eq "A")	{$thisRepeatScore=ASCORE;}
				if($_ eq "T")	{$thisRepeatScore=TSCORE;}
				if($_ eq "G")	{$thisRepeatScore=GSCORE;}
				if($_ eq "C")	{$thisRepeatScore=CSCORE;}
				if(&WEIGHTPOSITION()) {
					my $checkflag=0;
					while(1) {
						if(&WEIGHTLENGTH()>$currentRepeatPosition) {$thisRepeatScore*=&EXTERNALMOD(); last;}
						if(($currentRepeatCount-&WEIGHTLENGTH())=>$currentRepeatPosition) {$thisRepeatScore*=&INTERNALMOD(); last;}
						if(($currentRepeatCount-&WEIGHTLENGTH())<$currentRepeatPosition) {$thisRepeatScore*=&EXTERNALMOD(); last;}
					}
					$currentRepeatPosition++;
				}
				$currentScore+=$thisRepeatScore;
			}
			if(&NORMALIZE()) {$currentScore=$currentScore/($end-$start+1);}
			if ($currentScore>$score) { #Throw out the score unless it's higher then 
				$score = $currentScore;
				@highScore=($start,$end);
			}
		}	
	}
	
	if(&MINSCORE()<$score && &MINREPEAT()<=($highScore[1]-$highScore[0]+1)) {
		my @htmlMarker=@repeatLoc[$highScore[0]..$highScore[1]];
		my $last=0;
		my $retSeq="";
		foreach(@htmlMarker) {
			if(($_-$last)<0) {print "Overlap Frame Error: $_ $last\n"; next;}
			$retSeq.=substr($aSeq,$last,($_-$last))."<span class=\"repeat\">".substr($aSeq,$_,3)."</span>";
			$last=$_+3;
			
		}
		$retSeq.=substr($aSeq,$last);
		return (1,$retSeq,$score,($highScore[1]-$highScore[0]+1),$repeatLoc[$highScore[0]],$repeatLoc[$highScore[1]]);
	}
	
	return (0,"",0,0,0,0);
}

1;