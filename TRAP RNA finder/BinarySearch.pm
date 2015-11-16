package BinarySearch;

use strict;

BEGIN {

####################################################################################################
########																					########
########					Binary Search from adapted from PerlMonks						########
########																					########
####################################################################################################

sub cmpFunc {
	my ($index, $arrayRef, $target) = @_;
	my $item = $$arrayRef[$index];
	
	return $item <=> $target;		# returns -1 is $item < $target; 1 if $item > $target and 0 if $item == $target
}


sub BinSearch {
	my $target = $_[0];
	my @array = @{$_[1]};
	
	my $posmin = 0;
	my $posmax = $#array;
	
	return -0.5 if &cmpFunc (0, \@array, $target) > 0;					#checks if target is below minimum value in array
	return $#array + 0.5 if &cmpFunc ($#array, \@array, $target) < 0;	#checks if target is above maximum value in array
	
	while (1)
	  {
	  my $mid = int (($posmin + $posmax) / 2);
	  my $result = &cmpFunc ($mid, \@array, $target);
	  
	  if ($result < 0)
		{
		$posmin = $posmax, next if $mid == $posmin && $posmax != $posmin;
		return $mid + 0.5 if $mid == $posmin;
		$posmin = $mid;
		}
	  elsif ($result > 0)
		{
		$posmax = $posmin, next if $mid == $posmax && $posmax != $posmin;
		return $mid - 0.5 if $mid == $posmax;
		$posmax = $mid;
		}
	  else
		{
		return $mid;
		}
	  }
	}
####################################################################################################
########																					########
########						End of Binary Search from PerlMonks							########
########																					########
####################################################################################################
 
 }
1;