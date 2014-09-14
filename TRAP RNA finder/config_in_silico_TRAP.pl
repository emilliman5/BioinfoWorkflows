#!/usr/bin/perl -w

use strict;

print "Welcome to the in silico TRAP binding site indentification configuration script...\n
This script will check for the necessary dependnecies and set relevant paths to said dependencies as well as set some analysis parameters\n";

print "First we will check for EMBOSS' fuzznuc...\n";

system('./fuzznuc') == 0 or die "Fuzznuc failed: $?\n";

print "Next we will check for TransTermHP...\n";

eval {

	
			}
			
			
print "Next we will test for the appropriate Perl Modules...\n";

eval	{


			}
			
			
