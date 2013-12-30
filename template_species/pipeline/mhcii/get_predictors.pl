#!/usr/bin/perl
use strict;
#use warnings;

#get the command line arguments
my $input_file = $ARGV [0];
my $out_dir = $ARGV [1] . "/";

## hard coding
my $output_file = $out_dir . "ml_predictors.R";


#Global variables
my $first = 0;

#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;

print OUTPUT_FILE "input <- c(";

#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);
		
	my ($allele) = split (',',$line);
	
	#substitute -
	$allele =~ s/-/./g;
	
	#substitute :
	$allele =~ s/:/./g;
	
	#substitute *
	$allele =~ s/\*/./g;

	#substitute /
	$allele =~ s/\//./g;
	
	if ($first == 0) {
	
		print OUTPUT_FILE  '"' . $allele . '"';
		$first = 1;
	} else {
	
		print OUTPUT_FILE  ',"' . $allele . '"';
	
	}

			
} #end of LINE:while

print OUTPUT_FILE ")";

print OUTPUT_FILE "\n";

#####################SUB ROUTINES ################################

# sub-routine to remove  whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


#Check for empty string
sub  check_for_empty_string ($)
{
	my $string = shift;
	if ($string =~ m/^\s*$/) {
		return 1;
	} else {
		return 0;
	}
}
