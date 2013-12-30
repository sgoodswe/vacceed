#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $algorithm = $ARGV [0];
my $evidence_dir = $ARGV [1] . "/";


## hard coding
my $input_file = $evidence_dir . $algorithm . "_predictions.txt";
my $output_file = $evidence_dir . $algorithm . "_ml";

#Global variables
my $evidence = 0;

#create hashes
my %h_yes;
my %h_no;
my %h_count;


#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;

#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	if ($line =~ m/^(\d+)/) { 
	
		my ($no, $id, $candidate, $class) = split ('\s+',$line);

		if ($class eq "YES") {
			$h_yes {$id} = $h_yes {$id} + 1;
		} else {

			$h_yes {$id} = $h_yes {$id} + 0;
		}

		$h_count {$id} = $h_count {$id} + 1;
	}
			
} #end of LINE:while


foreach my $id ( sort keys %h_count) {
	$evidence = $h_yes {$id} / $h_count {$id};
	print OUTPUT_FILE $id . "," . $evidence . "\n";
}



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

## multiple parameters
sub  add_to_cluster ($)
{
	my ($query_id, $subject_id,$c) = (@_);

}
