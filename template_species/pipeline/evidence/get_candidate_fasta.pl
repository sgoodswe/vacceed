#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $input_file = $ARGV [0];
my $prefix = $ARGV [1];
my $proteome_dir = $ARGV [2] . "/";
my $ml_threshold = $ARGV [3];
my $existence_threshold = $ARGV [4];


## hard coding
my $candidate_file = $proteome_dir . "vaccine_candidates";
my $output_file = $proteome_dir . "vaccine_candidates.fasta";


#Global variables
my $required_id = 0;

#create hashes
my %h_final_score; 
my %h_existence_score;

#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open INPUT_FILE,$input_file ;
open CANDIDATE_FILE,$candidate_file ;
open OUTPUT_FILE,'>' . $output_file;

#LOOP through Candidate file 
while (defined (my $line = <CANDIDATE_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	my @scores = split (',',$line);
	my $id = $scores [0];

	$h_final_score {$id} = $scores [scalar @scores -1];
	$h_existence_score {$id} = $scores [scalar @scores -2];
			
} #end of LINE:while


#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);

	#Check for # 
	if ($line =~ m/^>$prefix\|(\w+)\|/) {
	
		my $id = $1;
		
		if (($h_final_score {$id} >= $ml_threshold) && ($h_existence_score {$id} >=$existence_threshold)) {
			print OUTPUT_FILE $line . "\n";
			$required_id = 1;
		} else {
			$required_id = 0;
		}
	} elsif ($required_id == 1) {
		print OUTPUT_FILE $line . "\n";
	}
	
			
} #end of LINE:while


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


