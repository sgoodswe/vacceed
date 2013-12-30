#!/usr/bin/perl
use strict;
#use warnings;

#get the command line arguments
my $evidence_dir = $ARGV [0] . "/";
my $proteome_dir = $ARGV [1] . "/";


## hard coding
my $output_file = $proteome_dir . "vaccine_candidates";
my $proteome_info = $proteome_dir . "proteome_info.txt";

#Global variables
my $header;

#create hashes
my %h_scores;
my %h_count;
my %h_total;
my %h_avg_ml_score;
my %h_proteome_score;


#Open debug file
#open DEBUG_FILE,'>debug.txt';

## Read the proteome file if it exists
open PROTEOME_FILE,$proteome_info ;

while (defined (my $line = <PROTEOME_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	my @scores = split (',',$line);
	my $id = $scores [0];
	$h_proteome_score {$id} = $scores [scalar @scores -1];
}


#Get files in evidence output directory
opendir(DIR, $evidence_dir);
my @FILES= readdir(DIR); 

#File handling

open OUTPUT_FILE,'>' . $output_file;

foreach my $file (@FILES) {

	if ($file !~ m/(\w+)_ml/) { next;}
	
	if (&check_for_empty_string ($header)) {
		$header =  "#ID," . $1;
	} else {
		$header = $header . "," . $1;
	}
	
	open INPUT_FILE, $evidence_dir . $file ;

	#LOOP 
	while (defined (my $line = <INPUT_FILE>)) {

		#Check for # 
		if ($line =~ m/^#/) { next;}
		
		$line = &trim ($line);
		
		my ($id,$score) = split (',',$line);
		
		$h_total {$id} = $h_total {$id} + $score;
		$h_count {$id} = $h_count {$id} + 1;

		$score = sprintf ("%.3f",$score);
		
		if (&check_for_empty_string ($h_scores {$id})) {
			$h_scores {$id} = $score;
		} else {
			$h_scores {$id} = $h_scores {$id} . "," . $score;
		}
				
	} #end of LINE:while

}

#Determine average scores
foreach my $id (sort keys %h_total) {

	$h_avg_ml_score {$id} = $h_total {$id}/$h_count {$id};
}

#Print out headers
print OUTPUT_FILE $header . ",existence_score,average_ML_score\n";

#print out scores in descending order for final score
foreach my $id (sort { $h_avg_ml_score {$b} <=> $h_avg_ml_score {$a}} keys %h_avg_ml_score) {

	my $avg_ml_score = sprintf ("%.3f",$h_avg_ml_score {$id});
	my $proteome_score  = sprintf ("%.3f",$h_proteome_score {$id});
	
	
	print OUTPUT_FILE $id . "," . $h_scores {$id} . ","  . $proteome_score . "," .  $avg_ml_score . "\n";

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
