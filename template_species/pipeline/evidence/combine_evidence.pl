#!/usr/bin/perl
use strict;
#use warnings;

my $output_file = $ARGV [0];
my $file_suffix = $ARGV [1];
my $evidence_dir = $ARGV [2];

#hard coding;
my $predictor_file = $evidence_dir . "/" . "ml_predictors";


## debugging
#print $output_file . "\n";
#print $file_suffix . "\n"; 
#print $evidence_dir  . "\n"; 

#Global variables

#create hashes
my %h_evidence_profiles;
my $header;

#Open debug file
#open DEBUG_FILE,'>debug.txt';


#Get files in working directory
opendir(DIR, $evidence_dir );
my @FILES= readdir(DIR); 

#File handling
open OUTPUT_FILE,'>' . $output_file;
open PREDICTOR_FILE,'>' . $predictor_file;

foreach my $file (@FILES) {

	if ($file !~ m/$file_suffix/) { next;}

	open INPUT_FILE, $evidence_dir . "/" . $file ;
	
	print $file . "\n";

	#LOOP 
	while (defined (my $line = <INPUT_FILE>)) {
	
		#Check for # 
		if ($line =~ m/^#/) { next;}
		
		$line = &trim ($line);
		
		my ($id,$rest) = split (',',$line);
	
		
		#get the rest of the line
		$id = quotemeta $id;
		$line =~ m/$id,(.+)/;
		
		my $rest_of_line = $1;
		
		#combine the header
		if ($id eq "ID") {
			if (&check_for_empty_string ($header)) {
				$header = "ID," . $rest_of_line;
			} else {
				$header = $header . "," . $rest_of_line;
			}

			my @headers = split (',',$rest_of_line);
			foreach (@headers) {

				print PREDICTOR_FILE $_ . "\n";
			}

		}elsif (&check_for_empty_string ($h_evidence_profiles {$id})) {
				$h_evidence_profiles {$id} = $rest_of_line;
		} else {
				$h_evidence_profiles {$id} = $h_evidence_profiles {$id} . "," . $rest_of_line;
		}
				
	} #end of LINE:while
	
	close INPUT_FILE;

}

print OUTPUT_FILE $header . ",Candidate" . "\n";

foreach my $id (sort keys %h_evidence_profiles) {

	print OUTPUT_FILE $id . "," . $h_evidence_profiles {$id} . "\n";

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
