#!/usr/bin/perl
use strict;
#use warnings;

my $input_file = $ARGV [0];
my $proteome_dir = $ARGV [1] . '/';
my $output_dir = $ARGV [2] . '/';
my $resource = $ARGV [3];
my $prot_id_prefix =  $ARGV [4];


my $output_file = $output_dir . $resource . "_" . $input_file;
my $temp_file = $proteome_dir . "tmp_" . $input_file;
my $input_file = $proteome_dir . $input_file;


#Global variables
my $id;
my $invalid_detected = 0;

#create hashes
my %h_seq;


#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open INPUT_FILE,$input_file ;

#Open temporary proteome.fasta file to handle long gene identifiers
open TEMP_FILE, '>' . $temp_file ;

#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);
	
	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	if ($line =~ m/^>$prot_id_prefix\|(\w+)/) { 

		print TEMP_FILE ">" . $prot_id_prefix . "|" . $1 . "\n";
		$id = $line;
	} else {
		print TEMP_FILE $line . "\n";
		
		if ($line !~ /\A[ARNDBCEQZGHILKMFPSTWYV]+\z/ix) {
			$invalid_detected = 1;
			$h_seq {$id} =  $h_seq {$id} . $line;
		}
	}
			
} #end of LINE:while



if ($invalid_detected == 1) {

	open OUTPUT_FILE, '>' . $output_file ;
	print OUTPUT_FILE "The following protein sequences contains characters other than [ARNDBCEQZGHILKMFPSTWYV]:\n\n";
	
	foreach my $id (sort keys %h_seq) {
		print OUTPUT_FILE $id . "\n";
		print OUTPUT_FILE $h_seq {$id} . "\n";
	}

	print "\n!! Program terminated: Invalid characters detected in $input_file \n\n";
	print "Please check the file $output_file for more details\n";
	exit 1;
	
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
