#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $input_file = $ARGV [0];
my $split_by = $ARGV [1];
my $output_dir = $ARGV [2] . "/";


## hard coding
my $output_file = $output_dir . "file_";


#Global variables
my $no_of_proteins = 0;
my $protein_count = 0;
my $protein_limit = 0;
my $file_count = 0;
my $id;

#create hashes
my %h_sequences;


#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;


#Count the number of proteins in fasta file
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);
	
	if ($line =~ m/^>/) {
		$no_of_proteins++;
		$id = $line;
	} else {
		$h_sequences {$id} = $h_sequences {$id} . $line;
	}
		
} #end of LINE:while

print "number of proteins "  . $no_of_proteins . "\n";

if ($split_by > $no_of_proteins) {
	$protein_limit = $no_of_proteins; 
} elsif ($split_by != 0) {
	$protein_limit = $no_of_proteins /$split_by;
} else {
	$protein_limit = $no_of_proteins /10;
}

my $proteins_per_file = sprintf ("%d",$protein_limit);

print "Approximate number of proteins per file "  . $proteins_per_file . "\n"; 

if (($no_of_proteins !=0) && ($protein_limit != 0)) {

	$file_count++;
	open OUTPUT_FILE,'>' . $output_file . $file_count . ".fasta";

	#LOOP 
	foreach my $seq_id (sort keys %h_sequences) {
			
		$protein_count++;

		if ($protein_count > $protein_limit) {

			$file_count++;
			open OUTPUT_FILE,'>' . $output_file . $file_count . ".fasta";
			print OUTPUT_FILE $seq_id . "\n";
			print OUTPUT_FILE $h_sequences {$seq_id} . "\n";
			$protein_count = 1;
	
		} else {
			print OUTPUT_FILE $seq_id . "\n";
			print OUTPUT_FILE $h_sequences {$seq_id} . "\n";
		}

		
		#change protein limit for last file
		if ($file_count == $split_by) {
			$protein_limit = $no_of_proteins - ($proteins_per_file * ($split_by - 1));
			
		}

	
	} #end of LINE:while
}

print "number of files to run " . $file_count . "\n";
exit $file_count;

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
