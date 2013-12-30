#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $input_file = $ARGV [0];
my $evidence_dir =  $ARGV [1] . "/";
my $resource = $ARGV [2];
my $start_prefix = $ARGV [3];
my $split_char = $ARGV [4];
my $id_info = $ARGV [5];
my $cols_required = $ARGV [6];
my $evd_headers = $ARGV [7];
my $ignore_lines = $ARGV [8];
my $extract_value = $ARGV [9];

## DEBUGGING
#print  $input_file . "\n";
#print $evidence_dir. "\n";
#print $resource . "\n";
#print $start_prefix . "\n";
#print $split_char . "\n";
#print $id_info . "\n";
#print $cols_required . "\n";
#print $evd_headers . "\n";
#print $ignore_lines . "\n";
#print $extract_value . "\n";


## hard coding
my $output_file = $evidence_dir . $resource . "_evd";

#split the ID info. Extract 1: Which column in the input line contains the identifier, the position of the id in the identifier, the separating character betweens parts of the identifier.
my ($id_col,$id_pos,$id_sep) = split (',',$id_info);

#Global variables

#create hashes

#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;


print OUTPUT_FILE "ID," . $evd_headers . "\n";   

#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	#reset variables
	my $col_count = 0;
	my $id;
	my @evidence = ();

	#Check for #  or blank line
	if (($line =~ m/^#/) || ($line =~ m/^\s*$/)) { next;}
	
	$line = &trim ($line);

	#Check if we need to ignore the line
	if (&ignore_line ($line) ) { next;}

	
	#check if we are searching for a specific but consistent start_prefix
	if (!&check_for_empty_string ($start_prefix)) {
		if ($line !~ m/^$start_prefix/) { next;}
	}
	
	#split the line input
	my @line_input = split ($split_char,$line);

	#Loop for each column in the line input
	foreach my $col (@line_input) {
	
		$col_count++;
		
		if ($col_count == $id_col) {
		
			$id = &get_ID ($col);
		
			if ($id eq "0") {
				print "ERROR: ID not found on line $line\n";
				exit 1;
			}
		}
		if (&column_check ($col_count)) {

			#check if we are extracting a  value within the column text itself e.g.  ExpAA=43.71
			if (!&check_for_empty_string ($extract_value)) {

				$col =~ m/$extract_value(.+)/;
				$col = $1;
			}

			push (@evidence,$col);
		}
	}
	
	#Print the required evidence  to the output file;
	print OUTPUT_FILE $id;

	foreach my $col (@evidence) {
		print OUTPUT_FILE "," . $col;
	}

	print OUTPUT_FILE "\n";
		
	
			
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

## Get the ID
sub  get_ID ($)
{
	my $id_line = shift;
	my $split_count = 0;
	
	my @id_input = split ($id_sep,$id_line);
		
	foreach my $col (@id_input) {
		
		$split_count++;

		if ($split_count == $id_pos) {
			return $col;
		}
	}

	return 0;
}

## Get the ID
sub  column_check ($)
{
	my $col_to_check = shift;
		
	my @columns_input = split (',',$cols_required);
		
	foreach my $col (@columns_input) {
	
		if ($col == $col_to_check) {
			return 1;
		}
	}

	return 0;
}


## Get the ID
sub  ignore_line ($)
{
	my $line = shift;
		
	my @ignore_input = split (',',$ignore_lines);
		
	foreach my $ignore (@ignore_input) {
	
		if ($line =~ m/^$ignore/) {
			return 1;
		}
	}

	return 0;
}
