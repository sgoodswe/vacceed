#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $proteome_dir = $ARGV [0] . "/";
my $prot_dir =  $ARGV [1] . "/";
my $prot_file =  $ARGV [2];
my $new_prot_file =  $ARGV [3];

#hard coding
my $output_file = $proteome_dir . "proteome.fasta";


my $protein_file = $prot_dir . $prot_file;
my $new_protein_file = $proteome_dir .$new_prot_file;



#Global variables

#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open PROT_FILE, $protein_file;
open NEW_PROT_FILE, $new_protein_file;

open OUTPUT_FILE,'>' . $output_file;


#LOOP  through protein file
while (defined (my $line = <PROT_FILE>)) {

	print  OUTPUT_FILE $line;
			
}

#LOOP  through new protein file
while (defined (my $line = <NEW_PROT_FILE>)) {

	print  OUTPUT_FILE $line;
			
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


