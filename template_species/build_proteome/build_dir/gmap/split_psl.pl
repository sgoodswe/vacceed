#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAM "split_psl.pl"

#get the command line arguments
my $est_file = $ARGV [0];
my $output_dir = $ARGV [1] . "/";
	

## hard coding
my $extention = "_est.psl";
my $input_file = $output_dir . $est_file;


#Global variables

#create hashes
my %h_chrs;

#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open INPUT_FILE,$input_file ;


## PSL format
#matches - Number of bases that match that aren't repeats
#misMatches - Number of bases that don't match
#repMatches - Number of bases that match but are part of repeats
#nCount - Number of 'N' bases
#qNumInsert - Number of inserts in query
#qBaseInsert - Number of bases inserted in query
#tNumInsert - Number of inserts in target
#tBaseInsert - Number of bases inserted in target
#strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
#qName - Query sequence name
#qSize - Query sequence size
#qStart - Alignment start position in query
#qEnd - Alignment end position in query
#tName - Target sequence name
#tSize - Target sequence size
#tStart - Alignment start position in target
#tEnd - Alignment end position in target
#blockCount - Number of blocks in the alignment (a block contains no gaps)
#blockSizes - Comma-separated list of sizes of each block
#qStarts - Comma-separated list of starting positions of each block in query
#tStarts - Comma-separated list of starting positions of each block in target' 


#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}

	$line = &trim ($line);
	
	my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,
	$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts) = split ('\s+',$line);
	
	my $chr;
	
	#get chromosome number from target name
	if ($tName =~ m/^chr(\w+)/){	
		 $chr = $1;
	} elsif ($tName =~ m/(\w+)/) {
		$chr = $1;
	} 
	
	
	#check if we have opened the file
	if (&check_for_empty_string ($h_chrs {$chr})) {
	
		$h_chrs {$chr} = $chr;
		open OUTPUT_FILE,'>' . $output_dir . "chr" . $chr . $extention;
		print OUTPUT_FILE $line . "\n";
		
	} else {
	
		open OUTPUT_FILE,'>>' . $output_dir . "chr" . $chr . $extention;
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

