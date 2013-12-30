#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAM "psl_to_fasta.pl"
#hard-coding
my $coverage_threshold = 80;

#get the command line arguments
my $prefix = $ARGV [0];
my $est_dir =  $ARGV [1] . "/";
my $chr_dir = $ARGV [2] . "/";
my $output_dir = $ARGV [3] . "/";
	

## hard coding
my $input_file = $output_dir . "chr" .  $prefix . "_est.psl";
my $chr_dna = $chr_dir . "chr" . $prefix . ".fasta";
my $output_file = $output_dir . "chr". $prefix . ".seq";


#Global variables
my $chr_seq;


#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';



#File handling
open INPUT_FILE,$input_file ;
open FASTA_FILE,$chr_dna;
open OUTPUT_FILE,'>' . $output_file;


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



#read the chromosome sequence
while (defined (my $line = <FASTA_FILE>)) {

	$line = &trim ($line);
	
	#Check for # 
	if ($line =~ m/^>/) { next;}
		
	$chr_seq = $chr_seq . $line;
		
} #end of LINE:while


my $count = 1;
#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}


	$line = &trim ($line);
	
	my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,
	$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts) = split ('\s+',$line);#Check for query coverage
	

	my $query_coverage = ((($qEnd - $qStart) / $qSize)) * 100;
	
	if ($query_coverage >= $coverage_threshold) {
	
		my $length = $tEnd - $tStart;
		my $seq = substr ($chr_seq,$tStart -1,$length + 1);
	
		my $id = "gene_" . $count++;
		print  OUTPUT_FILE ">$id|$qName|$tStart|$tEnd\n";
	
	
		#Reverse strand
		if ($strand eq "-") {
			$seq = &revdnacomp ($seq);
		}	
		
		print  OUTPUT_FILE $seq . "\n";
		 
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
sub revdnacomp {
  my $dna = shift @_;
  
  my $revcomp = reverse($dna);
  
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}
