#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAMconvert_dan_to_aa.pl

# Get the comand arguments
my $prefix = $ARGV[0];
my $output_path = $ARGV[1] . "/";
my $file_to_convert = $ARGV[2];

## hard coding
my $input_file = $output_path . $file_to_convert;
my $output_file = $output_path . $prefix . ".aa";


#Global variables
my $id;

#create hashes
my %h_dna;

#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;

#Read the DNA fasta file
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);
	#Check for # 
	if ($line =~ m/^>/) {

			$id = $line;
	} else {
	
		$h_dna {$id} .= $line; 
	}
		
} #end of LINE:while

#convert to amino acids
foreach my $id (sort keys %h_dna) {

	my $protein;
	my $codon;
	my $dna = $h_dna {$id};
	
	for(my $i=0;$i<(length($dna)-2);$i+=3)
	{
		$codon=substr($dna,$i,3);
		$protein .= &get_codon($codon);
	}

	print OUTPUT_FILE $id . "\n";
	print OUTPUT_FILE $protein . "\n";
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

sub get_codon{

	my($codon)=@_;
	
	$codon=uc $codon;
		
	my(%genetic_code) = (
		 'TCA'=>'S', #Serine
		 'TCC'=>'S', #Serine
		 'TCG'=>'S', #Serine
		 'TCT'=>'S', #Serine
		 'TTC'=>'F', #Phenylalanine
		 'TTT'=>'F', #Phenylalanine
		 'TTA'=>'L', #Leucine
		 'TTG'=>'L', #Leucine
		 'TAC'=>'Y', #Tyrosine
		 'TAT'=>'Y', #Tyrosine
		 'TAA'=>'', #Stop
		 'TAG'=>'', #Stop
		 'TGC'=>'C', #Cysteine
		 'TGT'=>'C', #Cysteine
		 'TGA'=>'', #Stop
		 'TGG'=>'W', #Tryptophan
		 'CTA'=>'L', #Leucine
		 'CTC'=>'L', #Leucine
		 'CTG'=>'L', #Leucine
		 'CTT'=>'L', #Leucine
		 'CCA'=>'P', #Proline
		 'CAT'=>'H', #Histidine
		 'CAA'=>'Q', #Glutamine
		 'CAG'=>'Q', #Glutamine
		 'CGA'=>'R', #Arginine
		 'CGC'=>'R', #Arginine
		 'CGG'=>'R', #Arginine
		 'CGT'=>'R', #Arginine
		 'ATA'=>'I', #Isoleucine
		 'ATC'=>'I', #Isoleucine
		 'ATT'=>'I', #Isoleucine
		 'ATG'=>'M', #Methionine
		 'ACA'=>'T', #Threonine
		 'ACC'=>'T', #Threonine
		 'ACG'=>'T', #Threonine
		 'ACT'=>'T', #Threonine
		 'AAC'=>'N', #Asparagine
		 'AAT'=>'N', #Asparagine
		 'AAA'=>'K', #Lysine
		 'AAG'=>'K', #Lysine
		 'AGC'=>'S', #Serine
		 'AGT'=>'S', #Serine
		 'AGA'=>'R', #Arginine
		 'AGG'=>'R', #Arginine
		 'CCC'=>'P', #Proline
		 'CCG'=>'P', #Proline
		 'CCT'=>'P', #Proline
		 'CAC'=>'H', #Histidine
		 'GTA'=>'V', #Valine
		 'GTC'=>'V', #Valine
		 'GTG'=>'V', #Valine
		 'GTT'=>'V', #Valine
		 'GCA'=>'A', #Alanine
		 'GCC'=>'A', #Alanine
		 'GCG'=>'A', #Alanine
		 'GCT'=>'A', #Alanine
		 'GAC'=>'D', #Aspartic Acid
		 'GAT'=>'D', #Aspartic Acid
		 'GAA'=>'E', #Glutamic Acid
		 'GAG'=>'E', #Glutamic Acid
		 'GGA'=>'G', #Glycine
		 'GGC'=>'G', #Glycine
		 'GGG'=>'G', #Glycine
		 'GGT'=>'G', #Glycine
	 );
	
	if(exists $genetic_code{$codon}){
		return $genetic_code{$codon};
	} else {
		print  "Bad codon \"$codon\"!!\n";
	}
}
