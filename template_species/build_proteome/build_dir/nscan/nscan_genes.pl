#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAM nscan_genes.pl

# Get the comand arguments
my $prefix = $ARGV[0];
my $input_file = $ARGV[1] . "/" . $prefix . ".fasta";
my $output_path = $ARGV[2] . "/";

#input files
my $gtf_file_suffix = $prefix . ".fasta.masked.gtf";
my $gene_id_prefix = $prefix . ".fasta.";
my $i_gtf_file = $output_path . $gtf_file_suffix;


#output files
my $seq_file = $output_path . $prefix . ".seq";
my $mrna_file = $output_path . $prefix . "_mrna.seq";
my $check_file =  $output_path . $prefix . "_check.seq";


#Debug
#$prefix = "chrIa";
#$input_file = "chrIa.fasta";
#$i_gtf_file =  $prefix . ".fasta.masked.gtf";
#$seq_file = $prefix . ".seq";
#$mrna_file = $prefix . "_mrna.seq";
#$check_file =  $prefix . "_check.seq";

#Global variables
my $chromosome_seq;

#create hashes
my %h_mrna_seq;
my %h_strand;
my %h_start;
my %h_end;

#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open SEQ_FILE,'>'  . $seq_file;
open MRNA_FILE,'>'  . $mrna_file;
open CHECK_FILE,'>' .  $check_file;


#LOOP through input file
while (defined (my $line = <INPUT_FILE>)) {

	#Check for >
	if ($line =~ m/>/) { next;}

	$line = &trim ($line);
	
	$chromosome_seq = $chromosome_seq . $line;
			
} #end of LINE:while


open GFF_FILE,$i_gtf_file ;

#LOOP through gff file
while (defined (my $line = <GFF_FILE>)) {

	#Check for #
	if ($line =~ m/^#/) { next;}

	$line = &trim ($line);

	#GTF format
	#[seqname] [source] [feature] [start] [end] [score] [strand] [frame] [attributes]

			
	my ($seqname,$source, $feature,$start,$end, $score, $strand, $frame,$attribute) = split ('\t',$line);
	
	#get the gene number from the attribute 
	#chrIa.fasta.001
	$attribute =~ m/$gene_id_prefix(\d+)/;
	$attribute =  $1;

	#assemble the mrna seq
	if  (($feature eq "CDS") || ($feature eq "stop_codon")){

		#record the strand
		$h_strand {$attribute} = $strand;
		
		#build mrna sequence
		my $length = $end - $start;
		
		my $seq = substr ($chromosome_seq,$start -1,$length + 1);
		
		$h_mrna_seq {$attribute} = $h_mrna_seq {$attribute} . $seq;
	} 
	
	#Assemble the gene sequence
	if (($feature eq "CDS") && (&check_for_empty_string ($h_start {$attribute}))  && ($strand eq "+")){
			$h_start {$attribute} = $start;
	}elsif (($feature eq "stop_codon") && ($strand eq "+")) {
			$h_end {$attribute} = $end;
	}elsif (($feature eq "stop_codon") && ($strand eq "-")) {
			$h_start {$attribute} = $start;
	}elsif (($feature eq "start_codon") && ($strand eq "-")) {
			$h_end {$attribute} = $end;
	}elsif (($feature eq "CDS") && ($strand eq "-")) {
			$h_end {$attribute} = $end;
	}
	
} #end of LINE:while

#print gene sequence
foreach my $id (sort keys %h_start) {

	my $start = $h_start {$id};
	my $end = $h_end {$id};
	my $strand = $h_strand {$id};
	
	my $length = $end - $start;
		
	my $seq = substr ($chromosome_seq,$start -1,$length + 1);
	
	#check for reverse strand
	if ($h_strand {$id} eq "-") {
		$seq = &revdnacomp ($seq);
	}	

	#print sequence
	print SEQ_FILE ">gene_" . $id . "\n";
	print SEQ_FILE $seq . "\n";
				
	#check if sequence starts with ATG
	if ($seq !~ m/^ATG/) {

		print  CHECK_FILE $id . "\t" . "Start=" . (substr ($seq,0,3)) . "\n";
	}

	#check if  end  sequence is either  TAG, TGA or TAA
	if ($seq !~ m/TAG$|TGA$|TAA$/) {

		print  CHECK_FILE  $id . "\t" . "End=" . (substr ($seq,$length - 2,3)) . "\n";
	}
}


#print mRNA sequence
foreach my $id (sort keys %h_mrna_seq) {

	print MRNA_FILE ">gene_" . $id . "\n";
	my $seq = $h_mrna_seq {$id};

		
	#check for reverse strand
	if ($h_strand {$id} eq "-") {
		$seq = &revdnacomp ($h_mrna_seq {$id});
	}	

	#print sequence
	print MRNA_FILE $seq . "\n";
	
	my $length = length ($seq);
		
	#check if sequence starts with ATG
	if ($seq !~ m/^ATG/) {

		print  CHECK_FILE "mrna " . $id . "\t" . "Start=" . (substr ($seq,0,3)) . "\n";
	}

	#check if  end  sequence is either  TAG, TGA or TAA
	if ($seq !~ m/TAG$|TGA$|TAA$/) {

		print  CHECK_FILE  "mrna " . $id . "\t" . "End=" . (substr ($seq,$length - 2,3)) . "\n";
	}
	
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

sub revdnacomp {
  my $dna = shift @_;
  
  my $revcomp = reverse($dna);

  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}

