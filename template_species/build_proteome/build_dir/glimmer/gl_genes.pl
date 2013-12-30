#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAM gl_genes.pl

# Get the comand arguments
my $prefix = $ARGV[0];
my $input_file = $ARGV[1] . "/" . $prefix . ".fasta";
my $output_path = $ARGV[2] . "/";

#input files
my $i_gff_file = $output_path . $prefix . ".gff";


#output files
my $seq_file = $output_path . $prefix . ".seq";
my $check_file =  $output_path . $prefix . "_check.seq";
my $exon_file = $output_path . $prefix . ".exons";
my $mrna_file = $output_path . $prefix . "_mrna.seq";
my $o_start_end_file = $output_path . $prefix  . "_gene_start_end.txt";


#Global variables
my $chromosome_seq;
my $correction;
my @exons = ();
my $not_first_time = 0;
my $exon_strand;
my $mrna_seq_id ;
my $mrna_seq;

#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open SEQ_FILE,'>'  . $seq_file;
open CHECK_FILE,'>' .  $check_file;
open EXON_FILE,'>'  . $exon_file;
open MRNA_FILE,'>'  . $mrna_file;
open  O_START_END_FILE,'>'  . $o_start_end_file;

#LOOP through input file
while (defined (my $line = <INPUT_FILE>)) {

	#Check for >
	if ($line =~ m/>/) { next;}

	$line = &trim ($line);
	
	$chromosome_seq = $chromosome_seq . $line;
			
} #end of LINE:while

open GFF_FILE,$i_gff_file ;

#LOOP through gff file
while (defined (my $line = <GFF_FILE>)) {

	

	#Check for #
	if ($line =~ m/^#/) { next;}

	$line = &trim ($line);

	#GFF format
	#TGME49_chrVIIb	GlimmerHMM	mRNA	650	7732	.	+	.	ID=TGME49_chrVIIb.path1.gene1;Name=TGME49_chrVIIb.path1.gene1

			
	my ($seqname,$source, $feature,$start,$end, $score, $strand, $frame,$attribute) = split ('\t',$line);

	if ($feature eq "mRNA"){

		#Print out the exons for previous gene
		if ($not_first_time == 1) {

			#print out exons
			&print_exons;
			
			#print mRNA seq
			&print_mrna_seq ($mrna_seq_id);

		}else {
			$not_first_time = 1;
		
		} # end of if ($not_first_time == 1)

		@exons = ();

		#get the gene number from the attribute 
		#ID=TGME49_chrVIIb.path1.gene1;Name=TGME49_chrVIIb.path1.gene1

		$attribute =~ m/gene(\d+)/;
		$attribute = "gene_" . $1;
	
		print  SEQ_FILE  ">" . $attribute  . "\n";
		print  MRNA_FILE  ">" . $attribute  . "\n";
		print O_START_END_FILE  $attribute . "\t" . $start . "\t" . $end . "\t" . $strand . "\n";
		print  EXON_FILE  ">" . $attribute  . "\t" . $start . "\t" . $end . "\n";
		
		my $length = $end - $start;
		my $seq = substr ($chromosome_seq,$start -1,$length + 1);

		if ($strand eq "+") {
				
			print  SEQ_FILE  $seq  . "\n";

		} elsif ($strand eq "-") {
	
			$seq = &revdnacomp ($seq);
			print  SEQ_FILE "$seq"  . "\n";

		}else {
			print  SEQ_FILE "Unknown direction\n";
		}

		#check if sequence starts with ATG
		if ($seq !~ m/^ATG/) {

			print  CHECK_FILE $attribute . "\t" . "Start=" . (substr ($seq,0,3)) . "\n";
		}
	
		#check if  end  sequence is either  TAG, TGA or TAA
		if ($seq !~ m/TAG$|TGA$|TAA$/) {

			print  CHECK_FILE  $attribute . "\t" . "End=" . (substr ($seq,$length - 2,3)) . "\n";
		}

		#Record current value of strand
		$exon_strand = $strand;
		$mrna_seq_id = $attribute;


	} elsif  ($feature eq "CDS") {

		push (@exons, $start . "," . $end );
		
		
		#build mrna sequence
		my $length = $end - $start;
		
		my $seq = substr ($chromosome_seq,$start -1,$length + 1);
	
		if ($strand eq "+") {
				
			$mrna_seq .= $seq;
	
		} elsif ($strand eq "-") {
	
			$mrna_seq .= $seq;
		}else {
			print  MRNA_FILE "#Unknown direction\n";
		}
	
	} 
			
} #end of LINE:while


#print out last set of exons and mRNA sequence
&print_exons;
&print_mrna_seq ($mrna_seq_id);

 

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

sub print_exons ($) {

	my $no_of_exons = scalar (@exons);
	my $count = 0;
	my $first = 1;
	my $correction = 0;


	#Reverse the array for complementary strand
	if ($exon_strand eq "-") {
		@exons = reverse(@exons);
	}

	foreach (@exons) {

		my $start_coords;
		my $end_coords;


		#Complementary strand
		if ($exon_strand eq "-") {

			($end_coords,  $start_coords ) = split (',',$_);

			if ($first == 1) {
				$correction = $start_coords + 1;
				$first = 0;
			}

			$start_coords = $correction - $start_coords;
			$end_coords = $correction - $end_coords;

		} else {

			($start_coords,$end_coords) = split (',',$_);

			if ($first == 1) {
				$correction = $start_coords - 1;
				$first = 0;
			}

			$start_coords = $start_coords - $correction;
			$end_coords = $end_coords - $correction;
		}

		$count ++;


		if ($no_of_exons == 1) {
			print  EXON_FILE "single exon" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 

		} elsif ($count == 1) {

			print  EXON_FILE "start" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
		
		} elsif ($count == $no_of_exons) {

			print  EXON_FILE "stop" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
		
		}else {
			print  EXON_FILE "exon" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
		
		}
	}

}

#print mRNA sequence
sub print_mrna_seq ($) {

	my $attribute = shift;
			
	if ($exon_strand eq "-") {
		$mrna_seq = &revdnacomp ($mrna_seq);
	}
	print  MRNA_FILE $mrna_seq . "\n";
		
	#check if sequence starts with ATG on mRNA
	if ($mrna_seq !~ m/^ATG/) {

		print  CHECK_FILE "mRNA_" . $attribute . "\t" . "Start=" . (substr ($mrna_seq,0,3)) . "\n";
	}

	my $seq_length = length ($mrna_seq);
	
	#check if  end  sequence is either  TAG, TGA or TAA on MRNA
	if ($mrna_seq !~ m/TAG$|TGA$|TAA$/) {

		print  CHECK_FILE  "mRNA_" . $attribute . "\t" . "End=" . (substr ($mrna_seq,$seq_length - 2,3)) . "\n";
	}
	
	$mrna_seq = "";
}
