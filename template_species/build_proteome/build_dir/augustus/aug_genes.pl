#!/usr/bin/perl
use strict;
#use warnings;
### PROGRAM aug_genes.pl

# Get the comand arguments
my $prefix = $ARGV[0];
my $input_file = $ARGV[1] . "/" . $prefix . ".fasta";
my $output_path = $ARGV[2] . "/";


#input files
my $gff_file = $output_path . $prefix . ".gff";

#output files
my $seq_file = $output_path . $prefix . ".seq";
my $exon_file = $output_path . $prefix . ".exons";
my $check_file = $output_path . $prefix . "_check.seq";


#Global variables
my $chromosome_seq;
my $correction;
my @exons = ();
my $not_first_time = 0;
my $first = 1;
my $gene_id ;
my $exon_strand;

#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';

# To run program enter: 

#File handling
open INPUT_FILE,$input_file ;
open SEQ_FILE,'>'  .  $seq_file;
open CHECK_FILE,'>'  . $check_file;
open EXON_FILE,'>' .   $exon_file;

#LOOP through input file
while (defined (my $line = <INPUT_FILE>)) {

	#Check for >
	if ($line =~ m/>/) { next;}

	$line = &trim ($line);
	
	$chromosome_seq = $chromosome_seq . $line;

	
			
} #end of LINE:while



open GFF_FILE,$gff_file ;

#LOOP through gff file
while (defined (my $line = <GFF_FILE>)) {

	#Check for #
	if ($line =~ m/^#/) { next;}

	$line = &trim ($line);

	#GFF format
	#TGME49_chrVIIb	AUGUSTUS	gene	54055	62358	0.08	+	.	g1

	print DEBUG_FILE $line . "\n";
	
	my ($seqname,$source, $feature,$start,$end, $score, $strand, $frame,$attribute) = split ('\t',$line);
	

	if ($feature eq "gene"){


		#Print out the exons for previous gene
		if ($not_first_time == 1) {

			my $no_of_exons = scalar (@exons);
			my $count = 0;

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
					print  EXON_FILE "exon" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
	
				} elsif ($count == 1) {

					print  EXON_FILE "start" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
					
				} elsif ($count == $no_of_exons) {

					print  EXON_FILE "stop" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
					
				}else {
					print  EXON_FILE "exon" . "\t" . $start_coords  . "\t" .  $end_coords . "\n"; 
					
				}
		
			}



		}else {
			$not_first_time = 1;
		}

		my $seq;
		my $length;
		@exons = ();
		$first = 1;
		$correction = 0;
	
		
		print  SEQ_FILE  ">" . $attribute  . "\n";
		print  EXON_FILE  ">" . $attribute  . "\t" . $start . "\t" . $end . "\n";
		

		$length = $end - $start;
		$seq = substr ($chromosome_seq,$start -1,$length + 1);

		if ($strand eq "+") {
				
			print  SEQ_FILE  $seq  . "\n";

		} elsif ($strand eq "-") {
	
			$seq = &revdnacomp ($seq);
			# $seq = reverse ($seq);
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

			print  CHECK_FILE $strand . "\t" . $attribute . "\t" . "End=" . (substr ($seq,$length - 2,3)) . "\n";
		}

	
		#Record current values
		$exon_strand = $strand;


	} elsif  ($feature eq "CDS") {

		push (@exons, $start . "," . $end );

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
