#!/usr/bin/perl
use strict;
#use warnings;

#get the command line arguments
my $output_dir = $ARGV [0] . "/";


## hard coding
my $input_file = $output_dir ."extract.txt";
my $output_file = $output_dir . "summary_of_extract.txt";
my $full_report = 1;


#Global variables
my $highest_score = 0;
my $lowest_score = 9999;
my $highest_avg_score = 0;
my $lowest_avg_score = 9999;
my $average;
my $average_consensus;
my $no_of_proteins = 0;
my $total_no_of_proteins = 0;
my $alleles_available = 0;


my $total_number_of_epitopes = 0;
my $max_epitopes_per_protein = 0;
my $min_epitopes_per_protein = 99999;

my $total_number_of_alleles_used = 0;
my $max_allele_per_protein = 0;
my $min_allele_per_protein = 99999;


#create hashes
my %h_score;
my %h_avg_score; 
my %allele_count;
my %prom_allele_count;
my %method_used;
my %peptide_used;
my %h_total_no_of_bs;
my %h_max_no_bs;
my %h_list_max_no_bs;
my %h_comb_max_no_bs;
my %h_no_of_alleles_used;
my %h_ID_allele;
my %h_prom_allele;
my %h_order_prom_allele;



#Open debug file
#pen DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file;
open OUTPUT_FILE,'>' .  $output_file;


#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

	#Check for # or Accession
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);

	$line =~ s/\x0D//g; #strips ^M characters
	
	#Count the number of proteins
	$total_no_of_proteins++;
	
		

#1=UniProt ID :2=Min score 3=score average:4=MHC2 allele :5=Method used :6=Peptide sequence :7=Allele with most binding sites
#8=Max. number of binding sites  :9=Total number of binding sites for all alleles :10=Total number of MHC2 alleles available
#11=No. of MHC2 alleles used :12=No. of MHC2 alleles not used
#O45145	0.01	 215.57	HLA-DPA1*01/DPB1*04:01 	CombLib	LFHLLAGLA	HLA-DRB1*01:02	258	317	5	3	2	 PEPTIDE

	
	my ($ID, $score, $avg_score, $allele, $method, $peptide, $prom_allele, $max_no_bs, $total_no_of_bs,
	$total_no_of_alleles, $no_of_alleles_used,$no_of_alleles_not_used) = split ('\s+',$line);

	
	
	if ($allele eq '-') {
	
		print DEBUG_FILE $ID;
	
	} else {
		$allele_count {$allele} = $allele_count {$allele} + 1;
		$method_used {$method} = $method_used {$method} + 1;
		$h_ID_allele {$ID} = $allele;
	}
	
	if ($peptide ne '-') {
		$peptide_used {$peptide} = $peptide_used {$peptide} + 1;
	}
	
	if ( $prom_allele ne '-') {
	
		$prom_allele_count {$prom_allele} = $prom_allele_count {$prom_allele} + 1;
		
		$h_max_no_bs {$prom_allele} = $h_max_no_bs {$prom_allele} + $max_no_bs;
						
	
	}
	
	#ignore 9999 
	if  ($score != 9999) {
		$h_score {$ID} =  $score;
		$h_avg_score {$ID} =  $avg_score;
		$h_total_no_of_bs {$ID} = $total_no_of_bs;
		
		if ( $total_no_of_bs > $max_epitopes_per_protein) {
			$max_epitopes_per_protein = $total_no_of_bs;
		}
		
		if ( $total_no_of_bs < $min_epitopes_per_protein) {
			$min_epitopes_per_protein = $total_no_of_bs;
		}
		
		$total_number_of_epitopes = $total_number_of_epitopes + $total_no_of_bs;
		$h_list_max_no_bs {$ID} = $max_no_bs;
		$h_no_of_alleles_used {$ID} = $no_of_alleles_used;
		
		$total_number_of_alleles_used = $total_number_of_alleles_used + $no_of_alleles_used;
		
		if ( $no_of_alleles_used > $max_allele_per_protein) {
			$max_allele_per_protein = $total_no_of_bs;
		}
		
		if ( $no_of_alleles_used < $min_epitopes_per_protein) {
			$min_allele_per_protein = $total_no_of_bs;
		}		
		
			
		if ($score > $highest_score) {
			$highest_score = $score;
		}
		
		if ($score < $lowest_score) {
			$lowest_score = $score;
		}
		
		if ($avg_score > $highest_avg_score) {
			$highest_avg_score = $avg_score;
		}
		
		if ($avg_score < $lowest_avg_score) {
			$lowest_avg_score = $avg_score;
		}
		
		
		$average = $average + $score;
		$average_consensus = $average_consensus + $avg_score;
	}
	
	$no_of_proteins++;
	$alleles_available = $total_no_of_alleles;
			
} #end of LINE:while

my $allele_count = 0;

print OUTPUT_FILE "Total number of proteins input = " . $total_no_of_proteins . "\n";
print OUTPUT_FILE "Total number of MHC1 alleles available = " . $alleles_available . "\n";
print OUTPUT_FILE "Total number of epitopes = " . $total_number_of_epitopes . "\n";
print OUTPUT_FILE "\nAverage number of epitopes per protein = " . ($total_number_of_epitopes/$total_no_of_proteins) . "\n";
print OUTPUT_FILE "Max number of epitopes per protein = " . $max_epitopes_per_protein . "\n";
print OUTPUT_FILE "Min number of epitopes per protein = " . $min_epitopes_per_protein . "\n";
print OUTPUT_FILE "\Average number of alleles used per protein = " . ($total_number_of_alleles_used/$total_no_of_proteins) . "\n";
print OUTPUT_FILE "Max number of alleles per protein = " . $max_allele_per_protein . "\n";
print OUTPUT_FILE "Min number of alleles per protein = " . $min_allele_per_protein . "\n";


if ($full_report == 1) {

	print OUTPUT_FILE "\n### Frequency of MHC1 alleles used that bind to the lowest scoring epitope###\n";
	print OUTPUT_FILE "#1=Allele :2=Number of proteins with epitopes that bind to allele\n";
	foreach my $key ( sort  { $allele_count {$b} <=> $allele_count {$a}} keys  %allele_count) {

		print OUTPUT_FILE  $key . "\t" . $allele_count {$key} . "\n";
	
		$allele_count++;
	
	}
	print OUTPUT_FILE "Number of MHC1 alleles used = $allele_count\n";

	$allele_count = 0;
	
	

	print OUTPUT_FILE "\n\n### Frequency of prediction method used ###\n";
	foreach my $key ( sort  { $method_used {$b} <=> $method_used {$a}} keys  %method_used) {

		print OUTPUT_FILE  $key . "\t" . $method_used {$key} . "\n";
	
	}

	$average = $average/$no_of_proteins;
	$average_consensus = $average_consensus/$no_of_proteins;


	print OUTPUT_FILE "\n###IEDB IC50 scores based on mininum\n";
	print OUTPUT_FILE  "Highest score = " . $highest_score . " Lowest score = " . $lowest_score . " Average = " . $average . "\n";
	print OUTPUT_FILE "\n###IEDB IC50 scores based on average score of minimum percentile\n";
	print OUTPUT_FILE  "Highest score = " . $highest_avg_score . " Lowest score = " . $lowest_avg_score . " Average = " . $average_consensus . "\n";


	$allele_count = 0;

	print OUTPUT_FILE "\n### Frequency and average number of epitopes for promiscuous MHC1 alleles used ###\n";
	print OUTPUT_FILE "#1=Allele :2=Number of proteins with epitopes that bind to allele :3=average number of epitopes that bind to allele\n";
	foreach my $key ( sort  { $prom_allele_count {$b} <=> $prom_allele_count {$a}} keys  %prom_allele_count) {

		my $average = $h_max_no_bs {$key}/ $prom_allele_count {$key};

		print OUTPUT_FILE  $key . "\t" . $prom_allele_count {$key} . "\t";
		printf OUTPUT_FILE  ("%.1f\n"), $average;
	
		$allele_count++;
	
	}
	print OUTPUT_FILE "Number of promiscuous MHC1 alleles used = $allele_count\n";

	


	print OUTPUT_FILE "\n\n### Frequency of peptide used ###\n";
	foreach my $key ( sort  { $peptide_used {$b} <=> $peptide_used {$a}} keys  %peptide_used) {

		print OUTPUT_FILE  $key . "\t" . $peptide_used {$key} . "\n";
	
	}

	
} #end of full report


 

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
