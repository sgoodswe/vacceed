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
my $allele_length_available = 0;


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
my %allele_length_count;
my %prom_allele_count;
my %prom_allele_length_count;
my %method_used;
my %peptide_used;
my %h_total_no_of_bs;
my %h_max_no_bs;
my %h_list_max_no_bs;
my %h_comb_max_no_bs;
my %h_no_of_alleles_used;
my %h_ID_allele;
my %h_ID_allele_length;
my %h_ID_prom_allele_length;
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
	
		
	#1=UniProt ID :2=Min score :3=MHC1 allele :4=Peptide length :5=species :6=Method used :7=Peptide sequence :8 & 9=Allele-length combination with most binding sites
	#10=Max. number of binding sites  :11=Total number of binding sites for all alleles :12=Total number of MHC1 allele-length combinations available
	#13=No. of MHC1 allele-length combinations used :14=No. of MHC1 allele-length combinations not used
	#O00933	1.3099089214	 HLA-A*31:01 9 human	SMM	RTYRHFSPR	HLA-A*30:01	9	6	65	303	40	207	 PEPTIDE
	
	
	my ($ID, $score, $avg_score, $allele, $length, $species, $method, $peptide, $prom_allele, $prom_length, $max_no_bs, $total_no_of_bs,
	$total_no_of_alleles, $no_of_alleles_used,$no_of_alleles_not_used,$status) = split ('\s+',$line);
	
	#print DEBUG_FILE "$ID $score $allele $length $species $method $peptide $prom_allele $prom_length $max_no_bs $total_no_of_bs $total_no_of_alleles $no_of_alleles_used $no_of_alleles_not_used $status\n";
	
	
	if ($allele eq '-') {
	
		print DEBUG_FILE $ID;
	
	} else {
		$allele_count {$allele} = $allele_count {$allele} + 1;
		$allele_length_count {$allele . "\t" . $length} = $allele_length_count {$allele . "\t" . $length} + 1;
		$method_used {$method} = $method_used {$method} + 1;
		$h_ID_allele {$ID} = $allele;
		$h_ID_allele_length {$ID} = $allele . "\t" . $length;
	}
	
	if ($peptide ne '-') {
		$peptide_used {$peptide} = $peptide_used {$peptide} + 1;
	}
	
	if ( $prom_allele ne '-') {
	
		$prom_allele_count {$prom_allele} = $prom_allele_count {$prom_allele} + 1;
		
		$h_max_no_bs {$prom_allele} = $h_max_no_bs {$prom_allele} + $max_no_bs;
		
		$prom_allele_length_count {$prom_allele . "\t" . $prom_length} = $prom_allele_length_count {$prom_allele . "\t" . $prom_length} + 1;
		
		$h_comb_max_no_bs {$prom_allele . "\t" . $prom_length} = $h_comb_max_no_bs {$prom_allele . "\t" . $prom_length} + $max_no_bs;
		
		$h_ID_prom_allele_length {$ID} = $prom_allele . "\t" . $prom_length;
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
	$allele_length_available = $total_no_of_alleles;
			
} #end of LINE:while

my $allele_count = 0;

print OUTPUT_FILE "Total number of proteins input = " . $total_no_of_proteins . "\n";
print OUTPUT_FILE "Total number of MHC1 allele combinations available = " . $allele_length_available . "\n";
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
	
	print OUTPUT_FILE "\n\n### Frequency of MHC1 alleles and length combinations used  that bind to the lowest scoring epitope###\n";
	print OUTPUT_FILE "#1=Allele :2=length :3=Number of proteins with epitopes that bind to allele\n";
	foreach my $key ( sort  { $allele_length_count {$b} <=> $allele_length_count {$a}} keys  %allele_length_count) {

		print OUTPUT_FILE  $key . "\t" . $allele_length_count {$key} . "\n";
	
		$allele_count++;
	
	}
	print OUTPUT_FILE "Number of MHC1 allele and length combinations used = $allele_count out of $allele_length_available\n";


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

	$allele_count = 0;
	
	print OUTPUT_FILE "\n\n### Frequency and average number of epitopes for promiscuous MHC1 alleles and length used ###\n";
	print OUTPUT_FILE "#1=Allele :2=length :3=Number of proteins with epitopes that bind to allele :4=average number of epitopes that bind to allele\n";
	foreach my $key ( sort  { $prom_allele_length_count {$b} <=> $prom_allele_length_count {$a}} keys  %prom_allele_length_count) {

	
		my $average = $h_comb_max_no_bs {$key}/ $prom_allele_length_count {$key};
		print OUTPUT_FILE  $key . "\t" . $prom_allele_length_count {$key} . "\t";
		printf OUTPUT_FILE  ("%.1f\n"), $average;
	
		$allele_count++;
	
	}

	print OUTPUT_FILE "Number of promiscuous MHC1 allele and length combinations used = $allele_count out of $allele_length_available\n";



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
