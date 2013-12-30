#!/usr/bin/perl
use strict;
#use warnings;

### PROGRAM "e_blast.pl"
#hard-coding
my $l_query_coverage_threshold = 75;
my $u_query_coverage_threshold = 125;
my $est_query_coverage_threshold = 100;
my $ident_threshold = 95;
my $gene_file_prefix = "genes_";

#get the command line arguments
my $prefix = $ARGV [0];
my $gene_dir =  $ARGV [1] . "/";
my $output_dir = $ARGV [2] . "/";
my $assembly_dir = $ARGV [3] . "/";
my $sum_dir = $ARGV [4] . "/";
my $resource = $ARGV [5];
my $EST = $ARGV [6];

#Check if we are analysing a blast file generated from ESTs  ($EST=1 for ESTs)
if (&check_for_empty_string ($EST)){
	$EST = 0;
}


#input files
my $i_test_file = $gene_dir . $gene_file_prefix . $prefix .  ".fasta";
my $input_file = $output_dir . $prefix . "_blastn.txt";
my $i_fasta_file =  $output_dir . $prefix . ".seq";


#output_files
my $o_map_file = $output_dir . $prefix . "_map.txt";
my $o_summary_file = $output_dir  . $prefix . "_summary.txt";
my $o_genome_summary =  $sum_dir . "genome_summary.txt";
my $o_gene_id_file = $assembly_dir . $prefix . "_" . $resource . "_gene_IDs.txt";
my $o_new_gene_file = $assembly_dir . $prefix . "_". $resource . "_new_gene.fasta";


#Global variables
my $count = 0;
my $total_count = 0;
my $no_of_test_genes = 0;
my $no_of_new_genes = 0;
my $no_of_non_matched_NCBI_genes = 0;


#create hashes
my %h_predicted_genes;
my %h_predicted_fasta;
my %h_test_genes;
my %h_match;
my %h_queries;
my %h_sorted_query_file;
my %h_bitscores;

# E-value = 0
my %h_all_0;
my %h_query_coverage_0;
my %h_ident_0;
# E-value not 0
my %h_all;
my %h_query_coverage;
my %h_ident;

#not matched
my %h_non_match;

my %h_e_value;
my %h_new_gene;
my %h_test_gene_match;
my %h_partial_match;



#Open debug file
#open DEBUG_FILE,'>debug.txt';

# To run program enter: 

#Open and read the fasta sequence file for the predicted genes
open FASTA_FILE,$i_fasta_file ;

my $id;

#LOOP through input file
while (defined (my $line = <FASTA_FILE>)) {

	$line = &trim ($line);

	if ($line =~ m/^>(\w+)/) {
		
		$id = $1;
		$h_predicted_genes {$id} = $1;

	} else {
	
		$h_predicted_fasta {$id} =  $line;
	}

}

#Open and read the test sequence file for the validated genes
open TEST_FILE,$i_test_file ;

#LOOP through input file
while (defined (my $line = <TEST_FILE>)) {

	$line = &trim ($line);

	if ($line =~ m/^>(\w+)/) {
		
		$no_of_test_genes++;

		$h_test_genes {$1} = $1;
	
	}
}

#File handling
open INPUT_FILE,$input_file ;
open SUMMARY_FILE,'>'  . $o_summary_file;

open MAP_FILE,'>'  . $o_map_file;
open GENE_ID_FILE,'>'  . $o_gene_id_file;
open NEW_GENE_FILE,'>'  . $o_new_gene_file;


#check ifsummary file exists
if (-e $o_genome_summary){
	#add to existing file
	open GENOME_SUMMARY_FILE,'>>'  . $o_genome_summary;
} else {
	open GENOME_SUMMARY_FILE,'>'  . $o_genome_summary;
	print GENOME_SUMMARY_FILE "#1:chromosome 2:100% query coverage and ident (E-value = 0) 3:Partial match \n"; 
	print GENOME_SUMMARY_FILE "#4:Total number of matched genes #5: number of duplicate matches 6:number of predicted genes 7:number of NCBI genes 8: matches - duplicates\n";
	print GENOME_SUMMARY_FILE "#9:number of new genes #10: number of non-matched NCBI genes\n";
	print GENOME_SUMMARY_FILE "#11:%Same length,coverage,matches 12:%Same length and coverage 13:%Same length and  matches 14:%Same coverage and matches\n"; 
	print GENOME_SUMMARY_FILE "#15:%Total number of matched genes #16:%number of duplicate matches 17:%matches - duplicates\n";
}

print GENOME_SUMMARY_FILE $prefix . ",";

#LOOP through the blast file and extract the highest bit score for each query
LOOP:while (defined (my $line = <INPUT_FILE>)) {
	
	$line = &trim ($line);

	#blast output
	#0.0,g1.t1,803,NCLIV_000010,818,1670,4325,100.00,803,803
	#"evalue qseqid qlen sseqid slen bitscore score pident nident mismatches positive"
	
	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $line);


	#extract the number from the query id
	$query_id =~ m/(\d+)/;
	my $id = $1;
	
	#Get the highest bit score 
	if (&check_for_empty_string ($h_bitscores {$id})) {
			$h_sorted_query_file {$id} = $line;
			$h_bitscores {$id} = $bitscore;
	} elsif ($bitscore > $h_bitscores {$id}) {
			$h_sorted_query_file {$id} = $line;
			$h_bitscores {$id} = $bitscore;
	}
	
}

foreach my $id (sort keys %h_sorted_query_file) {

	#print DEBUG_FILE $h_sorted_query_file {$id} . "\n";
}


#LOOP through query hash
foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});

	my $query_coverage;
	
	if ($EST == 1) {
		#Query length - $matches = query covergae
		#Ident = mismatches, - matches
		$query_coverage = (($query_length - $matches)/$query_length) * 100;
		$query_coverage = sprintf("%.2f",100 - $query_coverage);
		
	} else {
		#Query length - subject length = query covergae
		#Ident = mismatches, - matches
		$query_coverage = (($subject_length - $query_length)/$subject_length) * 100;
		$query_coverage = sprintf("%.2f",100 - $query_coverage);
	
	}
	
	#extract the first part from the query id
	$query_id =~ m/(\w+)/;
	$query_id = $1;
	
	#Record the Query coverage and ident
	$h_queries {$query_id} = $query_coverage . "\t" . $ident;
	
	if ($evalue == 0) {
	
		#100% query coverage and 100% ident	
		if (($query_length == $subject_length) && ($ident >= $ident_threshold))  {

			$h_all_0 {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
						
		#Partial match if query coverage > than lower threshold and < upper threshold
		} elsif (($query_coverage >= $l_query_coverage_threshold) && ($query_coverage <= $u_query_coverage_threshold) && ($EST == 0) ) {
		
			$h_query_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
		
		## Perfect  ESTs (100% coverage and 100% ident)
		} elsif (($EST == 1) &&($query_coverage >= $est_query_coverage_threshold) && ($ident >= $ident_threshold)) {
	
			$h_all_0 {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
			
		#For ESTs check for 100% coverage
		} elsif (($EST == 1) &&($query_coverage >= $est_query_coverage_threshold)) {
			$h_query_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
			
		#For ESTs check for 100% ident
		} elsif (($EST == 1) && ($ident >= $ident_threshold) ) {
	
			$h_ident_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;

		} else {
		
			$h_non_match {$query_id} = $subject_id;
		}

	#Matches an existing gene  but e-value not 0
	} else {
	
		#ALL	
		if (($query_length == $subject_length) && ($ident >= $ident_threshold))  {

			$h_all {$query_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_partial_match {$query_id} = $subject_id;
		
		#Partial match if query coverage > than lower threshold and < upper threshold
		} elsif (($query_coverage >= $l_query_coverage_threshold) && ($query_coverage <= $u_query_coverage_threshold) && ($EST == 0) ) {
		
			$h_query_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
		
		##  ESTs (100% coverage and 100% ident)		
		} elsif (($EST == 1) &&($query_coverage >= $est_query_coverage_threshold) && ($ident >= $ident_threshold)) {
	
			$h_all {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
		
		#For ESTs check for 100% coverage
		} elsif (($EST == 1) &&($query_coverage >= $est_query_coverage_threshold)) {
	
			$h_query_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;
			
		#For ESTs check for 100% ident
		} elsif (($EST == 1) && ($ident >= $ident_threshold) ) {
	
			$h_ident {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_gene_match {$subject_id} = $query_id;
			$h_e_value {$query_id} = $subject_id;

		}else {
		
			$h_non_match {$query_id} = $subject_id;
		}

	}
			
} #end of LINE:while

### E-Value = 0

my $all = 0;

print  SUMMARY_FILE "# E-value 0 - 100% query coverage and ident\n";
foreach my $key (sort (keys %h_all_0)) {
	$all++;
	if ($key =~ m/gene/) {
		print  SUMMARY_FILE $key . "\t" . $h_all_0 {$key}  ."\n";
	} else {
		print  SUMMARY_FILE "gene_" . $key . "\t" . $h_all_0 {$key}  ."\n";
	}
}
print  SUMMARY_FILE "#\tNumber of genes: " . $all . "\n";


$total_count = $total_count + $all;

if ($EST == 0) {
	print  SUMMARY_FILE "# E-value 0 - query coverage > $l_query_coverage_threshold% and < $u_query_coverage_threshold%\n";
} else {
	print  SUMMARY_FILE "# E-value 0 - query coverage = $est_query_coverage_threshold%\n";
}
my $query_coverage = 0;

foreach my $key (sort (keys %h_query_coverage_0)) {
	$query_coverage++;
	if ($key =~ m/gene/) {
		print  SUMMARY_FILE  $key . "\t" . $h_query_coverage_0 {$key} . "\n"
	} else {
		print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_query_coverage_0 {$key} . "\n"	
	}
}

print  SUMMARY_FILE "#\tNumber of genes: " . $query_coverage . "\n";

$total_count = $total_count + $query_coverage;

#Check for ESTs
if ($EST == 1) {

	print  SUMMARY_FILE "# E-value 0 - >= $ident_threshold% ident (but part of EST is not matched)\n";
	my $ident_coverage = 0;

	foreach my $key (sort (keys %h_ident_0)) {
		$ident_coverage++;
		
		if ($key =~ m/gene/) {
			print  SUMMARY_FILE  $key . "\t" . $h_ident_0 {$key} . "\n"
		} else {
			print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_ident_0 {$key} . "\n"		
		}
	}

	print  SUMMARY_FILE "#\tNumber of genes: " . $ident_coverage . "\n";

	$total_count = $total_count + $ident_coverage;

}


### E-Value NOT 0

$all = 0;

print  SUMMARY_FILE "# E-value NOT 0 - 100% query coverage and ident\n";
foreach my $key (sort (keys %h_all)) {
	$all++;
	if ($key =~ m/gene/) {
		print  SUMMARY_FILE  $key . "\t" . $h_all {$key} . "\n";
	} else {
		print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_all {$key} . "\n";	
	}
}
print  SUMMARY_FILE "#\tNumber of genes: " . $all . "\n";
$total_count = $total_count + $all;

if ($EST == 0) {
	print  SUMMARY_FILE "# E-value NOT 0 - query coverage > $l_query_coverage_threshold% and < $u_query_coverage_threshold%\n";
} else {
	print  SUMMARY_FILE "# E-value NOT 0 - query coverage = $est_query_coverage_threshold%\n";
}
$query_coverage = 0;

foreach my $key (sort (keys %h_query_coverage)) {
	$query_coverage++;
	if ($key =~ m/gene/) {
		print  SUMMARY_FILE  $key . "\t" . $h_query_coverage {$key} . "\n"
	} else {
		print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_query_coverage {$key} . "\n"	
	}
}

print  SUMMARY_FILE "#\tNumber of genes: " . $query_coverage . "\n";

$total_count = $total_count + $query_coverage;

#Check for ESTs
if ($EST == 1) {

	print  SUMMARY_FILE "# E-value NOT 0 - >= $ident_threshold% ident (but part of EST is not matched)\n";
	my $ident_coverage = 0;

	foreach my $key (sort (keys %h_ident)) {
		$ident_coverage++;
		if ($key =~ m/gene/) {
			print  SUMMARY_FILE  $key . "\t" . $h_ident {$key} . "\n"
		} else {
			print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_ident {$key} . "\n"
		}
	}

	print  SUMMARY_FILE "#\tNumber of genes: " . $ident_coverage . "\n";

	$total_count = $total_count + $ident_coverage;

}

print  SUMMARY_FILE "#NON matched genes\n";
my $non_match = 0;

foreach my $key (sort (keys %h_non_match)) {
	$non_match++;
	if ($key =~ m/gene/) {
		print  SUMMARY_FILE  $key . "\t" . $h_non_match {$key} . "\n"
	} else {
		print  SUMMARY_FILE  "gene_" . $key . "\t" . $h_non_match {$key} . "\n"	
	}
}

print  SUMMARY_FILE "#\n#\n# Total number of matched genes: " . $total_count . "\n";
print  SUMMARY_FILE "# Total number of non matched genes: " . $non_match . "\n";


#Create map file
my $no_of_predicted_genes  = 0;
my $no_of_matched_genes  = 0;
my $no_of_duplicated_genes  = 0;

foreach my $predicted_gene (sort (keys %h_predicted_genes)) {

	 $no_of_predicted_genes++;

	if (&check_for_empty_string ($h_e_value {$predicted_gene})) {

		print  MAP_FILE $predicted_gene . "\t" . "NO_MATCH" . "\n";
		$h_new_gene {$predicted_gene} = $predicted_gene;

		#print possible new genes to a fasta file
		print NEW_GENE_FILE ">" . $resource . "-" . $predicted_gene . "\n";
		print NEW_GENE_FILE $h_predicted_fasta {$predicted_gene} . "\n";
		
		
	} else {

				
		my $test_gene = $h_e_value {$predicted_gene};

		
		#check for duplicate match
		if (my $duplicate_gene = &check_for_duplicate_prediction ($test_gene)) {
		
			print  MAP_FILE $predicted_gene . "\t" . $h_e_value {$predicted_gene} . "\t" . "Duplicate:" . $duplicate_gene . "\n";
			$no_of_duplicated_genes++;

		} else {

			$h_match {$test_gene} = $predicted_gene;
			print  MAP_FILE $predicted_gene . "\t" . $h_e_value {$predicted_gene} . "\n";
			$no_of_matched_genes++;	

		}
	}
}

print  MAP_FILE "#\n#\n# Total number of predicted genes: " .  $no_of_predicted_genes . "\n";
print  MAP_FILE "# Total number of test genes: " .  $no_of_test_genes . "\n";
print  MAP_FILE "# Total number of matched genes: " .  $no_of_matched_genes . "\n";
print  MAP_FILE "# Total number of duplicated genes: " .  $no_of_duplicated_genes . "\n";


### Gene_id file
my $perfect = 0;
my $partial = 0;
my $no_of_non_matched_NCBI_genes = 0;

print  GENE_ID_FILE "#1:UniProt ID 2:Query coverage (100 = perfect) 3:Ident (100 = perfect) 4:Query ID\n";
print  GENE_ID_FILE "# PERFECT MATCH\n";
foreach my $key (sort keys %h_all_0) {
	print GENE_ID_FILE $h_all_0 {$key} . "\t" . $h_queries {$key} . "\t" . $key . "\n";
	$perfect++;
}

print  GENE_ID_FILE "# PARTIAL MATCH\n";
foreach my $key (sort (keys %h_partial_match)) {
	print  GENE_ID_FILE   $h_partial_match {$key} . "\t" . $h_queries {$key} . "\t" .  $key . "\n";
	$partial++;
}

print  GENE_ID_FILE "# NON-MATCHED NCBI GENES\n";
foreach my $key (sort (keys %h_test_genes)) {

if (&check_for_empty_string ($h_test_gene_match {$key})) {
		print  GENE_ID_FILE $key  . "\n";
		$no_of_non_matched_NCBI_genes++;
	}
}

print  GENE_ID_FILE "# NEW Genes\n";
foreach my $key (sort (keys %h_new_gene)) {
	print GENE_ID_FILE  $resource . "-" . $h_new_gene {$key}  . "\t" . $h_queries {$key} . "\n";
	$no_of_new_genes++;
}


### Genome summary file
if ($no_of_predicted_genes == 0) {
	print "ERROR: Number of predictions = 0\n";
	print "e-blastn.pl script was terminated\n";
	exit 1;
}

my $percent_perfect = ($perfect/$no_of_predicted_genes) * 100;
my $percent_partial = ($partial/$no_of_predicted_genes) * 100;
my $percent_total_count = ($total_count/$no_of_predicted_genes) * 100;
my $percent_no_of_duplicated_genes = ($no_of_duplicated_genes/$no_of_predicted_genes) * 100;
my $percent_no_of_matched_genes = ($no_of_matched_genes/$no_of_predicted_genes) * 100;



print GENOME_SUMMARY_FILE "$perfect,$partial,$total_count,$no_of_duplicated_genes,$no_of_predicted_genes,$no_of_test_genes,$no_of_matched_genes,";
print GENOME_SUMMARY_FILE "$no_of_new_genes,$no_of_non_matched_NCBI_genes,";

printf GENOME_SUMMARY_FILE ("%.1f,%.1f,%.1f,%.1f,%.1f\n",$percent_perfect,$percent_partial,$percent_total_count,$percent_no_of_duplicated_genes ,$percent_no_of_matched_genes); 


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


#Check for empty string
sub  check_for_duplicate_prediction ($)
{

	my $test_gene = shift;
	
	foreach my $key (sort (keys %h_match)) {

		if ($key eq $test_gene) {

			return $h_match {$key};
		}

	}

	return;

	
}
