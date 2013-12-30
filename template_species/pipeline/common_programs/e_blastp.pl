#!/usr/bin/perl
use strict;
#use warnings;


## hard coding
my $l_query_coverage_threshold = 75;
my $u_query_coverage_threshold = 125;
my $ident_threshold = 75;


#get the command line arguments
my $prefix = $ARGV [0];
my $prot_dir =  $ARGV [1] . "/";
my $sum_dir = $ARGV [2] . "/";
my $output_dir = $ARGV [3] . "/";
my $assembly_dir = $ARGV [4] . "/";
my $map_dir =  $ARGV [5] . "/";
my $prot_file =  $ARGV [6];
my $resource = $ARGV [7];



### More hard coding
## Warning the following will be different with a different test fasta file
## The fasta header maybe for example >tr|UniProt_ID OR >gi|NCBI_ID
my $subject_prefix_id = "tr";

### mapping files should be in mapping directory
my $i_chr_map =  $map_dir . "map_uniprot_gene_chr.txt";
my $i_map_uniprot_gi = $map_dir . "map_uniprot_gi.txt";

my $i_test_file =  $prot_dir . $prot_file;
my $input_file = $output_dir . $prefix . "_blastp.txt";
my $i_fasta_file =  $output_dir .  $prefix . ".aa";

my $o_map_file = $output_dir .  $prefix . "_prot_map.txt";
my $o_summary_file = $output_dir . $prefix .  "_prot_summary.txt";

my $o_prot_summary = $sum_dir .   "protein_summary.txt";

my $o_prot_id_file = $assembly_dir .  $prefix . "_" . $resource . "_prot_IDs.txt";
my $o_new_prot_file = $assembly_dir . $prefix . "_" .  $resource . "_new_prot.fasta";


#Global variables
my $count = 0;
my $total_count = 0;
my $no_of_test_proteins = 0;
my $no_of_predicted_proteins = 0;
my $no_of_new_proteins = 0;
my $no_of_non_matched_NCBI_genes = 0;


#create hashes
my %h_predicted_proteins;
my %h_test_proteins;
my %h_match;
my %h_sorted_query_file;
my %h_bitscores;
my %h_queries;
my %h_protein_chr_map;
my %h_gi_uniprot_map;

# E-value = 0
my %h_all_0;
my %h_query_coverage_0;
my %h_ident_0;
my %h_partial_coverage_0;
# E-value not 0
my %h_all;
my %h_query_coverage;
my %h_ident;
my %h_partial_coverage;

#not matched
my %h_non_match;

my %h_e_value;
my %h_new_protein;
my %h_test_protein_match;
my %h_partial_match;
my %h_predicted_fasta;

#Open debug file
#open DEBUG_FILE,'>debug.txt';

#Open and read the gene to chromosome map file
open CHR_MAP_FILE,$i_chr_map ;



#LOOP through  file
while (defined (my $line = <CHR_MAP_FILE>)) {

	$line = &trim ($line);

	
	my ($uniprot_id,$toxodb_id,$chr) = split ('\s+',$line);
	
	$h_protein_chr_map {$uniprot_id} = $chr;

	#Get chromosome number from prefix
	$prefix =~ m/chr(\w+)/;

	my $chr_no = $1;
	
	if ($chr eq $chr_no) {
		$no_of_test_proteins++;
	}

}

#only do the following if a subject ID prefix is a gi for NCBI GI number
if ($subject_prefix_id eq (lc ("gi"))) {

	open MAP_UNIPROT_GI_FILE,$i_map_uniprot_gi ;

	#LOOP through  file
	while (defined (my $line = <MAP_UNIPROT_GI_FILE>)) {

		$line = &trim ($line);

		my ($uniprot_id,$gi) = split ('\s+',$line);
		
		$h_gi_uniprot_map {$gi} = $uniprot_id;
	
	}
}

#Open and read the fasta sequence file for the predicted proteins
open FASTA_FILE,$i_fasta_file ;

my $id;


#LOOP through input file
while (defined (my $line = <FASTA_FILE>)) {

	$line = &trim ($line);
	
	if ($line =~ m/^>(\w+)/) {
	
		$id = $1;
		$h_predicted_proteins {$id} = $1;

	} else {
	
		$h_predicted_fasta {$id} = $h_predicted_fasta {$id} . $line;
	}


}

#Open and read the test sequence file for test proteins
open TEST_FILE,$i_test_file ;

#LOOP through input file
while (defined (my $line = <TEST_FILE>)) {

	$line = &trim ($line);

	if ($line =~ m/^>tr.(\w+)/) {
	
		$h_test_proteins {$1} = $1;
	
	}
}

#File handling
open INPUT_FILE,$input_file ;
open SUMMARY_FILE,'>'  . $o_summary_file;
open MAP_FILE,'>'  . $o_map_file;
open PROT_ID_FILE,'>'  . $o_prot_id_file;
open NEW_PROT_FILE,'>'  . $o_new_prot_file;


#check ifsummary file exists
if (-e $o_prot_summary) {
	#add to existing file
	open PROTEIN_SUMMARY_FILE,'>>'  . $o_prot_summary;
} else {	
	open PROTEIN_SUMMARY_FILE,'>'  . $o_prot_summary;
	print PROTEIN_SUMMARY_FILE "#1:chromosome 2:100% query coverage and ident (E-value = 0) 3:Partial match 4:Total number of matched proteins #5: number of duplicate matches\n"; 
	print PROTEIN_SUMMARY_FILE "#6:number of predicted proteins #7:number of expected poroteins for chr $prefix #8: matches - duplicates #9:number of possible new proteins\n";
	print PROTEIN_SUMMARY_FILE "#10: number of non-matched UniProt proteins #11:% of 100% query coverage and ident (E-value = 0) 12:%partial\n"; 
	print PROTEIN_SUMMARY_FILE "#13:%Total number of matched proteins #14:%number of duplicate matches 15:%matches - duplicates\n";
	
} 
print PROTEIN_SUMMARY_FILE $prefix . ",";


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


foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});
	
	#extract only alpha numberic chaharacters from the query id
	$query_id =~ m/(\w+)/;
	$query_id = $1;
	

	#Query length - subject length = query covergae
	#Ident = mismatches, - matches
	my $query_coverage = (($subject_length - $query_length)/$subject_length) * 100;
	$query_coverage = sprintf("%.2f",100 - $query_coverage);
	
	#Record the Query coverage and ident
	$h_queries {$query_id} = $query_coverage . "\t" . $ident;
	
	#extract the subject_id e.g.  tr|F0V712|F0V712_NEOCL or gi|401395063|ref|XP_003879545.1|
	if ($subject_id =~ m/$subject_prefix_id.(\w+)/) {
		$subject_id = $1;
	}
	
	#Map GI to UniProt if $subject_prefix_id is equal to gi
	if ($subject_prefix_id eq (lc ("gi"))) {
			
			if (&check_for_empty_string ($h_gi_uniprot_map {$subject_id})) {
				$subject_id = "NO_MAP";
			} else {
				$subject_id = $h_gi_uniprot_map {$subject_id};
			}
	}
		
	if ($evalue == 0) {
	
		#ALL	
		if (($query_length == $subject_length) && ($ident == 100))  {

			$h_all_0 {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
		
		} elsif ($query_length == $subject_length) {
	
			$h_query_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;

		} elsif ($ident == 100){

			$h_ident_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
			
		#Partial match if query coverage > than lower threshold and < upper threshold AND > ident threshold
		} elsif (($query_coverage > $l_query_coverage_threshold) && ($query_coverage < $u_query_coverage_threshold) && ($ident > $ident_threshold)) {

			$h_partial_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;

		} else {
		
			$h_non_match {$query_id} = $subject_id;
		}

	#Matches an existing protein  but e-value not 0
	} else {
	
		#ALL	
		if (($query_length == $subject_length) && ($ident == 100))  {

			$h_all {$query_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;;
		
		} elsif ($query_length == $subject_length) {

			$h_query_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;

		} elsif ($ident == 100){

			$h_ident {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			
		#Partial match if query coverage > than lower threshold and < upper threshold AND > ident threshold
		} elsif (($query_coverage > $l_query_coverage_threshold) && ($query_coverage < $u_query_coverage_threshold) && ($ident > $ident_threshold)) {

			$h_partial_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			$h_test_protein_match {$subject_id} = $subject_id;
			$h_e_value {$query_id} = $subject_id;

		} else {
		
			$h_non_match {$query_id} = $subject_id;
		}

	}
			
} #end of LINE:while


### E-Value = 0

my $all = 0;

print  SUMMARY_FILE "# E-value 0 - 100% query coverage and ident\n";
foreach my $key (sort (keys %h_all_0)) {

	$all++;

	if ($h_protein_chr_map {$h_all_0 {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_all_0 {$key} . "\t" . $h_protein_chr_map {$h_all_0 {$key}} ."\n";
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_all_0 {$key}  ."\n";
	}


}
print  SUMMARY_FILE "#\tNumber of proteins: " . $all . "\n";


$total_count = $total_count + $all;


print  SUMMARY_FILE "# E-value 0 - 100% query coverage\n";
my $query_coverage = 0;

foreach my $key (sort (keys %h_query_coverage_0)) {

	$query_coverage++;
	if ($h_protein_chr_map {$h_query_coverage_0 {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_query_coverage_0 {$key} . "\t" . $h_protein_chr_map {$h_query_coverage_0 {$key}} . "\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_query_coverage_0 {$key} . "\n"
	}

}

print  SUMMARY_FILE "#\tNumber of proteins: " . $query_coverage . "\n";

$total_count = $total_count + $query_coverage;

print  SUMMARY_FILE "# E-value 0 - 100% ident\n";
my $ident = 0;

foreach my $key (sort (keys %h_ident_0)) {

	$ident++;

	if ($h_protein_chr_map { $h_ident_0 {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_ident_0 {$key} . "\t" . $h_protein_chr_map { $h_ident_0 {$key}} ."\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_ident_0 {$key} . "\n"
	}


}

print  SUMMARY_FILE "#\tNumber of proteins: " . $ident . "\n";

$total_count = $total_count + $ident;

print  SUMMARY_FILE "# E-value 0 - query coverage > $l_query_coverage_threshold and < $u_query_coverage_threshold AND ident > $ident_threshold\n";
my $partial_coverage = 0;

foreach my $key (sort (keys %h_partial_coverage_0)) {

	$partial_coverage++;

	if ($h_protein_chr_map { $h_partial_coverage_0 {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_partial_coverage_0 {$key} . "\t" . $h_protein_chr_map { $h_partial_coverage_0 {$key}} ."\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_partial_coverage_0 {$key} . "\n"
	}


}

print  SUMMARY_FILE "#\tNumber of proteins: " . $partial_coverage . "\n";

$total_count = $total_count + $partial_coverage;


### E-Value NOT 0

$all = 0;

print  SUMMARY_FILE "# E-value NOT 0 - 100% query coverage and ident\n";
foreach my $key (sort (keys %h_all)) {

	$all++;

	if ($h_protein_chr_map {$h_all {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_all {$key} . "\t" . $h_protein_chr_map {$h_all{$key}} ."\n";
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_all {$key}  ."\n";
	}


}
print  SUMMARY_FILE "#\tNumber of proteins: " . $all . "\n";
$total_count = $total_count + $all;


print  SUMMARY_FILE "# E-value NOT 0 - 100% query coverage\n";
$query_coverage = 0;

foreach my $key (sort (keys %h_query_coverage)) {

	$query_coverage++;
	if ($h_protein_chr_map {$h_query_coverage {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_query_coverage {$key} . "\t" . $h_protein_chr_map {$h_query_coverage {$key}} . "\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_query_coverage {$key} . "\n"
	}

}

print  SUMMARY_FILE "#\tNumber of proteins: " . $query_coverage . "\n";

$total_count = $total_count + $query_coverage;

print  SUMMARY_FILE "# E-value NOT 0 - 100% ident\n";
$ident = 0;

foreach my $key (sort (keys %h_ident)) {

	$ident++;

	if ($h_protein_chr_map { $h_ident {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_ident {$key} . "\t" . $h_protein_chr_map { $h_ident {$key}} ."\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_ident {$key} . "\n"
	}


}

print  SUMMARY_FILE "#\tNumber of proteins: " . $ident . "\n";

$total_count = $total_count + $ident;

print  SUMMARY_FILE "# E-value NOT 0 - query coverage > $l_query_coverage_threshold and < $u_query_coverage_threshold AND ident > $ident_threshold\n";
my $partial_coverage = 0;

foreach my $key (sort (keys %h_partial_coverage)) {

	$partial_coverage++;

	if ($h_protein_chr_map { $h_partial_coverage {$key}} ne $prefix) {
		print  SUMMARY_FILE $key . "\t" . $h_partial_coverage {$key} . "\t" . $h_protein_chr_map { $h_partial_coverage {$key}} ."\n"
	} else {
		print  SUMMARY_FILE $key . "\t" . $h_partial_coverage {$key} . "\n"
	}


}

print  SUMMARY_FILE "#\tNumber of proteins: " . $partial_coverage . "\n";

$total_count = $total_count + $partial_coverage;

print  SUMMARY_FILE "#NON matched proteins\n";
my $non_match = 0;

foreach my $key (sort (keys %h_non_match)) {
	$non_match++;
	print  SUMMARY_FILE $key . "\t" . $h_non_match {$key} . "\n"
}


print  SUMMARY_FILE "#\n#\n# Total number of matched proteins: " . $total_count . "\n";
print  SUMMARY_FILE "# Total number of non matched proteins: " . $non_match . "\n";


#Create map file
my $no_of_predicted_proteins  = 0;
my $no_of_matched_proteins  = 0;
my $no_of_duplicated_proteins  = 0;

foreach my $predicted_protein (sort (keys %h_predicted_proteins)) {

	 $no_of_predicted_proteins++;

	if (&check_for_empty_string ($h_e_value {$predicted_protein})) {

		print  MAP_FILE $predicted_protein . "\t" . "NO_MATCH" . "\n";
		$h_new_protein {$predicted_protein} = $predicted_protein;
		
		#Print the fasta sequence to a file
		print NEW_PROT_FILE ">" . $resource . "-" . $predicted_protein . "\n";
		print NEW_PROT_FILE $h_predicted_fasta {$predicted_protein} . "\n";
		
	} else {

		my $test_protein = $h_e_value {$predicted_protein};

		
		#check for duplicate match
		if (my $duplicate_protein = &check_for_duplicate_prediction ($test_protein)) {
		
			print  MAP_FILE $predicted_protein . "\t" . $test_protein . "\t" . "Duplicate:" . $duplicate_protein . "\n";
			$no_of_duplicated_proteins++;

		} else {

			$h_match {$test_protein} = $predicted_protein;
			print  MAP_FILE $predicted_protein . "\t" . $test_protein . "\n";
			$no_of_matched_proteins++;	

		}
		
	}


	

}

print  MAP_FILE "#\n#\n# Total number of predicted proteins encoded in chromosome $prefix: " .  $no_of_predicted_proteins . "\n";
print  MAP_FILE "# Total number of test proteins encoded in chromosome $prefix: " .  $no_of_test_proteins . "\n";
print  MAP_FILE "# Total number of matched proteins: " .  $no_of_matched_proteins . "\n";
print  MAP_FILE "# Total number of duplicated proteins: " .  $no_of_duplicated_proteins . "\n";

### Print to Prot_id file
my $perfect = 0;
my $partial = 0;
my $no_of_non_matched_UniProt_protiens = 0;

print  PROT_ID_FILE "#1:UniProt ID 2:Query coverage (100 = perfect) 3:Ident (100 = perfect) 4:Query ID\n";
print  PROT_ID_FILE "# PERFECT MATCH\n";
foreach my $key (sort (keys %h_all_0)) {
	print PROT_ID_FILE $h_all_0 {$key} . "\t" . $h_queries {$key} . "\t" . $key . "\n";
	$perfect++;
}

print  PROT_ID_FILE "# PARTIAL MATCH\n";
foreach my $key (sort (keys %h_partial_match)) {
	print  PROT_ID_FILE   $h_partial_match {$key} . "\t" . $h_queries {$key} . "\t" .  $key . "\n";
	$partial++;
}

print  PROT_ID_FILE "# NON-MATCHED UniProt proteins\n";
foreach my $key (sort (keys %h_test_proteins)) {

	if ((&check_for_empty_string ($h_test_protein_match {$key}) )&& ($h_protein_chr_map {$key} eq $prefix)) {
		print  PROT_ID_FILE $key . "\n";
		$no_of_non_matched_UniProt_protiens++;
	}
}

print  PROT_ID_FILE "# NEW Protein\n";
foreach my $key (sort (keys %h_new_protein)) {
	print PROT_ID_FILE  $resource . "-" . $h_new_protein {$key} . "\t" . $h_queries {$key} . "\n";
	$no_of_new_proteins++;
}

### protein summary file
my $percent_perfect = ($perfect/$no_of_predicted_proteins) * 100;
my $percent_partial = ($partial/$no_of_predicted_proteins) * 100;
my $percent_total_count = ($total_count/$no_of_predicted_proteins) * 100;
my $percent_no_of_duplicated_proteins = ($no_of_duplicated_proteins/$no_of_predicted_proteins) * 100;
my $percent_no_of_matched_proteins = ($no_of_matched_proteins/$no_of_predicted_proteins) * 100;



print PROTEIN_SUMMARY_FILE "$perfect,$partial,$total_count,$no_of_duplicated_proteins,$no_of_predicted_proteins,$no_of_test_proteins,$no_of_matched_proteins,";
print PROTEIN_SUMMARY_FILE "$no_of_new_proteins,$no_of_non_matched_UniProt_protiens,";

printf PROTEIN_SUMMARY_FILE ("%.1f,%.1f,%.1f,%.1f,%.1f\n",$percent_perfect,$percent_partial,$percent_total_count,$percent_no_of_duplicated_proteins ,$percent_no_of_matched_proteins); 



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

	my $test_protein = shift;
	
	foreach my $key (sort (keys %h_match)) {

		if ($key eq $test_protein) {

			return $h_match {$key};
		}

	}

	return;

	
}
