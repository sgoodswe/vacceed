#!/usr/bin/perl
use strict;
#use warnings;

#Program: e_blast_sim.pl


## hard coding
my $l_query_coverage_threshold = 25;
my $u_query_coverage_threshold = 175;
my $ident_threshold = 25;


#get the command line arguments
my $prefix = $ARGV [0];
my $prot_dir =  $ARGV [1] . "/";
my $prog_dir = $ARGV [2] . "/";
my $output_dir = $ARGV [3] . "/";
my $assembly_dir = $ARGV [4] . "/";
my $map_dir =  $ARGV [5] . "/";
my $resource = $ARGV [6];
my $blast_file = $ARGV [7];
my $fasta_file = $ARGV [8];
my $map_protein_gene =  $ARGV [9];
my $map_extra = $ARGV [10];
## Warning the following will be different with a different test fasta file
## The fasta header maybe for example >tr|UniProt_ID OR >gi|NCBI_ID
my $subject_prefix_id = $ARGV [11];

### mapping files should be in mapping directory
my $i_chr_map =  $map_dir . $map_protein_gene;
my $i_map_uniprot_gi = $map_dir . $map_extra;

my $input_file =  $blast_file;
my $i_fasta_file = $fasta_file;

my $o_summary_file = $output_dir . $prefix .  "_sim_summary.txt";
my $o_prot_id_file = $assembly_dir .  $prefix . "_" . $resource . "_sim_IDs.txt";
my $o_new_prot_file =$assembly_dir . $prefix . "_" .  $resource . "_new_sim.fasta";

## DEBUG
#$prefix = "chrIa";
#$resource= "gl";
#$input_file = "chrIa_homology_blastp.txt";
#$i_fasta_file = $prefix . ".aa";
#$o_summary_file = $prefix .  "_sim_summary.txt";
#$o_prot_id_file = $prefix . "_" . $resource . "_sim_IDs.txt";
#$i_chr_map =  "map_uniprot_gene_chr.txt";
#$i_map_uniprot_gi = "map_uniprot_gi.txt";
#$o_new_prot_file = $prefix . "_" .  $resource . "_new_sim.fasta";



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
	
	if ($line =~ m/^>(.+)/) {
	
		$id = $1;
		$h_predicted_proteins {$id} = $1;

	} else {
	
		$h_predicted_fasta {$id} = $h_predicted_fasta {$id} . $line;

	}

}

#File handling
open INPUT_FILE,$input_file ;
open SUMMARY_FILE,'>'  . $o_summary_file;
open PROT_ID_FILE,'>'  . $o_prot_id_file;
open NEW_PROT_FILE,'>'  . $o_new_prot_file;


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
				$subject_id = "NO_MAP_$subject_id";
			} else {
				$subject_id = $h_gi_uniprot_map {$subject_id};
			}
	}
			
	if ($evalue == 0) {
	
		#ALL	
		if (($query_length == $subject_length) && ($ident == 100))  {

			$h_all_0 {$query_id} = $subject_id;
					
		} elsif ($query_length == $subject_length) {
	
			$h_query_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
			
		} elsif ($ident == 100){

			$h_ident_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
		
			
		#Partial match if query coverage > than lower threshold and < upper threshold AND > ident threshold
		} elsif (($query_coverage > $l_query_coverage_threshold) && ($query_coverage < $u_query_coverage_threshold) && ($ident > $ident_threshold)) {

			$h_partial_coverage_0 {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
	
		} else {
		
			$h_non_match {$query_id} = $subject_id;
		}

	#Matches an existing protein  but e-value not 0
	} else {
	
		#ALL	
		if (($query_length == $subject_length) && ($ident == 100))  {

			$h_all {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;;
		
		} elsif ($query_length == $subject_length) {

			$h_query_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
	
		} elsif ($ident == 100){

			$h_ident {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
				
		#Partial match if query coverage > than lower threshold and < upper threshold AND > ident threshold
		} elsif (($query_coverage > $l_query_coverage_threshold) && ($query_coverage < $u_query_coverage_threshold) && ($ident > $ident_threshold)) {

			$h_partial_coverage {$query_id} = $subject_id;
			$h_partial_match {$query_id} = $subject_id;
	
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


### Print to Prot_id file

print  PROT_ID_FILE "#1:UniProt ID 2:Query coverage (100 = perfect) 3:Ident (100 = perfect) 4:Query ID\n";
print  PROT_ID_FILE "# PERFECT MATCH\n";
foreach my $key (sort (keys %h_all_0)) {
	print PROT_ID_FILE $h_all_0 {$key} . "\t" . $h_queries {$key} . "\t" . $key . "\n";
}

print  PROT_ID_FILE "# PARTIAL MATCH\n";
foreach my $key (sort (keys %h_partial_match)) {
	print  PROT_ID_FILE   $h_partial_match {$key} . "\t" . $h_queries {$key} . "\t" .  $key . "\n";
	
	if ($h_partial_match {$key} =~ m/NO_MAP_(.+)/) {
	
		my $subject_id = $1;
	
		#Print the fasta sequence to a file
		print NEW_PROT_FILE ">" . $key . '|' . $subject_id . "\n";
		print NEW_PROT_FILE $h_predicted_fasta {$key} . "\n";
	}
}


print  PROT_ID_FILE "# NEW Protein\n";
foreach my $key (sort (keys %h_non_match)) {
	print PROT_ID_FILE  $key . "\t" . $h_queries {$key} . "\t" . $h_non_match {$key} . "\n";
	
	my $subject_id;
	
	#get subject id
	if ($h_non_match {$key} =~ m/NO_MAP_(.+)/) {
		$subject_id = $1;
	} else {
		$subject_id = $h_non_match {$key};
	}
	
	
	#Print the fasta sequence to a file
	print NEW_PROT_FILE ">" . $key . '|' . $subject_id . "\n";
	print NEW_PROT_FILE $h_predicted_fasta {$key} . "\n";
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
