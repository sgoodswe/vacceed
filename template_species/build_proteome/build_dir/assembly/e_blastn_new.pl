#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $prefix = $ARGV [0];
my $assembly_dir = $ARGV [1] . "/";
my @EST = split (',',$ARGV [2]);

## hard coding
my $ident_threshold = 75;

#input files
my $input_file = $assembly_dir . $prefix . "_merged_blastn.txt";
my $i_fasta_file =  $assembly_dir . $prefix . "_merged_gene.fasta";

#output file
my $output_file = $assembly_dir . $prefix . "_scored_genes.txt";



## Hard coding for debugging ##
#@EST = split (',',"blat,gmap");
#$prefix = "chrIa";
#$input_file = "chrIa_merged_blastn.txt";
#$i_fasta_file =  "chrIa_merged_gene.fasta";
#$output_file =  "chrIa_scored_genes.txt";


#Global variables
my $cluster_index = 0;
my $id;
my $score = 1;


#create hashes
my %h_genes;
my %h_bitscores;
my %h_sorted_query_file;
my %h_cluster;
my %h_cluster_score;


#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open FASTA_FILE,$i_fasta_file ;
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' .  $output_file;

#LOOP 
while (defined (my $line = <FASTA_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	#Check for > and get ID
	if ($line =~ m/^>(.+)/) {
		$id = $1;
	} else {
		$h_genes {$id} = $line;
	}

} #end of LINE:while


#LOOP through the blast file and extract the highest bit score for each query _ subject
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);

	print DEBUG_FILE $line . "\n";

	#0.0,gmap-gene_1,547,gmap-gene_1,547,1011,547,100.00,547,547
	
	#0.0,gmap-gene_1,547,gmap-gene_111,569,1011,547,100.00,547,547
	
	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $line);

	if ($query_id ne $subject_id) {
	
		my $id = $query_id . $subject_id;
	
		#Get the highest bit score 
		if (&check_for_empty_string ($h_bitscores {$id})) {
				$h_sorted_query_file {$id} = $line;
				$h_bitscores {$id} = $bitscore;
		} elsif ($bitscore > $h_bitscores {$id}) {
				$h_sorted_query_file {$id} = $line;
				$h_bitscores {$id} = $bitscore;
		}
	}
			
} #end of LINE:while


#LOOP through query hash and check for full genes only
foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});

	my $EST = 0;
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
	
	my $query_EST = (&check_for_EST_program ($query_resource));
	
	#get the query resource program
	$subject_id =~ m/(\w+)-/;
	my $subject_resource = $1;
	
	my $subject_EST = (&check_for_EST_program ($subject_resource));
	
		
	## Do not add to cluster if the blast hit is from the same resource
	if ($query_resource ne $subject_resource) {
		
		if ( ($query_EST == 0) && ($subject_EST == 0)) {
			if (($query_length == $subject_length) && ($ident >= $ident_threshold))  {
				$score = 1 * ($ident/100);
				&add_to_cluster ($query_id,$subject_id,$score);
			}
		}
		
	}
}

#LOOP through query hash for Query = full gene and Subject = EST
foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});

	my $EST = 0;
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
	
	my $query_EST = (&check_for_EST_program ($query_resource));
	
	#get the query resource program
	$subject_id =~ m/(\w+)-/;
	my $subject_resource = $1;
	
	my $subject_EST = (&check_for_EST_program ($subject_resource));
			
	if (($query_EST == 0) && ($subject_EST == 1) && ($ident == 100)) {
		$score = 1;
		&add_to_cluster ($query_id,$subject_id,$score);
	}
}

#LOOP through query hash for Query = EST and Subject = FULL
foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});

	my $EST = 0;
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
	
	my $query_EST = (&check_for_EST_program ($query_resource));
	
	#get the query resource program
	$subject_id =~ m/(\w+)-/;
	my $subject_resource = $1;
	
	my $subject_EST = (&check_for_EST_program ($subject_resource));
			
	if (($query_EST == 1) && ($subject_EST == 0) && ($ident == 100)) {
		$score = 1;
		&add_to_cluster ($query_id,$subject_id,$score);
	}
}

#LOOP through query hash for Query = EST and Subject = EST
foreach my $query_search (sort { $a <=> $b } keys %h_sorted_query_file) {

	my ($evalue, $query_id, $query_length, $subject_id, $subject_length, $bitscore, $score,$ident,$mismatches, $matches) = split (',', $h_sorted_query_file {$query_search});

	my $EST = 0;
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
	
	my $query_EST = (&check_for_EST_program ($query_resource));
	
	#get the query resource program
	$subject_id =~ m/(\w+)-/;
	my $subject_resource = $1;
	
	my $subject_EST = (&check_for_EST_program ($subject_resource));
			
	if (($query_EST == 1) && ($subject_EST == 1)) {
	
		##  ESTs (100% coverage)		
		if (($query_length - $matches) == 0) {
			$score = 1;
			&add_to_cluster ($query_id,$subject_id,$score);
		}
		
	}
}

my $top_score = 1;
my $score_factor = 1;
my $gene_score = 0;

print OUTPUT_FILE "#1=cluster ID,2 to ?gene identifiers in cluster,2nd from last=number of elements in score,last=score for gene cluster\n";

## Print out the clusters
foreach my $cluster_id (sort { $h_cluster_score {$b} <=> $h_cluster_score {$a}} keys  %h_cluster_score ) {


	if ($top_score == 1) {
		$score_factor = (1.0 / $h_cluster_score{$cluster_id});
		$top_score = 0;
		$gene_score = 1;
	} else {

		$gene_score = $h_cluster_score{$cluster_id} * $score_factor;
	}

	$gene_score = sprintf ("%.2f",$gene_score);

	print OUTPUT_FILE $cluster_id . ","  . $h_cluster{$cluster_id} . "," . $h_cluster_score{$cluster_id} . "," .  $gene_score  . "\n";
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

#Check for EST program
sub  check_for_EST_program ($)
{
	my $program = shift;
	
	foreach my $EST_program (@EST) {
	
		if ($EST_program eq $program) {
			return 1;
		} 
	}
	return 0;
}

#add to cluster
sub  add_to_cluster ($)
{
	my ($query_id, $subject_id,$score) = (@_);

	my $increment_cluster = 1;
	
	foreach my $cluster_id (sort { $a <=> $b } keys %h_cluster ) { 
	
		if (($h_cluster{$cluster_id} =~ m/\b$query_id\b/) &&  ($h_cluster{$cluster_id} !~ m/\b$subject_id\b/)){
		
				if (&check_seq_match ($subject_id,$h_cluster{$cluster_id})) {
					$h_cluster{$cluster_id} = $h_cluster{$cluster_id} . "," . $subject_id;
					$h_cluster_score{$cluster_id} = $h_cluster_score{$cluster_id} + $score;
				}
				$increment_cluster = 0;
				last;
		}elsif (($h_cluster{$cluster_id} !~ m/\b$query_id\b/) &&  ($h_cluster{$cluster_id} =~ m/\b$subject_id\b/)) {
		
				if (&check_seq_match ($query_id,$h_cluster{$cluster_id})) {
					$h_cluster{$cluster_id} = $h_cluster{$cluster_id} . "," . $query_id;
					$h_cluster_score{$cluster_id} = $h_cluster_score{$cluster_id} + $score;
				}
				$increment_cluster = 0;
				last;
		}elsif (($h_cluster{$cluster_id} =~ m/\b$query_id\b/) &&  ($h_cluster{$cluster_id} =~ m/\b$subject_id\b/)) {
				$increment_cluster = 0;
				last;
		}
	}
	
	if ($increment_cluster == 1) {
		$cluster_index++;
		$h_cluster{$cluster_index} = $query_id . "," . $subject_id;
		$h_cluster_score{$cluster_index} = $h_cluster_score{$cluster_index} + $score;
	}
	
	return 0;
}



#Check for matching sequence
sub  check_seq_match ($)
{
	my ($query_id, $cluster_genes) = (@_);
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
		
	my @ids = split (',',$cluster_genes);
	
	foreach $id (@ids) {
	
		#get the cluster gene resource program
		$id =~ m/(\w+)-/;
		my $cluster_gene_resource = $1;
		
		#only do the following if they are from the same resource
		if ($query_resource eq $cluster_gene_resource) {
				if ($h_genes {$id} eq $h_genes {$query_id}) {
					return 0;
			}
		}
	}

	return 1;
}
