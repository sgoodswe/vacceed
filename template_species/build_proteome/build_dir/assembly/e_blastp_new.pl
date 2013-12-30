#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $prefix = $ARGV [0];
my $assembly_dir = $ARGV [1] . "/";
my $proteome_dir = $ARGV [2] . "/";
my $prot_id_prefix =  $ARGV [3];


## hard coding
my $ident_threshold = 0.95;
my $l_query_coverage_threshold = 0.75;
my $u_query_coverage_threshold = 1.25;


#input files
my $input_file = $assembly_dir . $prefix . "_merged_blastp.txt";
my $i_fasta_file =  $assembly_dir . $prefix . "_merged_protein.fasta";
my $i_scored_genes = $assembly_dir . $prefix . "_scored_genes.txt";


#output file
my $output_file = $assembly_dir . $prefix . "_scored_proteins.txt";
my $new_seq_file = $proteome_dir . "new_prot.fasta";

## Hard coding for debugging ##
#$prefix = "chrIa";
#$input_file = "chrIa_merged_blastp.txt";
#$i_fasta_file =  "chrIa_merged_protein.fasta";
#$output_file =  "chrIa_scored_proteins.txt";
#$i_scored_genes = "chrIa_scored_genes.txt";


#Global variables
my $cluster_index = 0;
my $id;
my $score = 1;
my $new_seq_index = 0;


#create hashes
my %h_proteins;
my %h_proteins_lengths;
my %h_gene_scores;
my %h_bitscores;
my %h_sorted_query_file;
my %h_cluster;
my %h_cluster_score;


#Open debug file
#open DEBUG_FILE,'>debug.txt';

#File handling
open FASTA_FILE,$i_fasta_file ;
open SCORED_GENES_FILE,$i_scored_genes ;
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' .  $output_file;

#check if new prot fasta exists
if (-e $new_seq_file) {
	open NEW_FASTA_FILE,'>>' .  $new_seq_file;
} else {
	open NEW_FASTA_FILE,'>' .  $new_seq_file;
}

#LOOP  through protein FASTA file
while (defined (my $line = <FASTA_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	#Check for > anf get ID
	if ($line =~ m/^>(.+)/) {
		$id = $1;
	} else {
		$h_proteins {$id} = $line;
		
		$h_proteins_lengths {$id} = length ($line);
	
	}
			
} #end of LINE:while

#LOOP  through scored genes file
while (defined (my $line = <SCORED_GENES_FILE>)) {

	#Check for # 
	if ($line =~ m/^#/) { next;}
	
	$line = &trim ($line);
	
	$line =~ m/^(\d+),(.+)/;
	
	$h_gene_scores {$1} = $2;
			
} #end of LINE:while


#LOOP through the blast file and extract the highest bit score for each query _ subject
while (defined (my $line = <INPUT_FILE>)) {

	$line = &trim ($line);

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
	
	#get the query resource program
	$subject_id =~ m/(\w+)-/;
	my $subject_resource = $1;
	
			
	## Do not add to cluster if the blast hit is from the same resource
	if ($query_resource ne $subject_resource) {
	
	
			#Query length - subject length = query covergae
			my $query_coverage = (($subject_length - $query_length)/$subject_length);
			$query_coverage = sprintf("%.2f",1 - $query_coverage);
			
				#score if query if coverage is greater than subject gene
				if ($query_coverage > 1) {
					$query_coverage = $query_coverage - 1;
					$query_coverage = 1 - $query_coverage;
					
					#check if query coverage  < 0  (i.e. negative)
					if ($query_coverage < 0) {
						$query_coverage = 0;
					}
				}
				
										
				#convert ident to value less than 1
				$ident = $ident / 100;
			
			
			
			if (($query_coverage >= $l_query_coverage_threshold) && ($query_coverage <= $u_query_coverage_threshold) && ($ident >= $ident_threshold)) {
			#if (($query_coverage >= $coverage_threshold) && ($ident >= $ident_threshold))  {
				$score = ($query_coverage * $ident);
				&add_to_cluster ($query_id,$subject_id,$score);
			}
			
	}
}

## Print out the clusters

my $top_score = 1;
my $score_factor = 1;
my $protein_score = 0;



print OUTPUT_FILE "#1=New seq ID,2=cluster ID,3 to ?=protein identifiers in cluster,3rd from last=accumulative score for protein cluster\n";
print OUTPUT_FILE "#2nd last=score for gene clusters,last=score for new protein existence\n";

foreach my $cluster_id (sort { $h_cluster_score {$b} <=> $h_cluster_score {$a}} keys  %h_cluster_score ) {

	#print new sequence and return new sequence ID
	my $new_seq_ID = &print_new_seq ($h_cluster{$cluster_id});


	#get gene score if available
	my $gene_score = &get_gene_score ($h_cluster{$cluster_id});

	if ($top_score == 1) {
			$score_factor = (1.0 / $h_cluster_score{$cluster_id});
			$top_score = 0;
			$protein_score = 1;
		} else {

			$protein_score = $h_cluster_score{$cluster_id} * $score_factor;
		}

	$protein_score = $protein_score * $gene_score;
	
	$protein_score = sprintf ("%.2f",$protein_score);


	print OUTPUT_FILE $new_seq_ID . "," . $cluster_id . ",". $h_cluster{$cluster_id} . "," . $h_cluster_score{$cluster_id} . "," . $gene_score . ",";	 	 print OUTPUT_FILE $protein_score ."\n";
	
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
	my ($query_id, $cluster_proteins) = (@_);
	
	#get the query resource program
	$query_id =~ m/(\w+)-/;
	my $query_resource = $1;
		
	my @ids = split (',',$cluster_proteins);
	
	foreach my $id (@ids) {
	
		#get the cluster gene resource program
		$id =~ m/(\w+)-/;
		my $cluster_protein_resource = $1;
		
		#only do the following if they are from the same resource
		if ($query_resource eq $cluster_protein_resource) {
				if ($h_proteins {$id} eq $h_proteins {$query_id}) {
					return 0;
			}
		}
	}

	return 1;
}

#print new sequences
sub  print_new_seq ($)
{
	my ($cluster_proteins) = (@_);
	my $length = 0;
	my $new_id;
			
	my @ids = split (',',$cluster_proteins);
	
	#get the sequence with the longest length
	foreach my $id (@ids) {
	
		if ($h_proteins_lengths {$id} > $length) {

			$length = $h_proteins_lengths {$id};
		    $new_id = $id;
		}
		
	}

	$new_seq_index++;
	my $new_seq_ID = $prefix . "_new" . $new_seq_index;

	print NEW_FASTA_FILE ">" . $prot_id_prefix . "|" . $new_seq_ID . "|" . $cluster_proteins . "\n";
	print NEW_FASTA_FILE $h_proteins {$new_id} . "\n";

	#return new sequence ID
	return $new_seq_ID;
}

sub  get_gene_score ($)
{
	my ($cluster_proteins) = (@_);
	my $score = 0;
				
	my @ids = split (',',$cluster_proteins);
	
	#Find matching id
	foreach my $id (@ids) {
	
		foreach my $gene_cluster_id (sort keys %h_gene_scores) {
		
			if ($h_gene_scores{$gene_cluster_id} =~ m/\b$id\b/) {
			
				#get the last value from $h_gene_scores cluster which is the score
				my @gene_ids = split (',',$h_gene_scores{$gene_cluster_id});
				
				foreach (@gene_ids) {
						$score = $_;
				}
				return $score;
			}
		}
	}
	
	return 0;
}
