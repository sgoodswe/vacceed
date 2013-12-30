#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $proteome_dir = $ARGV [0] . "/";
my $assembly_dir = $ARGV [1];
my $map_dir =  $ARGV [2] . "/";
my $map_file = $ARGV [3];
my $prot_dir =  $ARGV [4] . "/";
my $prot_file =  $ARGV [5];
my $prefix =  $ARGV [6];

## hard coding
my $output_file = $proteome_dir . "proteome_info.txt";
$prefix = $prefix . '\|';

## DEBUG
#my $output_file =  "proteome_info.txt";
#my $protein_file = "UniProt_proteins.fasta";
#my $map_file = "map_gene_uniprot.txt";



my $protein_file = $prot_dir . $prot_file;
my $map_file = $map_dir . $map_file;


#Global variables
my $header;
my $no_of_evidence;


#create hashes
my %h_proteins;
my %h_evidence;
my %h_assemble;
my %h_map;
my %h_total_score;

#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open PROT_FILE, $protein_file;
open MAP_FILE, $map_file;
open OUTPUT_FILE,'>' . $output_file;


#LOOP  through protein file
while (defined (my $line = <PROT_FILE>)) {

	#Get the ID
	if ($line =~ m/>$prefix(\w+)/) { 
		$h_total_score {$1} = 0;
	}
			
} #end of LINE:while

#LOOP  through map file
while (defined (my $line = <MAP_FILE>)) {

	$line = &trim ($line);
	
	my ($gene_id,$uniprot_id) = split ('\s+',$line);
	$h_map {$gene_id} = $uniprot_id;
			
} #end of LINE:while


#Get files from assmebly directory
opendir(DIR, $assembly_dir );
my @FILES= readdir(DIR); 

foreach my $file (@FILES) {

	#split file name into components
	my ($chr,$resource,$type,$ext) = split ('_',$file);

	#process gene types or prot types
	if (($type eq "gene") || ($type eq "prot") || ($type eq "sim")  ){
	
			&process_resource ($file);
	}
}

## Assemble the evidence
foreach my $resource( sort keys %h_evidence ) {

	#build the header
	if (!$header) {
		$header = "#ID," . $resource;
	} else {
		$header = $header . "," . $resource;
	}

	#record the number of evidence items
	$no_of_evidence++;
	
	foreach my $id (sort keys   %{ $h_evidence {$resource}}) {

			if (&check_for_empty_string ($h_assemble {$id})) {
				$h_assemble {$id} =  $id . "," . $h_evidence {$resource} {$id};
			} else {
				$h_assemble {$id} = $h_assemble {$id} . "," . $h_evidence {$resource} {$id} ;
			}
	
	}
	
}



## Get the protein existence scores for new proteins
foreach my $file (@FILES) {

	#split file name into components
	my ($chr,$resource,$type,$ext) = split ('_',$file);

	#process only resource = "scored" and type = "proteins"
	if (($resource eq "scored") && ($type eq "proteins.txt")) {
	
		#open the file to read
		open INPUT_FILE,$assembly_dir . "/" . $file ;

			#Loop through file to read
			while (defined (my $line = <INPUT_FILE>)) {

				my $id;
				my $score;

				#Check for # 
				if ($line =~ m/^#/) { next;}

				$line = &trim ($line);

				#get first and last column
				my @columns = split (',',$line);

				foreach (@columns) {

					if (&check_for_empty_string ($id)) {

						$id = $_;
					} else {
						$score = $_;
			
					}	
		
				}

				my $line_to_add;
				
				my @headers = split (',',$header);

				foreach (@headers) {

					if (&check_for_empty_string ($line_to_add)) {

						$line_to_add = $id;

					} else {
				 		$line_to_add = $line_to_add . ",0.00"
					}

				}

				$h_assemble {$id} =  $line_to_add;

				#multiple score by number of evidence items (i.e. number of headers). This will
				# standardise the score since the average is calculated later
				$h_total_score {$id} = $score * $no_of_evidence;
				
		}
	}
}




#print the proteome information
print OUTPUT_FILE $header . ",Average score\n";

foreach my $id (sort { $h_total_score {$b} <=> $h_total_score {$a}} keys   %h_total_score) {


	print OUTPUT_FILE  $h_assemble{$id} . ",";

	my $average_score = $h_total_score {$id} / $no_of_evidence; 
	printf OUTPUT_FILE ("%.2f\n",$average_score);
	

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

sub  process_resource ($)
{
		my $file = shift;
		my $match_type;
		my $perfect_score = "1.00";
		my $partial_score;
			
		#split file name into components
		my ($chr,$resource,$type,$ext) = split ('_',$file);
		
		#initial the hash evidence for each protein 
		foreach my $id (sort keys %h_total_score) {

			#only do this if it is empty
			if (&check_for_empty_string ($h_evidence {$resource . "_" . $type} {$id})) {
				$h_evidence {$resource . "_" . $type} {$id} = "0.00";
			}
		}
	
		#open the file to read
		open INPUT_FILE,$assembly_dir . "/" . $file ;

			
		#Loop through file to read
		while (defined (my $line = <INPUT_FILE>)) {

			$line = &trim ($line);
			
			if ($line =~ m/PERFECT MATCH/) {
				$match_type = "perfect";
			}elsif ($line =~ m/PARTIAL MATCH/) {
				$match_type = "partial";
			}elsif ($line =~ m/#/) {
				$match_type = "";
			} else {
			
				#If match is either perfect or partial
				if ($match_type) {
				
					#expected format = #1:UniProt ID 2:Query coverage (100 = perfect) 3:Ident (100 = perfect) 4:Query ID
				
					my ($id,$query_coverage,$ident,$query_id) = split ('\s+',$line);
					
					#convert gene ID to UniProt ID
					if ($type eq "gene") {
						$id = $h_map {$id};
					}
					
					#Only do this once - ignore duplicates especially for est evidence
					if 	($h_evidence {$resource . "_" . $type} {$id} eq "0.00") {
					
						if ($match_type eq "perfect") {
							$h_evidence {$resource . "_" . $type} {$id}  = $perfect_score;
						} else {
						
							#score if query coverage is greater than test gene
							if ($query_coverage > 100) {
								$query_coverage = $query_coverage - 100;
								$query_coverage = 100 - $query_coverage;
								
								#check if query coverage  < 0  (i.e. negative)
								if ($query_coverage < 0) {
									$query_coverage = 0;
								}
							}
							#convert query coverage to value less than 1
							$query_coverage = $query_coverage / 100;
													
							#convert ident to value less than 1
							$ident = $ident / 100;
							
							#calculate partial score
							$partial_score = $query_coverage * $ident;
							
							#Convert to 2 decimal places
							$partial_score  = sprintf ("%.2f",$partial_score);
							
													
							$h_evidence {$resource . "_" . $type} {$id}  = $partial_score;
						}
				
						$h_total_score {$id} = $h_total_score {$id} + $h_evidence {$resource . "_" . $type} {$id};
					
						
					}
				}
			}
		} #end of while loop for reading file

}
