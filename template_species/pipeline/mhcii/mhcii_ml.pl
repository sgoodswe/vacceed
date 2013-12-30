#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $input_file = $ARGV [0];
my $output_dir = $ARGV [1] . "/";
my $train_dir = $ARGV [2] . "/";


## hard coding
my $train_file = $train_dir . "generic_training_set";
my $output_file = $output_dir . "mhcii_ml.txt";
my $threshold = 50;
my $ignore_CombLib = 0;
my $ignore_Sturniolo = 0;



#Global variables
my $ID;
my $start = 1;
my $get_headers = 1;
my $get_scores = 0;
my $add_required = 0;
my $min_Percentile = 999999;
my $min_percentile_scores;
my $no_of_alleles = 0;
my $no_of_cols = 0;
my $no_of_ids = 0;
my $header;


#create hashes
my %h_score;
my %h_candidate_status;
my %h_header_check;



#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open TRAIN_FILE,$train_file ;
open OUTPUT_FILE,'>'  . $output_file;




#LOOP through protein evidence set (the map file)
while (defined (my $line = <TRAIN_FILE>)) {

	$line = &trim ($line);

	my ($ID,$Candidate) = split (',',$line);

	$h_candidate_status {$ID} = $Candidate;

}


#LOOP 
while (defined (my $line = <INPUT_FILE>)) {
		
	$line = &trim ($line);
	
	#Check for # 
	if ($line =~ m/^#/) { 
	
			next;
			
	} elsif (($line =~ m/Allele/) || ($line =~ m/netMHCIIpan/)) { 
	
			next;
	#Get sequence ID		
	} elsif ($line =~ m/^>..\|(\w+)/) { 

	
		if ($start == 0) {
			$start = 1;
	
			if ($add_required == 1) {
				#Add average score to hash
				&add_average_score;
			}
	
			$get_headers = 0;
			$get_scores = 0;
			%h_header_check = ();
							
		} else {
			$start = 1;
		}
		
		$ID = $1;
		#print DEBUG_FILE $no_of_alleles . "\t" . $no_of_cols . "\n";
		$no_of_cols = 0;
		
	
	# Get binding info				
	} elsif ($line =~ m/^BINDING_TEST\s+(.+)/) {
	
		
		#check if the allele has already been used i.e. there may be duplicate alleles in the input file
		if (&check_for_empty_string ( $h_header_check {$1})) {
				$h_header_check {$1} = $1;
				$get_scores = 1;
		} else {
							
				$get_scores = 0;
				next;
		}


		if ($start == 0) {
		
			#Add average score to hash
			&add_average_score;
					
		} else {
			$start = 0;
		}
		
		if ($get_headers == 1) {
					
			$header = $header . "," . $1;
					
			$no_of_alleles++;
		}
		
		#reset 
		$min_Percentile = 9999999;
		
			
	# Ignore line if no matching species				
	} elsif ($line =~ m/Could not find tools matching species/) {
	
		next;
	
	} elsif ($get_scores == 1) {
	
			
		my ($MHC,$seq_no,$start,$end,$Method_used,$Peptide,$Percentile_Rank,$comblib_core,$CombLib_IC50,$comblib_percentile,$smm_core,$SMM_IC50,$smm_percentile,$nn_core,$ANN_IC50,$nn_percentile,				$netMHCIIpan_core,$NetMHCpan_IC50,$netMHCIIpan_percentile,$Sturniolo_core,$Sturniolo_IC50,$Sturniolo_percentile) = split ('\t',$line);
			
		
		#Get minimum percentile
		if ($Percentile_Rank < $min_Percentile) {
			$min_Percentile = $Percentile_Rank;
			$min_percentile_scores = $ANN_IC50 . "," . $SMM_IC50 . "," . $CombLib_IC50 . "," . $NetMHCpan_IC50 . "," . $Sturniolo_IC50;
		}

		$add_required = 1;	
			
	}
		
} #end of LINE:while



#get the average score for the last allele
if ($add_required == 1) {
	#Add average score to hash
	&add_average_score;
}



#Print out
print OUTPUT_FILE "ID" . $header . ",Candidate\n";


foreach my $id ( sort keys %h_score ) {


	if (&check_for_empty_string ( $h_candidate_status {$id})) {

		$h_candidate_status {$id} = 0;
	}

	print OUTPUT_FILE $id . $h_score {$id} . "," . $h_candidate_status {$id} . "\n";
	$no_of_ids++;

}


print "# Number of proteins = $no_of_ids\n";
print "# Number of allele combinations = $no_of_alleles\n";




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

sub add_average_score ()
{
	#get the scores for the min percentile
	my ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score,$Sturniolo_score) = split (',', $min_percentile_scores);
	
	#print OUTPUT_FILE "Check_in_abs " . $ID . "\t" .  $ANN_score . "\t" . $SMM_score . "\t" .	$CombLib_score . "\t" .	$NetMHCpan_score. "\t" . $binding_info  . "\n";
	
	#Ignore CombLib score
	if ($ignore_CombLib == 1) {
		$CombLib_score = "-";
	}
	
	#Ignore CombLib score
	if ($ignore_Sturniolo == 1) {
		$Sturniolo_score = "-";
	}
					
	my @scores = ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score,$Sturniolo_score);
	
	my $sum_of_scores = 0;
	my $no_of_methods = 0;
			
	foreach my $score (@scores) {
	
		#check for a dash
		if ($score ne "-" ) {
			$sum_of_scores = $sum_of_scores + $score;
			$no_of_methods++;
		}
	}
	
	my $av_score;
	
	if ($sum_of_scores != 0) {
		$av_score = ($sum_of_scores/$no_of_methods);
		
		#reduce number of decimal places
		$av_score = sprintf ("%.3f",$av_score);
	} else {
		$av_score = "9999";
	}
	
	$h_score {$ID} =  $h_score {$ID} . "," . $av_score;
	$no_of_cols++;
	
	$add_required = 0;
}
