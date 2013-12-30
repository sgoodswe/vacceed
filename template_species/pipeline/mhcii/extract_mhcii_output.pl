#!/usr/bin/perl
use strict;
#use warnings;


#get the command line arguments
my $input_file = $ARGV [0];
my $output_dir = $ARGV [1] . "/";


#hard coded output files
my $output_file = $output_dir . "extract.txt"; 
my $stats_file = $output_dir . "global_stats.txt";
my $peptide_file = $output_dir . "peptide_stats.txt";
my $id_peptide_file = $output_dir . "id_peptide_stats.txt";
my $id_peptide_allele_file = $output_dir . "id_peptide_allele.txt";
my $threshold = 50;
my $ignore_CombLib = 0;


#Global variables
my $ID;
my $binding_info;
my $peptide_used;
my $start = 1;
my $min_Percentile = 999999;
my $ANN_score = 9999999;
my $SMM_score = 9999999;
my $CombLib_score = 9999999;
my $NetMHCpan_score =9999999;
my $Sturniolo_score =9999999;
my $abs_min_score = 999999999;
my $avg_min_score = 0;
my $min_percentile_scores;
my $min_percentile_peptides;
my $total_no_of_peptides = 0;
my $no_of_mhc_alleles = 0;

#create hashes
my %h_ANN_min;
my %h_SMM_min;
my %h_CombLib_min;
my %h_NetMHCpan_min;
my %h_Sturniolo_min;
my %h_avg_min_score;

my %h_ANN_peptide;
my %h_SMM_peptide;
my %h_CombLib_peptide;
my %h_NetMHCpan_peptide;
my %h_Sturniolo_peptide;

my %h_no_of_peptides;

my %h_global_no_of_peptides;
my %h_not_used;

my %h_prom_peptides;
my %h_ID_peptide;
my %h_ID_peptide_allele;

#Open debug file
#open DEBUG_FILE,'>debug.txt';



#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;
open STATS_FILE,'>' . $stats_file;
open PEPTIDE_FILE,'>' . $peptide_file;
open ID_PEPTIDE_FILE,'>' . $id_peptide_file;
open ID_PEPTIDE_ALLELE_FILE,'>' . $id_peptide_allele_file;

print OUTPUT_FILE "#1=UniProt ID :2=Min score 3=score average:4=MHC2 allele :5=Method used :6=Peptide sequence :7=Allele with most binding sites\n";
print OUTPUT_FILE "#8=Max. number of binding sites  :9=Total number of binding sites for all alleles :10=Total number of MHC2 alleles available\n";
print OUTPUT_FILE "#11=No. of MHC2 alleles used :12=No. of MHC2 alleles not used\n"; 	

#LOOP through protein evidence set (the map file)
while (defined (my $line = <INPUT_FILE>)) {

#>tr|Q27298|Q27298_TOXGO P30 OS=Toxoplasma gondii GN=p30 PE=4 SV=1

	#Check for # 
	if ($line =~ m/^#/) { 
	
			next;
			
	} elsif (($line =~ m/Allele/) || ($line =~ m/netMHCIIpan/)) { 
	
			next;
	#Get sequence ID		
	} elsif ($line =~ m/>..\|(\w+)/) { 
	
		if ($start == 0) {
			&print_out;
			$start = 1;
		} else {
			$start = 1;
		}
		
		$ID = $1;
		$abs_min_score = 999999999;
		$avg_min_score = 0;
		$total_no_of_peptides = 0;
		$no_of_mhc_alleles = 0;
		
		#Reset
		%h_ANN_min = ();
		%h_SMM_min = ();
		%h_CombLib_min = ();
		%h_NetMHCpan_min = ();
		%h_Sturniolo_min = ();

		%h_ANN_peptide = ();
		%h_SMM_peptide = ();
		%h_CombLib_peptide = ();
		%h_NetMHCpan_peptide = ();
		%h_Sturniolo_peptide = ();
		%h_no_of_peptides = ();
		
	# Get binding info				
	} elsif ($line =~ m/BINDING_TEST\s+(.+)/) {
	
		if ($start == 0) {
		
			#get the abslotue minimum score
			&get_abs_min_score;
					
		} else {
			$start = 0;
		}
		
		$binding_info = $1;
		
		#reset the method scores
		$min_Percentile = 9999999;
		$no_of_mhc_alleles++;
		
	# Ignore line if no matching species				
	} elsif ($line =~ m/Could not find tools matching species/) {

		$h_not_used {$binding_info} = $h_not_used {$binding_info} + 1;
		
		next;
			
	} else {

		my ($MHC,$seq_no,$start,$end,$Method_used,$Peptide,$Percentile_Rank,$comblib_core,$CombLib_IC50,$comblib_percentile,$smm_core,$SMM_IC50,$smm_percentile,$nn_core,$ANN_IC50,$nn_percentile,		
		$netMHCIIpan_core,$NetMHCpan_IC50,$netMHCIIpan_percentile,$Sturniolo_core,$Sturniolo_IC50,$Sturniolo_percentile) = split ('\t',$line);

	
	#print DEBUG_FILE "#########################" . "\n";	
	#print DEBUG_FILE $line . "\n";
	#print DEBUG_FILE "allele " . $MHC . "\n";
	#print DEBUG_FILE "Sequence_no " . $seq_no . "\n"; 
	#print DEBUG_FILE "start " . $start . "\n"; 
	#print DEBUG_FILE "end " . $end . "\n"; 
	#print DEBUG_FILE "Mehthod_used " .$Method_used . "\n"; 
	#print DEBUG_FILE "peptide seq " .$Peptide . "\n"; 
	#print DEBUG_FILE "consensus_Rank " . $Percentile_Rank . "\n"; 
	#print DEBUG_FILE "comblib_core " . $comblib_core . "\n"; 
	#print DEBUG_FILE "comblib_score " . $CombLib_IC50 . "\n"; 
	#print DEBUG_FILE "comblib_percent " . $comblib_percentile . "\n"; 
	#print DEBUG_FILE "smm_core " . $smm_core . "\n"; 
	#print DEBUG_FILE "smm_score " . $SMM_IC50 . "\n"; 
	#print DEBUG_FILE "smm_percent " . $smm_percentile . "\n"; 
	#print DEBUG_FILE "nn_core " . $nn_core . "\n"; 
	#print DEBUG_FILE "nn_score " . $ANN_IC50 . "\n"; 
	#print DEBUG_FILE "nn_percent " . $nn_percentile . "\n"; 
	#print DEBUG_FILE "netMHCIIpan_core " . $netMHCIIpan_core . "\n"; 
	#print DEBUG_FILE "netMHCIIpan_score " . $NetMHCpan_IC50 . "\n"; 
	#print DEBUG_FILE "netMHCIIpan_prcent" . $netMHCIIpan_percentile . "\n"; 
	#print DEBUG_FILE "Sturniolo_core " . $Sturniolo_core . "\n"; 
	#print DEBUG_FILE "Sturniolo_score " . $Sturniolo_IC50 . "\n"; 
	#print DEBUG_FILE "Sturniolo_percent" . $Sturniolo_percentile . "\n"; 	
	
		
		#Assign zero if empty
		if (&check_for_empty_string ($h_no_of_peptides {$MHC})) {
			$h_no_of_peptides {$MHC} = 0;
		}
		
		#Assign zero if empty
		if (&check_for_empty_string ($h_global_no_of_peptides {$MHC})) {
			$h_global_no_of_peptides {$MHC} = 0;
		}
			
		if ($Method_used =~ m/Consensus\s+\((.+)\)/) {
		
			my @methods = split (',',$1);
			
			my $no_of_methods = scalar @methods;
			my $no_below_threshold = 0;
			
			foreach (@methods) {
			
				if (($_ eq "NN") &&  ($ANN_IC50 < $threshold)) {
					$no_below_threshold++;
				}elsif (($_ eq "SMM") && ($SMM_IC50 < $threshold)) {
					$no_below_threshold++;
				# check CombLib_Sidney2008
				}elsif (($_ =~ m/COMB/) && ($CombLib_IC50 < $threshold)) {
					$no_below_threshold++;
				}elsif (($_ =~ m/MHC/) && ($NetMHCpan_IC50 < $threshold)) {
					$no_below_threshold++;
				}elsif (($_ =~ m/Sturn/) && ($Sturniolo_IC50 < $threshold)) {
					$no_below_threshold++;
				}
				
				
			}
						
			#Check for consensus peptide
			if ($no_below_threshold ==  $no_of_methods) {
				
				$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
				$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
				$total_no_of_peptides++;
				
				$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
				$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
				
				$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
										
			}
		
		} elsif (($Method_used =~ m/MHC/) && ($NetMHCpan_IC50 ne "-")) {
				
				if ($NetMHCpan_IC50 < $threshold) {
					$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
					$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
				}
			
		} elsif (($Method_used eq "NN") && ($ANN_IC50 ne "-")) {
				
				if ($ANN_IC50 < $threshold) {
					$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
					$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
				}
			
		} elsif (($Method_used eq "SMM") && ($SMM_IC50 ne "-")) {
				
				if ($SMM_IC50 < $threshold) {
					$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
					$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
				}
			
		} elsif (($Method_used eq "Sturniolo") && ($Sturniolo_IC50 ne "-")) {
				
				if ($Sturniolo_IC50 < $threshold) {
					$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
					$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
				}
			
		} elsif (($Method_used =~ m/COMB/) && ($Sturniolo_IC50 ne "-")) {
				
				if ($CombLib_IC50 < $threshold) {
					$h_no_of_peptides {$MHC} = $h_no_of_peptides {$MHC} + 1;
					$h_global_no_of_peptides {$MHC} = $h_global_no_of_peptides {$MHC} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_peptide_allele {$ID} {$Peptide} {$MHC} = $h_ID_peptide_allele {$ID} {$Peptide} {$MHC} + 1;
				}
			
		}
		
		
		#Get minimum percentile
		if ($Percentile_Rank < $min_Percentile) {
			$min_Percentile = $Percentile_Rank;
			$min_percentile_scores = $ANN_IC50 . "," . $SMM_IC50 . "," . $CombLib_IC50 . "," . $NetMHCpan_IC50 . "," . $Sturniolo_IC50;
			$min_percentile_peptides = $nn_core . "," . $smm_core . "," . $comblib_core . "," . $netMHCIIpan_core . "," . $Sturniolo_core;
		}
	
	}

			
} #end of LINE:while

#print last prediction
&print_out;

my $allele_count = 0;

print STATS_FILE "### Frequency of MHC2 alleles used ###\n";
foreach my $key ( sort  { $h_global_no_of_peptides {$b} <=> $h_global_no_of_peptides {$a}} keys  %h_global_no_of_peptides) {

	print STATS_FILE  $key . "\t" . $h_global_no_of_peptides {$key} . "\n";
	
	$allele_count++;
	
}

print STATS_FILE "\nNumber of MHC2 alleles used = $allele_count\n";

$allele_count = 0;

print STATS_FILE "\n### No methods found for following MHC2 alleles ###\n";
#Print alleles not used
foreach my $key (sort { $h_not_used {$b} <=> $h_not_used {$a}} keys   %h_not_used) {

 print STATS_FILE $key . "\t" . $h_not_used {$key} . "\n";
 
 $allele_count++;
	
}

#Print alleles not used
print STATS_FILE "\nNumber of MHC2 alleles not used = $allele_count\n";


print ID_PEPTIDE_FILE "### Proteins with common peptides###\n";

foreach my $id ( sort keys %h_ID_peptide ) {

	print ID_PEPTIDE_FILE "ID:" . $id . "\n";
	
	my $child_hashref = $h_ID_peptide{$id};
   
		foreach my $peptide (sort { $child_hashref-> {$b} <=> $child_hashref->{$a}}   keys   %$child_hashref) {
	
		print ID_PEPTIDE_FILE $id . "\t" . $peptide . "\t" . $child_hashref-> {$peptide} . "\n";
	}
}

print PEPTIDE_FILE "### Common peptides###\n";

foreach my $key (sort { $h_prom_peptides{$b} <=> $h_prom_peptides{$a}} keys   %h_prom_peptides) {

 print PEPTIDE_FILE $key . "\t" . $h_prom_peptides {$key} . "\n";
 
	
}


foreach my $id ( sort keys %h_ID_peptide_allele ) {

	print ID_PEPTIDE_ALLELE_FILE "ID:" . $id . "\n";
	
	my $child_hashref = $h_ID_peptide{$id};
   
	foreach my $peptide (sort { $child_hashref-> {$b} <=> $child_hashref->{$a}}   keys   %$child_hashref) {
	
		print ID_PEPTIDE_ALLELE_FILE $id . "\t" . $peptide . "\t" . $child_hashref-> {$peptide} . "\n";
	
		foreach my $allele (sort keys   %{$h_ID_peptide_allele {$id} {$peptide}}) {
	
			print ID_PEPTIDE_ALLELE_FILE "\t" . $allele . "\t" . $h_ID_peptide_allele {$id} {$peptide} {$allele} . "\n";
			
		}
	}
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

sub print_out ()
{

	#get the abslotue minimum score
	&get_abs_min_score;
	
	#get  binding info
	my $method_used = &get_binding_info;
	
	#get peptide
	my $peptide_used = &get_peptide;
	
	my ($max_no_of_peptides,$mhc2,$no_of_allele_used,$no_of_allele_not_used) = &get_max_no_of_peptides;
	
	if (($abs_min_score < $threshold) && ($total_no_of_peptides > 0)){
		print OUTPUT_FILE $ID . "\t" . $abs_min_score . "\t " . $h_avg_min_score {$ID.$abs_min_score} . "\t" . $method_used  . "\t" . $peptide_used . "\t" . $mhc2 . "\t" . $max_no_of_peptides . 		"\t" . $total_no_of_peptides . "\t" . $no_of_mhc_alleles  . "\t" . $no_of_allele_used . "\t" . $no_of_allele_not_used . "\t" . " PEPTIDE\n";
	} else {
		print OUTPUT_FILE $ID . "\t " . $abs_min_score . "\t" . $h_avg_min_score {$ID.$abs_min_score} . "\t" .$method_used  . "\t" . $peptide_used . "\t" . $mhc2 . "\t" . $max_no_of_peptides . 			"\t" . $total_no_of_peptides  . "\t" . $no_of_mhc_alleles . "\t" . $no_of_allele_used . "\t" . $no_of_allele_not_used . "\n";
	}
}

sub get_abs_min_score ()
{
	#get the scores for the min percentile
	my ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score,$Sturniolo_score) = split (',', $min_percentile_scores);
	
	print DEBUG_FILE "Check_in_abs " . $ID . "\t" .  $ANN_score . "\t" . $SMM_score . "\t" .	$CombLib_score . "\t" .	$NetMHCpan_score. "\t" . $Sturniolo_score . "\t" . $binding_info  . "\n";
	
	#get the peptides for the min percentile
	my ($ANN_peptide,$SMM_peptide,$CombLib_peptide,$NetMHCpan_peptide,$Sturniolo_peptide) = split (',', $min_percentile_peptides);
	
	#Ignore CombLib score
	if ($ignore_CombLib == 1) {
		$CombLib_score = "-";
	}
		
	$h_ANN_min {$binding_info} = $ANN_score;
	$h_SMM_min {$binding_info} = $SMM_score;
	$h_CombLib_min {$binding_info} = $CombLib_score;
	$h_NetMHCpan_min {$binding_info} = $NetMHCpan_score;
	$h_Sturniolo_min {$binding_info} = $Sturniolo_score;
	
	#Record the peptide for each score
	$h_ANN_peptide {$ANN_peptide . "?" . $binding_info} = $ANN_score;
	$h_SMM_peptide {$SMM_peptide . "?" . $binding_info} = $SMM_score;
	$h_CombLib_peptide {$CombLib_peptide . "?" . $binding_info} = $CombLib_score;
	$h_NetMHCpan_peptide  {$NetMHCpan_peptide . "?" . $binding_info} = $NetMHCpan_score;
	$h_Sturniolo_peptide  {$Sturniolo_peptide. "?" . $binding_info} = $Sturniolo_score;
			
	my @min_score = ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score,$Sturniolo_score);
	
	my $no_of_scores = 0;
	my $abs_min_changed = 0;	
	foreach my $score (@min_score) {
		#check for a dash
		if ($score ne "-" ) {
		
			$avg_min_score = $avg_min_score + $score;
			$no_of_scores++;
			
			if ($score < $abs_min_score) {
				$abs_min_score = $score;
				$abs_min_changed = 1;
			}
		}
	}
	
	if ( $abs_min_changed == 1) {
		#print DEBUG_FILE $avg_min_score . "\t" . $no_of_scores . "\n";
		$avg_min_score =  $avg_min_score/$no_of_scores;
		$h_avg_min_score {$ID.$abs_min_score} = sprintf ("%.2f",$avg_min_score);
		}
	print DEBUG_FILE "##MIN" . $ID . "\t" . $abs_min_score . "\t" . $h_avg_min_score {$ID.$abs_min_score} . "\n";
	$avg_min_score = 0;
}

sub get_binding_info ()
{
	#ANN
	foreach my $binder (sort keys  %h_ANN_min) {
	
		#print DEBUG_FILE $binder  . "\t" . $h_ANN_min {$binder} . "\t" . $abs_min_score . "\tANN\n";
		if ($abs_min_score  == $h_ANN_min {$binder} ) {
			return  $binder . "\tANN" ;
		}
	}
	
	#SMM
	foreach my $binder  (sort keys  %h_SMM_min) {
		#print DEBUG_FILE $binder . "\t" . $h_SMM_min {$binder} . "\t" . $abs_min_score . "\tSMM\n";
		if ($abs_min_score  == $h_SMM_min {$binder} ) {
			return  $binder . "\tSMM" ;
		}
	}
	
	#CombLib
	foreach my $binder  (sort keys  %h_CombLib_min) {
		#print DEBUG_FILE $binder . "\t" . $h_SCombLib_min {$binder} . "\t" . $abs_min_score . "\tCombLib\n";
		if ($abs_min_score  == $h_CombLib_min {$binder} ) {
			return  $binder . "\tCombLib" ;
		}
	}
	
	#NetMHCpan
	foreach my $binder  (sort keys  %h_NetMHCpan_min) {
		#print DEBUG_FILE $binder . "\t" . $h_NetMHCpan_min {$binder} . "\t" . $abs_min_score . "\tNetMHCpan\n";
		if ($abs_min_score  == $h_NetMHCpan_min {$binder} ) {
			return  $binder . "\tNetMHCpan" ;
		}
	}
	
	#Sturniolo
	foreach my $binder  (sort keys  %h_Sturniolo_min) {
		#print DEBUG_FILE $binder . "\t" . $h_NetMHCpan_min {$binder} . "\t" . $abs_min_score . "\tNetMHCpan\n";
		if ($abs_min_score  == $h_Sturniolo_min {$binder} ) {
			return  $binder . "\tSturniolo" ;
		}
	}
}


sub get_peptide ()
{
	
	#ANN
	foreach my $peptide_mhc (sort keys  %h_ANN_peptide) {
	
		#print DEBUG_FILE $peptide_mhc  . "\t" . $h_ANN_peptide {$peptide_mhc} . "\t" . $abs_min_score . "\tANN\n";
		if ($abs_min_score  == $h_ANN_peptide {$peptide_mhc} ) {
		my ($peptide, $mhc) = split ('\?',$peptide_mhc);
			return  $peptide;
		}
	}
	
	#SMM
	foreach my $peptide_mhc (sort keys  %h_SMM_peptide) {
	
		#print DEBUG_FILE $peptide_mhc   . "\t" . $h_SMM_peptide {$peptide_mhc } . "\t" . $abs_min_score . "\tSMM\n";
		if ($abs_min_score  == $h_SMM_peptide {$peptide_mhc} ) {
			my ($peptide, $mhc) = split ('\?',$peptide_mhc);
			return  $peptide;
		}
	}
	
	#CombLib
	foreach my $peptide_mhc (sort keys  %h_CombLib_peptide) {
	
		#print DEBUG_FILE $peptide_mhc   . "\t" . $h_CombLib_peptide {$peptide_mhc } . "\t" . $abs_min_score . "\tCombLib\n";
		if ($abs_min_score  == $h_CombLib_peptide {$peptide_mhc} ) {
			my ($peptide, $mhc) = split ('\?',$peptide_mhc);
			#print DEBUG_FILE $peptide_mhc  . "\t" . $h_CombLib_peptide {$peptide_mhc} . "\t" . $abs_min_score . "\t" . $peptide ."ANN\n";
			return  $peptide;
		}
	}
	
	#NetMHCpan
	foreach my $peptide_mhc (sort keys  %h_NetMHCpan_peptide) {
	
		#print DEBUG_FILE $peptide_mhc   . "\t" . $h_NetMHCpan_peptide {$peptide_mhc } . "\t" . $abs_min_score . "\tCombLib\n";
		if ($abs_min_score  == $h_NetMHCpan_peptide {$peptide_mhc} ) {
			my ($peptide, $mhc) = split ('\?',$peptide_mhc);
			return  $peptide;
		}
	}
	
	#Sturniolo
	foreach my $peptide_mhc (sort keys  %h_Sturniolo_peptide) {
	
		#print DEBUG_FILE $peptide_mhc   . "\t" . $h_NetMHCpan_peptide {$peptide_mhc } . "\t" . $abs_min_score . "\tCombLib\n";
		if ($abs_min_score  == $h_Sturniolo_peptide {$peptide_mhc} ) {
			my ($peptide, $mhc) = split ('\?',$peptide_mhc);
			return  $peptide;
		}
	}

}

sub get_max_no_of_peptides ()
{

	my $max_no_of_peptides = 0;
	my $return_mhc = 0;
	
	my $no_of_allele_used = 0;
	my $no_of_allele_not_used = 0;
	
	foreach my $mhc (sort keys  %h_no_of_peptides) {
	
		#check if MHC allele has binding site
		if ($h_no_of_peptides {$mhc} > 0) {
			$no_of_allele_used++;
		} else {
			$no_of_allele_not_used++;
		}
		#print DEBUG_FILE  $mhc_length . "\t" . $h_no_of_peptides {$mhc_length} . "\n";
	
		if ($h_no_of_peptides {$mhc} > $max_no_of_peptides ) {
			$max_no_of_peptides = $h_no_of_peptides {$mhc};
			$return_mhc = $mhc;
		}
	}

	return ($max_no_of_peptides,$return_mhc,$no_of_allele_used,$no_of_allele_not_used);
	
	
}
