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
my $id_allele_peptide_file = $output_dir . "id_allele_peptide.txt";
my $id_allele_file = $output_dir . "id_allele.txt";
my $allele_peptide_file = $output_dir . "allele_peptide.txt";


## hard coding
my $program = "mhc1";
my $threshold = 50;
my $ignore_CombLib = 1;


#Global variables
my $ID;
my $binding_info;
my $start = 1;
my $min_Percentile = 999999;
my $ANN_score = 9999999;
my $SMM_score = 9999999;
my $CombLib_score = 9999999;
my $NetMHCpan_score =9999999;
my $abs_min_score = 999999999;
my $min_percentile_scores;
my $total_no_of_peptides = 0;
my $no_of_mhc_alleles = 0;
my $no_of_proteins = 0;
my $allele_used_with_min;
my $peptide_used_with_min;
my $method_used_with_min;
my $no_binding_data_flag = 0;

#create hashes
my %h_avg_min_score;
my %h_no_of_peptides;
my %h_ID_no_of_peptides;

my %h_global_no_of_peptides;
my %h_not_used;

my %h_prom_peptides;
my %h_ID_peptide;
my %h_ID_allele_peptide;
my %h_allele_peptide;

#Open debug file
#pen DEBUG_FILE,'>debug.txt';


#File handling
open INPUT_FILE,$input_file ;
open OUTPUT_FILE,'>' . $output_file;
open STATS_FILE,'>' . $stats_file;
open PEPTIDE_FILE,'>' . $peptide_file;
open ID_PEPTIDE_FILE,'>' . $id_peptide_file;
open ID_ALLELE_PEPTIDE_FILE,'>' . $id_allele_peptide_file;
open ID_ALLELE_FILE,'>' . $id_allele_file;
open ALLELE_PEPTIDE_FILE,'>' . $allele_peptide_file;

print OUTPUT_FILE "#1=UniProtID :2=Min score :3=score average:4=MHC1 allele :5=Peptide length :6=Method used :7=Peptide sequence :8&9=Allele-length combination with most binding sites\n";
print OUTPUT_FILE "#10=Max. number of binding sites  :11=Total number of binding sites for all alleles :12=Total number of MHC1 allele-length combinations available\n";
print OUTPUT_FILE "#13=No. of MHC1 allele-length combinations used :14=No. of MHC1 allele-length combinations not used\n"; 	

#LOOP through protein evidence set (the map file)
while (defined (my $line = <INPUT_FILE>)) {


#>tr|Q27298|Q27298_TOXGO P30 OS=Toxoplasma gondii GN=p30 PE=4 SV=1

	$line = &trim ($line);

	$line =~ s/\x0D//g; #strips ^M characters

	#Check for # or Accession
	if ($line =~ m/^#/) { 
	
			next;
			
	} elsif ($line =~ m/allele/) { 
	
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
		$no_of_proteins++;
		#Reset
		$abs_min_score = 999999999;
		$total_no_of_peptides = 0;
		$no_of_mhc_alleles = 0;
		$allele_used_with_min = 0; 
		$peptide_used_with_min = 0;
		$no_binding_data_flag = 0;
		
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
		
		my ($MHC,$PeptideLength) = split ('\s+',$binding_info);
		$h_no_of_peptides {$MHC . "\t" . $PeptideLength} = 0;
		
		#Assign zero if empty
		if (&check_for_empty_string ($h_global_no_of_peptides {$MHC . "\t" . $PeptideLength})) {
			$h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} = 0;
		}
		
		#Assign zero if empty
		if (&check_for_empty_string ($h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength})) {
			$h_ID_no_of_peptides {$ID}{$MHC . "\t" . $PeptideLength} = 0;
		}
				
		#reset the method scores
		$min_Percentile = 9999999;
		$no_of_mhc_alleles++;
		
	# Ignore line if no matching species				
	} elsif ($line =~ m/Could not find tools matching species/) {

		$h_not_used {$binding_info} = $h_not_used {$binding_info} + 1;
		
		next;
			
	} else {
	
		$no_binding_data_flag = 1;

		#allele	length	peptide	method	percentile_rank	ann_ic50	ann_rank	smm_ic50	smm_rank	comblib_sidney2008_score	comblib_sidney2008_rank	netmhcpan_ic50	netmhcpan_rank

		my ($MHC,$PeptideLength,$Peptide,$Percentile_Rank,$Method_used,
		$ANN_IC50,$ANN_Rank,$SMM_IC50,$SMM_Rank,$CombLib_IC50,$CombLib_Rank,$NetMHCpan_IC50,$NetMHCpan_Rank) = split ('\t',$line);


				
		if ($Method_used =~ m/Consensus\s+\((.+)\)/) {
		
			my @methods = split (',',$1);
			
			my $no_of_methods = scalar @methods;
			my $no_below_threshold = 0;
					
			my $score_test = $threshold;
			$Method_used = "NONE";
			
			foreach (@methods) {
			
				if (($_ eq "ANN") &&  ($ANN_IC50 < $threshold)) {
					$no_below_threshold++;
				}elsif (($_ eq "SMM") && ($SMM_IC50 < $threshold)) {
					$no_below_threshold++;
				# check CombLib_Sidney2008
				}elsif (($_ =~ m/CombLib/) && ($CombLib_IC50 < $threshold)) {
					## THIS IS not a true check and needs to be changed
					$no_below_threshold++;
				}elsif (($_ eq "NetMHCpan") && ($NetMHCpan_IC50 < $threshold)) {
					$no_below_threshold++;
				}
				
				#find lowest score
				if (($_ eq "ANN") && ($ANN_IC50 < $score_test)) {
					$score_test = $ANN_IC50;
					$Method_used = $_;
									
				}elsif (($_ eq "SMM") && ($SMM_IC50 < $score_test)) {
					$score_test = $SMM_IC50;
					$Method_used = $_;
									
				}elsif (($_ eq "NetMHCpan") && ($NetMHCpan_IC50 < $score_test)) {
					$score_test = $NetMHCpan_IC50;
					$Method_used = $_;
				} 
				
			}
						
			#Check for consensus peptide
			if ($no_below_threshold ==  $no_of_methods) {
				
				$h_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
				$h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
				$h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} = $h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} + 1;
				$total_no_of_peptides++;
				$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
				$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
				$h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide} = $h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide}+ 1;
				
				if (&check_for_empty_string ($h_allele_peptide  {$MHC . "\t" . $PeptideLength} )) {
					$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $Peptide ;
				} else {
					$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $h_allele_peptide  {$MHC . "\t" . $PeptideLength} . "\t" . $Peptide ;
				}
										
			}
		
		} elsif (($Method_used eq "NetMHCpan") && ($NetMHCpan_IC50 ne "-")) {
				
				if ($NetMHCpan_IC50 < $threshold) {
					$h_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} = $h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide} = $h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide}+ 1;
					
					if (&check_for_empty_string ($h_allele_peptide  {$MHC . "\t" . $PeptideLength} )) {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $Peptide ;
					} else {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $h_allele_peptide  {$MHC . "\t" . $PeptideLength} . "\t" . $Peptide ;
					}
				}
			
		} elsif (($Method_used eq "ANN") && ($ANN_IC50 ne "-")) {
				
				if ($ANN_IC50 < $threshold) {
					$h_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} = $h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide} = $h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide}+ 1;
					
					if (&check_for_empty_string ($h_allele_peptide  {$MHC . "\t" . $PeptideLength} )) {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $Peptide ;
					} else {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $h_allele_peptide  {$MHC . "\t" . $PeptideLength} . "\t" . $Peptide ;
					}
				}
			
		} elsif (($Method_used eq "SMM") && ($SMM_IC50 ne "-")) {
				
				if ($SMM_IC50 < $threshold) {
					$h_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} = $h_global_no_of_peptides {$MHC . "\t" . $PeptideLength} + 1;
					$h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} = $h_ID_no_of_peptides {$ID} {$MHC . "\t" . $PeptideLength} + 1;
					$total_no_of_peptides++;
					$h_prom_peptides {$Peptide} = $h_prom_peptides {$Peptide} + 1;
					$h_ID_peptide {$ID} {$Peptide} = $h_ID_peptide {$ID} {$Peptide} + 1;
					$h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide} = $h_ID_allele_peptide {$ID} {$MHC . "\t" . $PeptideLength} {$Peptide}+ 1;
					
					if (&check_for_empty_string ($h_allele_peptide  {$MHC . "\t" . $PeptideLength} )) {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $Peptide ;
					} else {
						$h_allele_peptide  {$MHC . "\t" . $PeptideLength} = $h_allele_peptide  {$MHC . "\t" . $PeptideLength} . "\t" . $Peptide ;
					}
					
				}
			
		}
		
		
		#Get minimum percentile
		if ($Percentile_Rank < $min_Percentile) {
			$min_Percentile = $Percentile_Rank;
			$min_percentile_scores = $ANN_IC50 . "," . $SMM_IC50 . "," . $CombLib_IC50 . "," . $NetMHCpan_IC50 . "," . $Peptide . "," . $Method_used;
		}
	
	}

			
} #end of LINE:while

#print last prediction
&print_out;

my $allele_count = 0;
my $not_used = 0;
my $used = 0;
my $count = 0;

print STATS_FILE "### Frequency of MHC1 alleles used ###\n";
foreach my $key ( sort  { $h_global_no_of_peptides {$b} <=> $h_global_no_of_peptides {$a}} keys  %h_global_no_of_peptides) {

	if ($h_global_no_of_peptides {$key} > 0) {
		$used++;
	} else {
		$not_used++;
	}

	$count = $count + $h_global_no_of_peptides {$key}; 
	print STATS_FILE  $key . "\t" . $h_global_no_of_peptides {$key} . "\n";
	
	$allele_count++;
	
}

print STATS_FILE "\nNumber of MHC1 alleles used = $used\n";
print STATS_FILE "Number of MHC1 alleles not used = $not_used\n";
print STATS_FILE "Number of MHC1 alleles available = $allele_count\n";
print STATS_FILE "Number of epitopes = $count\n";
print STATS_FILE "Number of proteins = $no_of_proteins\n";

$allele_count = 0;

print STATS_FILE "\n### No methods found for following MHC1 alleles ###\n";
#Print alleles not used
foreach my $key (sort { $h_not_used {$b} <=> $h_not_used {$a}} keys   %h_not_used) {

 print STATS_FILE $key . "\t" . $h_not_used {$key} . "\n";
 
 $allele_count++;
	
}

#Print alleles not used
print STATS_FILE "\nNumber of MHC1 alleles not used = $allele_count\n";


print ID_PEPTIDE_FILE "### Proteins with common peptides###\n";

foreach my $id ( sort keys %h_ID_peptide ) {

	print ID_PEPTIDE_FILE "ID:" . $id . "\n";
	
	my $child_hashref = $h_ID_peptide{$id};
   
		foreach my $peptide (sort { $child_hashref-> {$b} <=> $child_hashref->{$a}}   keys   %$child_hashref) {
	
		print ID_PEPTIDE_FILE $id . "\t" . $peptide . "\t" . $child_hashref-> {$peptide} . "\n";
	}
}

$count = 0;
$allele_count = 0;

print PEPTIDE_FILE "### Common sub-sequences used as epitopes###\n";

foreach my $key (sort { $h_prom_peptides{$b} <=> $h_prom_peptides{$a}} keys   %h_prom_peptides) {

$count = $count + $h_prom_peptides {$key};

 print PEPTIDE_FILE $key . "\t" . $h_prom_peptides {$key} . "\n";
  $allele_count++;
	
}

print PEPTIDE_FILE "\n#Number of Common peptides = $allele_count\n";
print PEPTIDE_FILE "\n#Number of epitopes = $count\n";

foreach my $id ( sort keys %h_ID_allele_peptide ) {

	print ID_ALLELE_PEPTIDE_FILE "ID:" . $id . "\n";
	
	my $child_hashref = $h_ID_allele_peptide {$id};
   
	foreach my $allele (sort { $child_hashref-> {$b} <=> $child_hashref->{$a}}   keys   %$child_hashref) {
	
		my $child2_hashref = $h_ID_allele_peptide {$id} {$allele};
	
		print ID_ALLELE_PEPTIDE_FILE $id . "\t" . $allele  . "\n";
			
		foreach my $peptide (sort { $child2_hashref-> {$b} <=> $child2_hashref->{$a}}   keys   %$child2_hashref) {
	
			print ID_ALLELE_PEPTIDE_FILE "\t" . $peptide . "\t" . $h_ID_allele_peptide {$id} {$allele} {$peptide}. "\n";
			
		}
	}
}

my $ID_count = 0;

foreach my $id ( sort keys %h_ID_no_of_peptides ) {

	$allele_count = 0;
	$not_used = 0;
	$used = 0;
	$count = 0;

	print ID_ALLELE_FILE "ID:" . $id . "\n";
	
	my $child_hashref = $h_ID_no_of_peptides {$id};
   
	foreach my $allele (sort { $child_hashref-> {$b} <=> $child_hashref->{$a}}   keys   %$child_hashref) {
	
			if ($child_hashref-> {$allele} > 0) {
				$used++;
			} else {
				$not_used++;
			}

			$count = $count + $child_hashref-> {$allele}; 
	
		print ID_ALLELE_FILE $id . "\t" . $allele . "\t" . $child_hashref-> {$allele} . "\n";
		$allele_count++;
	}

	print ID_ALLELE_FILE "\nNumber of MHC1 alleles used = $used\n";
	print ID_ALLELE_FILE "Number of MHC1 alleles not used = $not_used\n";
	print ID_ALLELE_FILE "Number of MHC1 alleles available = $allele_count\n";
	print ID_ALLELE_FILE "Number of epitopes = $count\n\n";
	$ID_count++;
	
}

print ID_ALLELE_FILE "\nNumber of proteins = $ID_count\n";

foreach my $mhc_length ( sort keys %h_allele_peptide ) {

	print ALLELE_PEPTIDE_FILE $mhc_length . "\n";
	
	my @peptide_array = split ('\s+',$h_allele_peptide {$mhc_length});
	
	foreach my $peptide (@peptide_array) {
	
		print ALLELE_PEPTIDE_FILE "\t" . $peptide . "\n";
	
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
	
	if ($no_binding_data_flag == 0) {
	
	print OUTPUT_FILE $ID . "\t9999\t9999\t-\t-\t-\tNO_method_used\t-\t-\t-\t0\t0\t" . $no_of_mhc_alleles  . "\t0\t" . $no_of_mhc_alleles . "\tNO_DATA\n";
		
	} else {
	
		#get the abslotue minimum score
		&get_abs_min_score;
		
		my ($max_no_of_peptides,$mhc_length,$no_of_allele_used,$no_of_allele_not_used) = &get_max_no_of_peptides;
			
		if (($abs_min_score < $threshold) && ($total_no_of_peptides > 0)){
				
			print OUTPUT_FILE $ID . "\t" . $abs_min_score . "\t " . $h_avg_min_score {$ID.$abs_min_score} . "\t" . $allele_used_with_min  . "\t";
			print OUTPUT_FILE $method_used_with_min . "\t" . $peptide_used_with_min . "\t" . $mhc_length . "\t" . $max_no_of_peptides . "\t"; 
			print OUTPUT_FILE $total_no_of_peptides . "\t" . $no_of_mhc_alleles  . "\t" . $no_of_allele_used . "\t" . $no_of_allele_not_used . "\t" . " PEPTIDE\n";
		} else {

			print OUTPUT_FILE $ID . "\t" . $abs_min_score . "\t " . $h_avg_min_score {$ID.$abs_min_score} . "\t" . $allele_used_with_min . "\t" . $method_used_with_min . "\t";
			print OUTPUT_FILE $peptide_used_with_min . "\t" . $mhc_length . "\t" . $max_no_of_peptides . "\t" . $total_no_of_peptides . "\t" . $no_of_mhc_alleles  . "\t";
			print OUTPUT_FILE $no_of_allele_used . "\t" . $no_of_allele_not_used . "\n";
		}
	
	}
	
	
}

sub get_abs_min_score ()
{
	#get the scores for the min percentile
	my ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score,$Peptide,$Method_used) = split (',', $min_percentile_scores);
	
		
	#Ignore CombLib score
	if ($ignore_CombLib == 1) {
		$CombLib_score = "-";
	}
	
			
	my @min_score = ($ANN_score,$SMM_score,$CombLib_score,$NetMHCpan_score);
	
	my $no_of_scores = 0;
	my $abs_min_changed = 0;
	my $avg_min_score = 0;
	
	foreach my $score (@min_score) {
		
		#check for a dash and empty score
		if (($score ne "-" ) && (!&check_for_empty_string ($score))) {
		
			$avg_min_score = $avg_min_score + $score;
			$no_of_scores++;
			
			if ($score < $abs_min_score) {
				$abs_min_score = $score;
				$allele_used_with_min = $binding_info;
				$peptide_used_with_min = $Peptide;
				$method_used_with_min = $Method_used;
				$abs_min_changed = 1;
			}
		}
	}
	
	if ( $abs_min_changed == 1) {
		$avg_min_score =  $avg_min_score/$no_of_scores;
		$h_avg_min_score {$ID.$abs_min_score} = sprintf ("%.2f",$avg_min_score);
	}
}


sub get_max_no_of_peptides ()
{

	my $max_no_of_peptides = 0;
	my $return_mhc_length = 0;
	
	my $no_of_allele_used = 0;
	my $no_of_allele_not_used = 0;
	
	foreach my $mhc_length (sort keys  %h_no_of_peptides) {
	
		#check if MHC allele has binding site
		if ($h_no_of_peptides {$mhc_length} > 0) {
			$no_of_allele_used++;
		} else {
			$no_of_allele_not_used++;
		}
			
		if ($h_no_of_peptides {$mhc_length} > $max_no_of_peptides ) {
			$max_no_of_peptides = $h_no_of_peptides {$mhc_length};
			$return_mhc_length = $mhc_length;
		}
	}

	if ($max_no_of_peptides == 0) {
	
		$return_mhc_length = "-\t-";
	
	}
	return ($max_no_of_peptides,$return_mhc_length,$no_of_allele_used,$no_of_allele_not_used);
	
	
}
