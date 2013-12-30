#!/usr/bin/perl
use strict;
#use warnings;

#get the command line arguments
my $prefix = $ARGV [0];
my $assembly_dir = $ARGV [1];

## hard coding
my $input_file = "";
my $gene_file =  $assembly_dir . "/" . $prefix . "_merged_gene.fasta";
my $prot_file = $assembly_dir . "/" . $prefix . "_merged_protein.fasta";

#Global variables


#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';


#File handling
open PROT_FILE, '>' . $prot_file;
open GENE_FILE, '>' . $gene_file;


#Get files from assmebly directory
opendir(DIR, $assembly_dir );
my @FILES= readdir(DIR); 

foreach my $file (@FILES) {

	if ($file =~ m/(.+).fasta/) { 
	
		#chrIa_aug_new_gene.fasta
		#split file name into components
		my ($chr,$resource,$new,$type) = split ('_',$1);
				
		if (($type eq "gene") || ($type eq "prot") || ($type eq "sim")  ){


			#check the chromosome #
			if ($chr eq $prefix) {
				&process_file ($file,$type);
			}
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

sub  process_file ($)
{
		my $file = $_[0];
       		my $type = $_[1];
		my $id;
			
		#open the file to read
		open INPUT_FILE, $assembly_dir . "/" . $file ;
		
		#Loop through file to read
		while (defined (my $line = <INPUT_FILE>)) {

			$line = trim($line);

			if ($line =~ m/>(.+)/) {
				$id = $1;
			} else {

			   	if (!&check_for_empty_string ($line)) {
					
					#process gene types or prot types
					if (($type eq "gene") || ($type eq "sim")){
						print GENE_FILE ">"  . $id . "\n";
						print GENE_FILE $line . "\n";
					} else {
						print PROT_FILE ">"  . $id . "\n";
						print PROT_FILE $line . "\n";
					}
				
				}
			}
		}
}
