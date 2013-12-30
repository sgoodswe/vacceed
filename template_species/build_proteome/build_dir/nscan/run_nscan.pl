#!/usr/bin/perl
use strict;
use File::HomeDir;
#use warnings;
### PROGRAM run_nscan.pl

# Get the comand arguments
my $prefix = $ARGV[0];
my $input_file = $ARGV[1];
my $informant_file = $ARGV[2];
my $output_dir = $ARGV[3] ;
my $chromosome = $ARGV[4];
my $log_file = $ARGV[5];

## hard coding
#$prefix = "chrIa";
#$input_file = "nscandriver.config";
#$chromosome = "/home/sgoodswe/vaccine_finder/neospora/build_proteome/chromosomes/chrIa.fasta";
#$informant_file = "/home/sgoodswe/vaccine_finder/toxoplasma/build_folder/chromosomes/chrIa.fasta";
#$output_dir = "/home/sgoodswe/vaccine_finder/neospora/build_proteome/build_folder/nscan/output";
my $informant_file_name =  "informant.fa";


#Global variables
my @config_lines;
my $tmp_dir;
my $output_file; 
my $output_informant_file;
my $system_command;

#create hashes


#Open debug file
#open DEBUG_FILE,'>debug.txt';
open LOG_FILE,'>>' . $log_file;

#get home directory
my $home_dir = File::HomeDir->my_home;


#delete any existing align file in output dir
$system_command = "rm -f $output_dir/$prefix*align >> $log_file 2>&1";
system("$system_command");


#File handling
open INPUT_FILE,$input_file ;


## Read the configuration file
#LOOP 
while (defined (my $line = <INPUT_FILE>)) {

			
	$line = &trim ($line);

	#get the tmp folder from the config file
	if ($line =~ m/tmpdir\=(.+)/) { 

		$tmp_dir = $1;

		#Check if temp dir contains $HOME
		if ($tmp_dir =~ m/\$HOME/){
			$tmp_dir =~ s/\$HOME/$home_dir/;
		}
		

		$tmp_dir = $tmp_dir . "tmp_$prefix";
		push (@config_lines, "tmpdir=" . $tmp_dir . "/");
	
	
	#Change the informant line
	} elsif ($line =~ m/^informant\=/) { 
	
		push (@config_lines, "informant=informant.fa");

	} else {
	
		push (@config_lines, $line);
	}
	
			
} #end of LINE:while


#Delete and remove tmp dir
$system_command = "mkdir -p $tmp_dir >> $log_file 2>&1";
system("$system_command");


# Change to new tmp dir
chdir  ($tmp_dir);

my $config_file = "nscan_" . $prefix . ".config";
$output_file = $tmp_dir . "/" . $config_file;

open OUTPUT_FILE,'>' . $output_file;

#print out the temporary config file to tmp directory
foreach my $line (@config_lines) {

	print OUTPUT_FILE $line . "\n";

}


#Read the informant file and change the header
open INFORMANT_FILE,$informant_file;

open OUTPUT_INFORMANT_FILE,'>' . $tmp_dir . "/" . $informant_file_name;


#LOOP 
while (defined (my $line = <INFORMANT_FILE>)) {

	#Check for #
	if ($line =~ m/^#/) { next;}
		
	$line = &trim ($line);
	
	#Check for the informant line
	if ($line =~ m/>/) { 
	
		print OUTPUT_INFORMANT_FILE ">informant\n";

	} else {
	
		print OUTPUT_INFORMANT_FILE $line . "\n";
	}
	
			
} #end of LINE:while


#### Run NSCAN
$system_command = "Nscan_driver.pl  --nomask -d $output_dir $chromosome $config_file";

print  LOG_FILE $system_command . "\n";

my $rc = system("$system_command >> $log_file 2>&1");

print  LOG_FILE "Return from running NSCAN " . $rc . "\n";

 

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
