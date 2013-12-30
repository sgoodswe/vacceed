#!/usr/bin/perl
use strict;
use Cwd;
use warnings;
use Config::Simple;
use File::HomeDir;


#####################################
# Version 1-00-00
#####################################
#Open debug file
#open DEBUG_FILE,'>debug.txt';

# To run program enter: 
#	cleanup.pl <species code>
#	 e.g species code = nc for neospora

#general global variables
my $program = "cleanup.pl";
my $run_type;
my $code_entered;
my $hash_search_variables;
my @delete;
my $resource_list="";

# variables for startup.ini
my $code;
my $species;
my $type;
my $config_file;
my $config_dir;

#variables for config file
my $cfg;
my $key;
my $hash;



#Get the number of  command-line arguments
my $numArgs = $#ARGV;

if ($numArgs == -1) {

	print "\nERROR: No command-line arguments entered\n";
	print "To clean up vaccine discovery pipeline directories enter species code\n";
	print "\tExample: perl cleanup.pl nc\n\n";
	print "To clean up build proteome directories enter build + species code \n";
	print "\tExample: perl cleanup.pl build nc\n";
	exit

} else {

   $code_entered = lc ($ARGV[0]);
   
   #Check if build + species code is entered
   if (($code_entered eq "build") && ($numArgs == 0)) { 
		$run_type = "BUILD";
		
		print "\nERROR: No species code entered in command-line\n";
		print "To clean up build proteome directories enter build + species code \n";
		print "\tExample: perl cleanup.pl build nc\n";
		exit
	} elsif ($code_entered eq "build") { 
		$run_type = "BUILD";
		
		#get the species code
		$code_entered = lc ($ARGV[1]);

	} else {
		$run_type = "PIPELINE";
	}
}


# Read the startup.ini
open INPUT_FILE,"startup.ini";

#   startup.ini is in the following format
#code<species<type e.g. build or pipeline<config<dir
#nc<Neospora caninum<build<neospora.ini<neospora

while (defined (my $line = <INPUT_FILE>)) 
{
		
	chomp ($line);
	($code,$species,$type,$config_file,$config_dir) = split (/</,$line);

	#conver code to lower case species code		
	$code  = lc ($code);
	
	#Only do the following  if the species code entered matches the species in the ini file
	if ($code_entered eq $code) {
	
		if (($run_type eq "BUILD") && (lc ($type) ne "build")) {
			$config_file="";
		} 
		last;
		
	} else {
		$code = "";
	}		
}

# Check if we have the required data to start
if (!$code)
{
	print "\nERROR: Species code not found in startup.ini\n";
	exit
}

#check for configuration file
if (!$config_file) {
	print "\nERROR: No $run_type configuration file in startup.ini\n";
	exit
}

#check if config file exists
my $config = $config_dir . "/" . $config_file;
if (-e $config) 
{
	#### Open and  read the configuration file
	$cfg = new Config::Simple();
	$cfg->read($config);
}
else
{
	print "\nERROR: configuration file $config_file cannot be found in directory $config_dir\n";
	exit
}

#get home directory
my $home_dir = File::HomeDir->my_home;

######################
# READ config file 
######################
$key = 'Resources.name';
my @resources = $cfg->param($key);

### Get the Main block
$key =  'Main';
my $hash_main_variables = $cfg->param(-block=>$key);

### Get the Variable block
$key =  'Variables';
my $hash_user_variables = $cfg->param(-block=>$key);

#get the proteome directory
$key =  'Variables.proteome_dir';
my $proteome_dir = $cfg->param($key);

##### END OF READING CONFIG GENERAL SECTION #####

#get home directory
my $home_dir = File::HomeDir->my_home;

# To run program enter: 
print "Starting Program $program \n";


print "\n###########################################################################\n";
print "# cleanup.pl deletes all files from the following resource sub-directories:\n";
print "#\toutput\n";
print "#\tsummary_files\n";
print "# All files from the proteome directory are also deleted\n";
print "#\n";
print "# You can either delete ALL files at once OR\n"; 
print "# be prompted to delete each directory in turn\n"; 
print "#\n";
print "############################################################################\n\n";
print "Enter A to delete all files OR P for prompts (A/P) [A]:" ; 
my $requested_action = <STDIN>;
chomp $requested_action;

#check for default
if (&check_for_empty_string ($requested_action)) {
	$requested_action = "A";
}

##########################################
##########################################
#   START LOOP for Resources
##########################################
##########################################
# Loop for each resource in config file
foreach (@resources) 
{
	
	#Get the resource name
	my $resource_name = $_ ;

	$resource_list = $resource_list . "\t" . $resource_name;
		
	#Get all from [resource]  
	$key = $resource_name;
	my $hash_resource_variables = $cfg->param(-block=>$key);
	
	#Add the correct resource name  to the user variable
	$hash_user_variables->{"resource_dir"} = lc ($resource_name);
	
	#Initialise hash search  variables
	$hash_search_variables = {};
	$hash_search_variables = {%$hash_main_variables, %$hash_user_variables,%$hash_resource_variables};
		
	#Get the  output directory from config file
	$key = $resource_name . '.out_dir';
	my $out_dir = $cfg->param($key);

	if (!&check_for_empty_string (defined $out_dir)) {
		
		#Check for variables
		my $dir = &check_for_variables ($out_dir);

		
	
		## User requested to delete ALL
		if (uc $requested_action eq "A") {

			print "direcotru " . $dir . "\n";
	
			push (@delete,$dir);
	
		} else {
	
			&delete_directory($dir);
		}
	}

	#Get the  summary_files directory from config file
	$key = $resource_name . '.sum_dir';
	my $sum_dir = $cfg->param($key);

	if (!&check_for_empty_string (defined $sum_dir)) {
		
		#Check for variables
		my $dir = &check_for_variables ($sum_dir);
	
		## User requested to delete ALL
		if (uc $requested_action eq "A") {
	
			push (@delete,$dir);
	
		} else {
	
			&delete_directory($dir);
		}
	}		
	
} # Loop for each resource  - foreach (@resources)

if (uc $requested_action eq "A") {

	print "\n" . 'You are about to delete contents of "output" and "summary_files" for resources:' ; 
	print "\n$resource_list\n\n" ;
	print "Are you sure (Y/N) [N]? " ;
	my $action = <STDIN>;
	chomp $action;
	
	#check for default
	if (uc $action eq "Y") {

		foreach my $dir (@delete) {

			#Check if directory contains $HOME
			if ($dir =~ m/\$HOME/){
				$dir =~ s/\$HOME/$home_dir/;
			}

			#check if directory exists and is not empty
			if ((-d $dir) && (&check_for_empty_dir ($dir))) {

				#construct system command
				my $system_command = "rm $dir/*";
		
				#run system command
				my $ret = system ($system_command);
				&check_for_errors ($ret,$dir);
			}
		}
	}

}	


## Code to delete contents of proteome directory

#Check for variables
$proteome_dir = &check_for_variables ($proteome_dir);

#Check if directory contains $HOME
if ($proteome_dir =~ m/\$HOME/){
	$proteome_dir =~ s/\$HOME/$home_dir/;
}

print "You are about to delete contents of proteome directory\n"; 
print "Are you sure (Y/N) [N]? " ;
$requested_action = <STDIN>;
chomp $requested_action;

#check for default
if (uc $requested_action eq "Y") {

	#check if directory exists and is not empty
	if ((-d $proteome_dir) && (&check_for_empty_dir ($proteome_dir))) {
		my $ret = system ("rm $proteome_dir/*");

		&check_for_errors ($ret,$proteome_dir);
	}

}

print "\nFinished program $program Successfully \n";

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


# sub-routine to check for user variables
sub check_for_main_variables($)
{
	my $string = shift;
	
	#check for $ sign
	if ($string =~ m/\$/){
	
		foreach my $key (sort keys %{$hash_main_variables}) 
		{
					
			$string =~ s/\$$key/$hash_main_variables->{$key}/;
					
		}		
	}
	
	return $string;
}




# sub-routine to check for variables
sub check_for_variables($)
{
	my $string = shift;
	

	if ($string !~ m/\$/){
		return $string;
	}

	
	foreach my $change_from (sort keys %{$hash_search_variables}) 
	{
		
		if ($string =~ m/\$\b$change_from\b/){

			my $change_to = $hash_search_variables->{$change_from};
			
			
			#check for $ sign
			if ($change_to =~ m/\$/){
				$change_to = &check_for_variables($change_to);
			}

			$string =~ s/\$$change_from/$change_to/;
			
		}
	}

	#check for $ sign
	if (($string =~ m/\$/) && ($string !~ m/\$HOME/)) {
		#invalid variable
		print "\nERROR: an unknown variable exists in command line $string. Check $config_file\n";
		exit;
	} else {	
		return $string;
	}
}


# delete directory
sub delete_directory ($)
{
	my $dir = shift;

	
	print "Delete contents of $dir (Y/N)[Y]? ";
	my $response = <STDIN>; 
	chomp $response; 

	if ((&check_for_empty_string ($response)) || (uc $response eq "Y")) {

		#Check if directory contains $HOME
		if ($dir =~ m/\$HOME/){
			$dir =~ s/\$HOME/$home_dir/;
		}

		#check if directory exists and is not empty
		if ((-d $dir) && (&check_for_empty_dir ($dir))) {

			#construct system command
			my $system_command = "rm $dir/*";
	
			#run system command
			my $ret = system ($system_command);
			&check_for_errors ($ret,$dir);
		}
	
	}
}


# Check for system error
sub check_for_errors($)
{
	my ($error,$dir) = @_;
	
	if ($error != 0) {
			print "Warning: Could not delete files in directory $dir\n";
	}

	return 0;	
}


# Check for empy directory
sub check_for_empty_dir($)
{

	my $dir = shift;
	
	my $i=0;
	
	opendir(DIR, $dir); 
	my @files = readdir(DIR);

	foreach my $file (@files) {
		unless ($file =~ /^[.][.]?\z/) {
		$i++;
		}
	}

	closedir(DIR);
	     
	if ($i != 0) {
		return 1;
	} else {
		return 0;
	}
	
}
