cat("Running R Wrapper for Support Vector Machines\n")


#get arguments
args <- commandArgs(TRUE)
prog_dir <- args[1]
out_dir <- args[2]
train_file <- args[3]
candidate_file <- args[4]

#hard coding
no_of_runs <- 1

#Source files
file_name <- paste(sep="/", prog_dir, "svm_runPred.R")
source(file_name)

file_name <- paste(sep="/", prog_dir, "svm_makePred.R")
source(file_name)


# load the predictors
file_name <- paste(sep="/", out_dir, "ml_predictors.R")
source(file_name)


#Load the data
candidates <- read.csv(candidate_file, na.strings=c(".", "NA", "", "?"), strip.white=TRUE, encoding="UTF-8")
dataset <- read.csv(train_file, na.strings=c(".", "NA", "", "?"), strip.white=TRUE, encoding="UTF-8")


#Run the predictions
svm_runPred(no_of_runs,out_dir)
