### Function calls


#### makee Prediction function  ####
nb_makePred <- function ()
{

	#training data
	train <- dataset
	
	#test data
	prediction <- 1:nrow (candidates)
	
	#Target variable
	target  <- "Candidate"
	
	nb_model <- naiveBayes (train[,input], train[,target] )
	

	# Make a prediction using the naiveBayes model
	pr_class <- predict(nb_model,(candidates[prediction, c(input, target)]),type="raw")
	
			
	# Extract the relevant variables from the dataset.
	sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))
	

	# Output the combined data.
	output <- cbind(sdata,pr_class)
	

}


