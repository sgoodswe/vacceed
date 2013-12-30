### Function calls


#### makee Prediction function  ####
nn_makePred <- function ()
{
	#training data
	train <- 1:nrow (dataset)
	
	#test data
	prediction <- 1:nrow (candidates)

	#Target variable
	target  <- "Candidate"
	

	# Build the NNet model.
	nnet <- nnet(as.factor(Candidate) ~ .,
		data= dataset[train,c(input,target)],
		size=10, skip=TRUE, MaxNWts=10000, trace=FALSE, maxit=100)

	#Obtain the response from the Neural network model.
	pr_class <- predict(nnet, candidates[prediction, c(input, target)], type="class")
		
	# Extract the relevant variables from the dataset.
	sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))

	# Output the combined data.
	output <- cbind(sdata,pr_class)
		
}

