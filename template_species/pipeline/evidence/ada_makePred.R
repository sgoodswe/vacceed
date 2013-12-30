### Function calls


#### makee Prediction function  ####
ada_makePred <- function ()
{

	#training data
	train <- 1:nrow (dataset)
	
	#test data
	prediction <- 1:nrow (candidates)

	#Target variable
	target  <- "Candidate"
	
		
	# Build the Ada Boost model.
	ada <- ada(Candidate ~ .,
	  data=dataset[train,c(input,target)],
	  control=rpart.control(maxdepth=30,
		   cp=0.010000,
		   minsplit=20,
		   xval=10),
			iter=50)

	# Make a prediction using the adaptive boosting model
	pr_class <- predict(ada,(candidates[prediction, c(input, target)]),type="prob")
	
		
	# Extract the relevant variables from the dataset.
	sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))

	# Output the combined data.
	output <- cbind(sdata,pr_class)

		
}

