### Function calls


#### make Prediction function  ####
knn_makePred <- function ()
{
				
	#training data
	train <- dataset [,input]


	#Target variable
	target <- "Candidate"

	#test data
	prediction <- candidates [,input]
	
			
	cl <- factor (dataset [,target])	

	knn_model <- knn(train, prediction, cl, k = 3, prob=TRUE)

	# Get the Class probabilities from the KNN model
	prob_attr <- attr (knn_model,"prob")
						

	#class
	no_of_elements <- length (knn_model)
	pr_class <- knn_model [1:no_of_elements]
	combined <- cbind(pr_class,prob_attr)
	
	# Extract the relevant variables from the dataset.
	prediction <- 1:nrow (candidates)
	sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))

	# Output the combined data.
	output <- cbind(sdata,combined)
	
}


