### Function calls


#### makee Prediction function  ####
svm_makePred <- function ()
{

		#training data
		train <- 1:nrow (dataset)
		
		#test data
		prediction <- 1:nrow (candidates)
	
		#Target variable
		target  <- "Candidate"
	
		svm <- ksvm ( formula (Candidate ~ .),
				data=dataset[train,c(input, target)], 
				kernel = "besseldot",
				type = "C-svc",
				prob.model=TRUE)
		
			
		# Make a prediction using the SVM model
		pr_class <- predict(svm, na.omit (candidates[prediction, input]), type="prob")
				
		
		# Extract the relevant variables from the dataset.
		sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))

		# Output the combined data.
		output <- cbind(sdata,pr_class)
	
		
}

