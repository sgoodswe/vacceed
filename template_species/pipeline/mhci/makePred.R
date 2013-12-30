makePred <- function ()
{

		#training data
		train <- 1:nrow (dataset)
		
		#test data
		prediction <- 1:nrow (candidates)
		

		#Target variable
		target  <- "Candidate"

		rf <- randomForest(as.factor(Candidate) ~ .,
			  data=dataset[train,c(input, target)], 
			  ntree=100,
			  mtry=5,
			  importance=TRUE,
			  na.action=na.roughfix,
			  replace=FALSE)
	

		# Make a prediction using the Random Forest model
		pr_class <- predict(rf,(candidates[prediction, c(input, target)]),type="prob")
		
			
		# Extract the relevant variables from the dataset.
		sdata <- subset(candidates[prediction,], select=c("ID","Candidate"))

		# Output the combined data.
		output <- cbind(sdata,pr_class)

}

