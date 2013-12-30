### Function calls

rf_runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running rf_runPred\n")
	
		file_name <- paste(sep="/", output_dir, "rf_predictions.txt")
		
		require(randomForest, quietly=TRUE)
		sink(file = file_name)
		
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- rf_makePred ()
		
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
