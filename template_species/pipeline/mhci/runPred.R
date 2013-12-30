### Function calls

runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running runPred\n")
	
		file_name <- paste(sep="/", output_dir, "predictions.txt")
		
		require(randomForest, quietly=TRUE)
		sink(file = file_name)
		
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- makePred ()
		
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
