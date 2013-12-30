### Function calls

nn_runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running nn_runPred\n")

		require(nnet, quietly=TRUE)
	
		file_name <- paste(sep="/", output_dir, "nn_predictions.txt")
			
		sink(file = file_name)
	
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- nn_makePred ()
			
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
