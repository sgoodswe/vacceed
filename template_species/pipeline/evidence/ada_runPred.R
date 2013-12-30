### Function calls

ada_runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running ada_runPred\n")

		require(ada, quietly=TRUE)
	
		file_name <- paste(sep="/", output_dir, "ada_predictions.txt")
			
		sink(file = file_name)
			
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- ada_makePred ()
			
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()
}
