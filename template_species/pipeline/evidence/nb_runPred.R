### Function calls

nb_runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running nb_runPred\n")

		library(e1071)
	
		file_name <- paste(sep="/", output_dir, "nb_predictions.txt")
			
		sink(file = file_name)
		
		
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- nb_makePred ()
			
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
