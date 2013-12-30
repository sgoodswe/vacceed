### Function calls

svm_runPred  <- function (no_of_repeats,output_dir)
{
		
		cat("Running svm_runPred\n")
	
		file_name <- paste(sep="/", output_dir, "svm_predictions.txt")
		
		library (kernlab)
		sink(file = file_name)

		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- svm_makePred ()
			
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
