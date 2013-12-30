### Function calls

knn_runPred <- function (no_of_repeats,output_dir)
{
		
		cat("Running knn_runPred\n")

		library ("class")
	
		file_name <- paste(sep="/", output_dir, "knn_predictions.txt")
			
		sink(file = file_name)
		
		increment <- 1
		while( increment <= no_of_repeats) { 

			out <- knn_makePred ()
			
			cat("Run:",increment,"\n")
	
			print (out)
			
			increment <- increment + 1
		}
		
		sink ()

}
