#############################################
## VolcanoPlot
#############################################

#' @export
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust
#' @importFrom ggrepel geom_text_repel
#' @importFrom marray maPalette

VolcanoPlot <- function(data             = data,
                        type             = type,
                        sig              = 0.05,
                        FCcutoff         = FALSE,
                        logBase.pvalue   = 10,
                        ylimUp           = FALSE,
                        ylimDown         = FALSE,
                        xlimUp           = FALSE,
                        x.axis.size      = 10,
                        y.axis.size      = 10,
                        dot.size         = 3,
                        text.size        = 4, 
                        legend.size      = 13,
                        ProteinName      = TRUE,
                        numProtein       = 100, 
                        clustering       = "both", 
                        width            = 10, 
                        height           = 10, 
                        which.Comparison = "all", 
                        address="") {
                        
    ## save process output in each step
    allfiles   <- list.files()
  	filenaming <- "msstats"
  
  	if (length(grep(filenaming, allfiles)) == 0) {
    
    	finalfile  <- "msstats.log"
    	processout <- NULL
    
  	} else {
    
    	num <- 0
    	finalfile <- "msstats.log"
    
    	while(is.element(finalfile, allfiles)) {
      		num <- num + 1
      		lastfilename <- finalfile ## in order to rea
      		finalfile    <- paste0(paste(filenaming, num, sep="-"), ".log")
    	}
    
    	finalfile <- lastfilename
    	processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
  	}	
  
  	processout <- rbind(processout, as.matrix(c(" ", " ", "MSstats - groupComparisonPlots function", " "), ncol=1))
  
  
  	## make upper letter
  	type <- toupper(type)
  
  	if (length(setdiff(type, c("HEATMAP", "VOLCANOPLOT", "COMPARISONPLOT"))) != 0) {
    
    	processout <- rbind(processout, c(paste("Input for type=", type, ". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".", sep="")))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop(paste("Input for type=", type, ". However,'type' should be one of \"Heatmap\", \"VolcanoPlot\",\"ComparisonPlot\".", sep=""))
  	}
  
  	## check logBase.pvalue is 2,10 or not
  	if (logBase.pvalue != 2 & logBase.pvalue != 10) {
    	processout <- rbind(processout, c("ERROR : (-) Logarithm transformation for adjusted p-values : log2 or log10 only - stop"))
    	write.table(processout, file=finalfile, row.names=FALSE)
    
    	stop("Only -log2 or -log10 for logarithm transformation for adjusted p-values are posssible.\n")
  	}
  
  	## choose comparison to draw plots
  
  	if (which.Comparison != "all") {
    	## check which.comparison is name of comparison
    	if (is.character(which.Comparison)) {
      
      		temp.name <- which.Comparison
      
      		## message if name of comparison is wrong.
      		if (length(setdiff(temp.name, unique(data$Label))) > 0) {
        
        		processout <- rbind(processout, paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "), sep=" "))
      			write.table(processout, file=finalfile, row.names=FALSE)
      
      			stop(paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse=", "), sep=" "))
      		}
    	}
     
    	## check which.comparison is order number of comparison
    	if (is.numeric(which.Comparison)) {
      
    		temp.name <- levels(data$Label)[which.Comparison]
      
      		## message if name of comparison is wrong.
      		if (length(levels(data$Label))<max(which.Comparison)) {
        		stop(paste("Please check your selection of comparisons. There are ", length(levels(data$Label)), " comparisons in this result.", sep=" "))
        	}
    	}  
    
   	 	## use only assigned proteins
    	data <- data[which(data$Label %in% temp.name), ]
    
   	 	data$Protein <- factor(data$Protein)
    	data$Label <- factor(data$Label)
  	} else {
    
    	data$Protein <- factor(data$Protein)
    	data$Label <- factor(data$Label)
  	}
                        }
