#############################################
## VolcanoPlot
#############################################

#' @export
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust
#' @importFrom ggrepel geom_text_repel
#' @importFrom marray maPalette

VolcanoPlot <- function(data             = data,
                        prob             = "qvalue"   
                        prob.name        = "Q-value" 
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
    
    	finalfile  <- lastfilename
    	processout <- as.matrix(read.table(finalfile, header = TRUE, sep = "\t"))
  	}	
  
  	processout <- rbind(processout, as.matrix(c(" ", " ", "MSstats - VolcanoPlot function", " "), ncol = 1))
  
  
  	## make upper letter
    type <- "VOLCANOPLOT"
  
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
        
        		processout <- rbind(processout, paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse = ", "), sep = " "))
      			write.table(processout, file = finalfile, row.names = FALSE)
      
      			stop(paste("Please check labels of comparions. Result does not have this comparison. -", paste(temp.name, collapse = ", "), sep = " "))
      		}
    	}
     
    	## check which.comparison is order number of comparison
    	if (is.numeric(which.Comparison)) {
      
    		temp.name <- levels(data$Label)[which.Comparison]
      
      		## message if name of comparison is wrong.
      		if (length(levels(data$Label)) < max(which.Comparison)) {
        		stop(paste("Please check your selection of comparisons. There are ", length(levels(data$Label)), " comparisons in this result.", sep = " "))
        	}
    	}  
    
   	 	## use only assigned proteins
    	data <- data[which(data$Label %in% temp.name), ]
    
   	 	data$Protein <- factor(data$Protein)
    	data$Label   <- factor(data$Label)
  	} else {
    
    	data$Protein <- factor(data$Protein)
    	data$Label   <- factor(data$Label)
  	}
  
  #######################
  ## VolcanoPlot
  #######################
    
    	## If there are the file with the same name, add next numbering at the end of file name		
    	if (address != FALSE) {
      		allfiles <- list.files()
      
      		num <- 0
      		filenaming <- paste(address, "VolcanoPlot",     sep="")
      		finalfile  <- paste(address, "VolcanoPlot.pdf", sep="")
      
      		while(is.element(finalfile, allfiles)) {
        		num       <- num + 1
        		finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
      		}	
      
      		pdf(finalfile, width = width, height = height)
    	}
    
    	if (logBase.pvalue == 2) {
      		y.limUp  <- 30
    	} else if (logBase.pvalue == 10) {
     	 	  y.limUp  <- 10
    	}
    	if (is.numeric(ylimUp)) y.limUp <- ylimUp 
    
    	## remove the result, NA
    	data <- data[!is.na(data[, prob]),]
    
    	## group for coloring dots
      fc.cutoff <- grepl("log[12][0]?FC", colnames(data), value = T)
      if(strsplit(fc.cutoff, split = "")[[1]][4] == 2 {
        fc.cutoff.base <- 2  
      } else if (strsplit(fc.cutoff, split = "")[[1]][4] == 1 &
                 strsplit(fc.cutoff, split = "")[[1]][5] == 10){
        fc.cutoff.base <- 10
        } else {
        processout <- rbind(processout, "Please check labels of Fold Changes. FC should be expressed as log2FC or log10FC")
      	write.table(processout, file = finalfile, row.names = FALSE)
      	stop("Please check labels of Fold Changes. FC should be expressed as log2FC or log10FC")
      }
    
    	if (is.numeric(FCcutoff) & FCcutoff > 0) {
      		data$colgroup <- "black"
      
        		data[data[, prob]  < sig & data[, fc.cutoff] >  log(FCcutoff, base = fc.cutoff.base), "colgroup"] <- "red"
        		data[data[, prob]  < sig & data[, fc.cutoff] < -log(FCcutoff, base = fc.cutoff.base), "colgroup"] <- "blue"
          	} else {
          data[data[, prob]  >= sig, "colgroup"] <- "black"
      		data[data[, prob]  < sig & data[, grepl("log[12][0]?FC", colnames(data))] > 0, "colgroup"] <- "red"    #
      		data[data[, prob]  < sig & data[, grepl("log[12][0]?FC", colnames(data))] < 0, "colgroup"] <- "blue" 
      }
    
     	data$colgroup <- factor(data$colgroup, levels=c("black", "blue", "red"))
    
    	## for multiple volcano plots, 
    	for(i in 1:nlevels(data$Label)) {
      
      		sub <- data[data$Label == levels(data$Label)[i], ]
      
      		if (logBase.pvalue == 2) {
        		sub[, prob] [sub[, prob]  < 2^(-y.limUp)] <- 2^(-y.limUp)
      		}
      
      		if (logBase.pvalue == 10) {
        		sub[, prob] [sub[, prob]  < 10^(-y.limUp)] <- 10^(-y.limUp)
      		}
      
      		sub <- as.data.frame(sub)
      
      		## ylimUp
      		if (logBase.pvalue == 2) {
        		y.limup <- ceiling(max(-log2(sub[!is.na(sub[, prob] ), prob])))
        		if (y.limup < (-log2(sig))) {
        			y.limup <- (-log2(sig) + 1) ## for too small y.lim
        		}
      		}
      
     	 	if (logBase.pvalue == 10) {
        		y.limup <- ceiling(max(-log10(sub[!is.na(sub[, prob] ), prob])))
        		if (y.limup < (-log10(sig))) {
        			y.limup <- (-log10(sig) + 1) ## for too small y.lim
        		}
      		}
       
      		## ylimDown
      		y.limdown <- 0 ## default is zero
      		if (is.numeric(ylimDown)) {
      			y.limdown <- ylimDown
      		}
      
      		## x.lim
      		x.lim <- ceiling(max(abs(sub[!is.na(sub[, 3]) & abs(sub[, 3]) != Inf , 3]))) ## log2FC or log10FC
      		if (x.lim < 3) {
      			x.lim <- 3
      		}
      		if (is.numeric(xlimUp)) {
      			x.lim <- xlimUp
      		}
      
      		## for assigning x in ggplot2
      		subtemp <- sub
      		colnames(subtemp)[3] <- "logFC"
      
      		if (logBase.pvalue == 2) {
        		subtemp$log2adjp <- (-log2(subtemp[, prob] ))
      		}
      
      		if (logBase.pvalue == 10) {
        		subtemp$log10adjp <- (-log10(subtemp[, prob] ))
        	}
        	
        	## for x limit for inf or -inf
        	subtemp$newlogFC <- subtemp$logFC
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing" & subtemp$logFC == Inf, "newlogFC"] <- (x.lim - 0.2)
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing" & subtemp$logFC == (-Inf), "newlogFC"] <- (x.lim - 0.2) *(-1)
        	
        	## add (*) in Protein name for Inf or -Inf
        	subtemp$Protein <- as.character(subtemp$Protein)
        	subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing", "Protein"] <- paste("*", subtemp[!is.na(subtemp$issue) & subtemp$issue == "oneConditionMissing", "Protein"], sep="")

      
      		## Plotting
      		if (logBase.pvalue == 2) {
        		ptemp <- ggplot(aes_string(x='logFC', y='log2adjp', color='colgroup', label='Protein'), data=subtemp)+
        				geom_point(size=dot.size)+
        				scale_colour_manual(values=c("gray65", "blue", "red"), 
        				                    limits=c("black", "blue", "red"), 
        				                    breaks=c("black", "blue", "red"), 
        				                    labels=c("No regulation", "Down-regulated", "Up-regulated"))+
        				scale_y_continuous('-Log2 (adjusted p-value)', 
        				                   limits=c(y.limdown, y.limup))+
        				labs(title=unique(sub$Label))
      		}
      
      		if (logBase.pvalue == 10) {
        		ptemp <- ggplot(aes_string(x='logFC', y='log10adjp', color='colgroup', label='Protein'), data=subtemp)+
        				geom_point(size=dot.size)+
        				scale_colour_manual(values=c("gray65", "blue", "red"), 
        				                    limits=c("black", "blue", "red"), 
        				                    breaks=c("black", "blue", "red"), 
        				                    labels=c("No regulation", "Down-regulated", "Up-regulated"))+
        				scale_y_continuous('-Log10 (adjusted p-value)', 
        				                   limits=c(y.limdown, y.limup))+
        				labs(title=unique(sub$Label))
        	}
      
      
     		## x-axis labeling
      		if (colnames(sub)[3] == "log2FC") {
      			ptemp <- ptemp+scale_x_continuous('Log2 fold change', limits=c(-x.lim, x.lim))
      		}
      		if (colnames(sub)[3] == "log10FC") {
      			ptemp <- ptemp+scale_x_continuous('Log10 fold change', limits=c(-x.lim, x.lim))
      		}
      
     		## add protein name
      		if (ProteinName) {
      			if(length(unique(subtemp$colgroup)) == 1 & any(unique(subtemp$colgroup) == 'black')){
      				message(paste("The volcano plot for ", unique(subtemp$Label), " does not show the protein names because none of them is significant.", sep=""))
      				
      			} else {
          			
          			ptemp <- ptemp + geom_text_repel(data=subtemp[subtemp$colgroup != "black", ], aes(label=Protein), size=text.size, col='black')
          		}
      		} 
      
      		## For legend of linetype for cutoffs
      		## first assign line type
      		ltypes <- c("type1"="twodash", "type2"="dotted")
      
      		## cutoff lines, FDR only
     		if (!FCcutoff) { 
        		if (logBase.pvalue == 2) {
          			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=20), log2adjp=(-log2(sig)), line='twodash')
          
          			pfinal <- ptemp + geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), 
          			                            colour="darkgrey", 
          			                            size=0.6, 
          			                            show.legend=TRUE)+
            			    scale_linetype_manual(values=c('twodash'=6), 
            			                          labels=c(paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
            			    		
       	 		}
        
        		if (logBase.pvalue == 10) {
          			sigcut <- data.frame(Protein='sigline', logFC=seq(-x.lim, x.lim, length.out=20), log10adjp=(-log10(sig)), line='twodash')
          
         	 		pfinal <- ptemp + geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), 
         	 		                            colour="darkgrey", 
         	 		                            size=0.6, 
         	 		                            show.legend=TRUE)+
            			    scale_linetype_manual(values=c('twodash'=6), 
            			                          labels=c(paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
				}				
      		}
      
      		## cutoff lines, FDR and Fold change cutoff
      		if (is.numeric(FCcutoff)) {
        		if (colnames(sub)[3] == "log2FC") {
          			if (logBase.pvalue == 2) {
            
           				## three different lines
            			sigcut <- data.frame(Protein='sigline', 
            			                     logFC=seq(-x.lim, x.lim, length.out=10), 
            			                     log2adjp=(-log2(sig)), 
            			                     line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', 
            			                       logFC=log2(FCcutoff), 
            			                       log2adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', 
            			                       logFC=(-log2(FCcutoff)), 
            			                       log2adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            			                          colour="darkgrey", 
            			                          size=0.6, 
            			                          show.legend=TRUE)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6, show.legend=TRUE)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), 
            			                          labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}
          
          			if (logBase.pvalue == 10) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', 
            			                     logFC=seq(-x.lim, x.lim, length.out=10), 
            			                     log10adjp=(-log10(sig)), line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', 
            			                       logFC=log2(FCcutoff), 
            			                       log10adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', 
            			                       logFC=(-log2(FCcutoff)), 
            			                       log10adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            			                          colour="darkgrey", 
            			                          size=0.6, 
            			                          show.legend=TRUE)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6, 
            				          show.legend=TRUE)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), 
            			                          labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}        
        		}
        
        		if (colnames(sub)[3] == "log10FC") {
          			if (logBase.pvalue == 2) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', 
            			                     logFC=seq(-x.lim, x.lim, length.out=10), 
            			                     log2adjp=(-log2(sig)), 
            			                     line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', 
            			                       logFC=log10(FCcutoff), 
            			                       log2adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', 
            			                       logFC=(-log10(FCcutoff)), 
            			                       log2adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            			                          colour="darkgrey", 
            			                          size=0.6, 
            			                          show.legend=TRUE)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6, 
            				          show.legend=TRUE)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log2adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), 
            			                          labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}
          
          			if (logBase.pvalue == 10) {
            
            			## three different lines
            			sigcut <- data.frame(Protein='sigline', 
            			                     logFC=seq(-x.lim, x.lim, length.out=10), 
            			                     log10adjp=(-log10(sig)), 
            			                     line='twodash')
            			FCcutpos <- data.frame(Protein='sigline', 
            			                       logFC=log10(FCcutoff), 
            			                       log10adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            			FCcutneg <- data.frame(Protein='sigline', 
            			                       logFC=(-log10(FCcutoff)), 
            			                       log10adjp=seq(y.limdown, y.limup, length.out=10), 
            			                       line='dotted')
            
            			## three lines, with order color first and then assign linetype manual
            			pfinal <- ptemp+geom_line(data=sigcut, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            			                          colour="darkgrey", 
            			                          size=0.6, 
            			                          show.legend=TRUE)+
            				geom_line(data=FCcutpos, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6, 
            				          show.legend=TRUE)+
            				geom_line(data=FCcutneg, aes_string(x='logFC', y='log10adjp', linetype='line'), 
            				          colour="darkgrey", 
            				          size=0.6)+
            			    scale_linetype_manual(values=c('dotted'=3, 'twodash'=6), 
            			                          labels=c(paste("Fold change cutoff (", FCcutoff, ")", sep=""), paste("Adj p-value cutoff (", sig, ")", sep="")))+
            			    guides(colour=guide_legend(override.aes=list(linetype=0)),
            			    		linetype=guide_legend())
          			}    
        		}
     		 }
      
      		pfinal <- pfinal+theme(
        		panel.background = element_rect(fill='white', colour="black"),
        		panel.grid.minor = element_blank(),
        		axis.text.x = element_text(size=x.axis.size, colour="black"),
        		axis.text.y = element_text(size=y.axis.size, colour="black"),
        		axis.ticks = element_line(colour="black"),
        		axis.title.x = element_text(size=x.axis.size+5, vjust=-0.4),
        		axis.title.y = element_text(size=y.axis.size+5, vjust=0.3),
        		title = element_text(size=x.axis.size+8, vjust=1.5),
        		legend.position="bottom",
        		legend.key = element_rect(fill='white', colour='white'),
        		legend.text = element_text(size=legend.size),
        		legend.title = element_blank()
        		)
      
      		print(pfinal)
    	} ## end-loop
    
    	if (address!=FALSE) dev.off()
}
