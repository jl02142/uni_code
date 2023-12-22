# J Ludwig 2018

### File should be output of samtools depth.  File name should be sample name without extensions (e.g. .txt). 
### File should be saved within a folder named after file.  
      # i.e. - Sample.name.directory/sample.name.samtools.output or e.g. P60002/P60002

if(file.exists(".Rdata")){unlink(".Rdata")}
rm(list = ls())
library(data.table)
args <- commandArgs(TRUE)
options(scipen = 999)

### If need code to install packages ###
# new.packages <- "data.table"[!("data.table" %in% installed.packages()[,"Package"])]
# if(length(new.packages)){
#   install.packages(new.packages)
#   library(data.table)
# }

path <- getwd()

# Option -i must be followed by a single column file with the directory name that contains the depth file.  
# Depth file must be unaltered from samtools output.  Samtools V 1.3.1 and V 1.6 depth confirmed to work with this code.  
# R V 3.5.0 and data.table V 1.11.4 confirmed to work. 
vector.of.options <- c("-i", "-o")

vector.of.args <- vector("list", length(args)/2)
if(length(args) > 0){
  for(k in 1:length(args)){
    if(k %% 2 == 1){ #if k == odd
      vector.of.args[[k]][1] <- args[k]
      vector.of.args[[k]][2] <- args[k+1]
    } else {
      next
    }
  }
}

if(length(args) == 0){ 
  stop("Error, at least one argument must be supplied")
} else if(is.logical(grep("-h|--help", args))){
  stop("Options listed here")
} else for(i in 1:length(vector.of.args)){
  if(is.na(charmatch(vector.of.args[[i]][1], vector.of.options))){
    stop("Error, unrecognized option supplied, see -h or --help for accepted options")
  } else if(is.na(charmatch(vector.of.args[[i]][1], "-i"))){
    stop("Error; input file, supplied in the form of -i input_file_name, required")
  }
} 

i <- vector.of.args[[grep("-i", vector.of.args, fixed = TRUE)]][2]

#Code below to obtain list of samples and list of chromosomes.
list.of.samples<-fread(i, header = FALSE, col.names = c("Sample.name"))
setwd(paste0(path, "/", list.of.samples$Sample.name[1]))
temp.chr.list.data <- fread(file.path(list.of.samples[1]), drop=2:3)
list.of.chr        <- sapply(temp.chr.list.data, subset, !duplicated(temp.chr.list.data))
list.chr.scale <- list()
for(p in list.of.chr){
  list.chr.scale[p] <- colSums(temp.chr.list.data == p)
}
rm(temp.chr.list.data)
setwd("..")

# Create data table called test.
start.time <- Sys.time()
test <- apply(list.of.samples, 1, function(x) fread(file.path(path, x, x), col.name = c("chromosome", "position", "depth")))
end.time <- Sys.time()
end.time - start.time

# Create single file containing per sample histogram of depth.
# Filtered data set on 99th quantile plotted.
# Mean and median are calculated from unfiltered data set.
start.time <- Sys.time()
pdf(paste0("all_hist.pdf"), width = 8.5, height = 11)
par(mfrow = c(3,3))
for(k in list.of.samples$Sample.name){
  k.num <- which(list.of.samples == k)
  hist.max <- round(test[[k.num]][, unname(quantile(depth, 0.99))])
  cut.value <- round(hist.max/30)
  temp.table <- test[[k.num]][depth < quantile(depth, 0.99)]
  hist(as.numeric(temp.table[, cut(depth, 30)]), xaxt = 'n', ylab = "", xlab = 'Depth (bp)', main = k, cex.main = 1, cex.lab = 1)
  axis(1, at = 0:29, labels = c(round(seq(0, max(temp.table[, depth]), length.out = 30))), line = 0, cex.axis = 0.75)  ### Added because at and labels were coming out of diff lengths by 1
  usr <- par("usr")
  clip(usr[1], usr[2], 0, usr[4])
  abline(v = test[[k.num]][, mean(depth)] / cut.value, col = 'blue', lwd = 1.5)
  abline(v = test[[k.num]][, median(depth)] / cut.value, col = 'red', lwd = 1.5)
  legend("topright", cex = 0.75, legend = c(paste("Mean ", round(test[[k.num]][, mean(depth)])), paste("Median ", round(test[[k.num]][, median(depth)]))), col=c("blue","red",bty="n"), lwd = 1, 
         bg = "white", bty = "n")
}
dev.off()
end.time <- Sys.time()
end.time - start.time

# Individual hist per chromosome code below
# Filtered data set on 99th quantile plotted
# Mean and median are calculated from unfiltered data set
start.time <- Sys.time()
for(k in list.of.samples$Sample.name){
  k.num <- which(list.of.samples == k)
  pdf(paste0(k,".pdf"), width = 8.5, height = 11)
  par(mfrow = c(3,3))
  hist.max <- round(test[[k.num]][, unname(quantile(depth, 0.99))])
  cut.value <- round(hist.max/30)
  temp.table <- test[[k.num]][depth < quantile(depth, 0.99)]
  setkey(temp.table)
  setkey(test[[k.num]])
  for (n in list.of.chr){
    hist(as.numeric(temp.table[n, cut(depth, 30)]), xaxt = 'n', ylab = "", xlab = 'Depth (bp)', main = paste0(k,"_",n), cex.main = 1, cex.lab = 1)
    axis(1, at = 0:29, labels = c(round(seq(0, max(temp.table[, depth]), length.out = 30))), line = 0, cex.axis = 0.75)  # Added because at and labels were coming out of diff lengths by 1
    usr <- par("usr")
    clip(usr[1], usr[2], 0, usr[4])
    abline(v = test[[k.num]][n, mean(depth)] / cut.value, col = 'blue', lwd = 1.5)
    abline(v = test[[k.num]][n, median(depth)] / cut.value, col = 'red', lwd = 1.5)
    legend("topright", cex = 0.75, legend = c(paste("Mean ", round(test[[k.num]][n, mean(depth)])), paste("Median ", round(test[[k.num]][n, median(depth)]))), col=c("blue","red",bty="n"), lwd = 1, 
           bg = "white", bty = "n")
  }
  dev.off()
}
end.time <- Sys.time()
end.time - start.time

# No filter non-overlapping 1kb sliding window.
start.time <- Sys.time()
pdf(paste0("all_window_raw.pdf"), width = 8.5, height = 11)
layout(matrix(1:32, 4, 8, byrow = TRUE), widths = list.chr.scale)
for(k in list.of.samples$Sample.name){
  k.num <- which(list.of.samples == k)
  setkey(test[[k.num]], chromosome)
  for(l in list.of.chr){
    i <- 1
    length.of.chr <- test[[k.num]][.(l), (.N)]
    temp.table <- test[[k.num]][.(l), depth]
    hold <- vector("list", length.of.chr/1000)
    for(m in seq(1, length.of.chr, 1000)){
      if((m+1000) < length.of.chr){ 
        hold[[i]] <- mean(temp.table[m:(m + 1000)])
      } else {
        hold[[i]] <- mean(temp.table[m:length.of.chr])
      }
      i <- i + 1
    }
    if(l == list.of.chr[1]){
      par(mar = c(3, 3, 1, 0))
    } else if (l == list.of.chr[length(list.of.chr)]){
      par(mar = c(3, 0, 1, 1))
    } else {
      par(mar = c(3, 0, 1, 0))
    }
    plot(seq(1, ceiling(length.of.chr), 1000), hold, xlab = '', ylab = '', pch = 20, cex = 0.75, ylim = c(0 , max(test[[k.num]][, max(depth)])), xaxt = 'n', yaxt = 'n')
    abline(h = test[[k.num]][l, mean(depth)], col = 'blue', lwd = 1.5)
    abline(h = test[[k.num]][l, median(depth)], col = 'red', lwd = 1.5)
    abline(h = test[[k.num]][, mean(depth)], col = 'green', lwd = 1.5)
    abline(h = test[[k.num]][, median(depth)], col = 'orange', lwd = 1.5)
    axis.at <- seq(0, length.of.chr, 500000)
    axis(1, at = c(axis.at[2:(length(axis.at))]), cex.axis = 0.75)
    if(l == list.of.chr[1]){
      axis(side = 2, cex.axis = 1)
      legend("topright", cex = 0.75, legend = c(
        paste("Mean ", round(test[[k.num]][l, mean(depth)])), 
        paste("Median ", round(test[[k.num]][l, median(depth)])), 
        paste("Sample mean ", round(test[[k.num]][, mean(depth)])), 
        paste("Sample median ", round(test[[k.num]][, median(depth)]))),
        col = c("blue", "red", "green", "orange"),
        lwd = 1, bg = "white", bty = "n"
      )
      title(l, line = 0.2, cex.main = 0.75)
      mtext(k, 3, adj = 0, cex = 0.75)
    } else {
      legend("topright", cex = 0.75, legend = c(
        paste("Mean ", round(test[[k.num]][l, mean(depth)])), 
        paste("Median ", round(test[[k.num]][l, median(depth)]))),
        #col = c("blue", "red"),
        bg = "white", bty = "n"
      )
      title(l, line = 0.2, cex.main = 0.75)
    }
  }
  mtext("Position (bp)", 1, cex = 0.75, outer = TRUE, line = -1)
  mtext("Mean Read Depth", 2, cex = 0.75, outer = TRUE, line = -1)
}
dev.off()
end.time <- Sys.time()  
end.time - start.time


# Filtered on 5X the median of sample non-overlapping 1kb sliding window. 
# Mean and median calculated from unfiltered data set
Rprof("file.out.2")
start.time <- Sys.time()
pdf(paste0("all_window_filter.pdf"), width = 8.5, height = 11)
layout(matrix(1:32, 4, 8, byrow = TRUE), widths = list.chr.scale)
for(k in list.of.samples$Sample.name){
  k.num <- which(list.of.samples == k)
  temp.data.table <- data.table(test[[k.num]][depth < (median(depth)) * 5])
  setkey(temp.data.table, chromosome)
  for(l in list.of.chr){
    i <- 1
    length.of.chr <- temp.data.table[.(l), (.N)]
    temp.table <- temp.data.table[.(l), depth]
    hold <- vector("list", length.of.chr/1000)
    for(m in seq(1, length.of.chr, 1000)){
      if((m+1000) < length.of.chr){
        hold[[i]] <- mean(temp.table[m:(m+1000)])
      } else {
        hold[[i]] <- mean(temp.table[m:length.of.chr])
      }
      i <- i + 1
    }
    if(l == list.of.chr[1]){
      par(mar = c(3, 3, 1, 0))
    } else if (l == list.of.chr[length(list.of.chr)]){
      par(mar = c(3, 0, 1, 1))
    } else {
      par(mar = c(3, 0, 1, 0))
    }
    plot(seq(1, ceiling(length.of.chr), 1000), hold, xlab = '', ylab = '', pch = 20, cex = 0.75, ylim = c(0 , max(temp.data.table[, max(depth)])), xaxt = 'n', yaxt = 'n')
    #abline(h = test[[k.num]][l, mean(depth)], col = 'blue', lwd =1.5)  ## Comment mean out
    abline(h = test[[k.num]][l, median(depth)], col = 'red', lwd = 1.5)
    #abline(h = test[[k.num]][, mean(depth)], col = 'green', lwd = 1.5)  ## Comment mean out
    #abline(h = test[[k.num]][, median(depth)], col = 'orange', lwd = 1.5) ## Comment sample-wide median out
    axis.at <- seq(0, length.of.chr, 500000)
    axis(1, at = c(axis.at[2:(length(axis.at))]), cex.axis = 0.75)
    if(l == list.of.chr[1]){
      axis(side = 2, cex.axis = 1)
      legend("topright", cex = 0.75, legend = c(
        ## The final paste statement below require a total of 3 ")" ; others require 2
        #paste("Mean ", round(test[[k.num]][l, mean(depth)])), ## Comment mean out
        paste("Median ", round(test[[k.num]][l, median(depth)]))), 
        #paste("Sample mean ", round(test[[k.num]][, mean(depth)])), ## Comment mean out
        #paste("Sample median ", round(test[[k.num]][, median(depth)]))), ## Comment sample-wide median out
        #col = c("blue", "red", "green", "orange"),
        col = "red", ## Comment this and uncomment above line in addition to Mean lines above to see mean. Also uncomment sample-wide median
        lwd = 1, bg = "white", bty = "n"
      )
      title(l, line = 0.2, cex.main = 0.75)
      mtext(k, 3, adj = 0, cex = 0.75)
    } else {
      legend("topright", cex = 0.75, legend = c(
        #paste("Mean ", round(test[[k.num]][l, mean(depth)])), 
        paste("Median ", round(test[[k.num]][l, median(depth)]))),
        #col = c("blue", "red"),
        col = "red",
        bg = "white", bty = "n"
      )
      title(l, line = 0.2, cex.main = 0.75)
    }
  }
  mtext("Position (bp)", 1, cex = 0.75, outer = TRUE, line = -1)
  mtext("Mean Read Depth", 2, cex = 0.75, outer = TRUE, line = -1)
}

dev.off()
end.time <- Sys.time()  
end.time - start.time
Rprof(NULL)


Rprof("file.out.4")
start.time <- Sys.time()

# Create pairwise comparision table. 
for(k in list.of.samples$Sample.name){
  k.num <- which(list.of.samples == k)
  setkey(test[[k.num]], chromosome)
  result <- data.table(matrix(nrow = length(list.of.chr), ncol = length(list.of.chr) + 1))
  setnames(result, c(k, list.of.chr))
  for(l in list.of.chr){
    l.num <- which(list.of.chr == l)
    result[[k]][l.num] <-  l
    for(m in list.of.chr){
      m.num <- which(list.of.chr == m)
      result[[l]][m.num] <- ((test[[k.num]][l, median(depth)]) / test[[k.num]][m, median(depth)])
    }
  }
  fwrite(result, file ="pairwise", append = TRUE, sep ="\t", col.names = TRUE)
}


end.time <- Sys.time()  
end.time - start.time
Rprof(NULL)
