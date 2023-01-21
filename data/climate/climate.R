### Beijing climate data set


### DATA PRE-PROCESSING

filenames <- list.files()
filenames <- filenames[grep("csv", filenames)]

stations <- unlist(lapply(sapply(filenames, strsplit, "_"), function(z) z[[3]]))
names(stations) <- rep("", length(stations))

#sapply(filenames, function(z) length(readLines(z)))-1
# each file has 35,064 records
n0 <- 35064

data_combined <- as.data.frame(matrix(nrow = length(filenames)*n0, ncol = 17))
colnames(data_combined) <- colnames(read.csv(filenames[1], nrows = 1, sep = ",")[,-1])
data_combined[,"station"] <- factor(data_combined[,"station"], levels = stations)

for(i in 1:length(filenames)){

    data_combined[((i-1)*n0+1):(i*n0),] <- read.csv(filenames[i], header = TRUE)[,-1]

}

is_c <- apply(data_combined, 1, function(z) !any(is.na(z)))

data_combined_c <- data_combined[is_c,]
row_nam <- as.numeric(rownames(data_combined_c))

write.csv(data_combined, file = "data_all.csv")

                                        #X <- data_combined_c[,c("TEMP","PRES","DEWP","RAIN","WSPM", "CO")]
library(xtable)
xtable(data_combined[c(1,410819),])
