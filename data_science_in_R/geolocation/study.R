rm(list = ls())

setwd("C:/Users/netmike/Documents/Training/data_science/data_science_in_R/geolocation")

txt = readLines("offline.final.trace.txt")


processLine = function (x) {
  tokens = strsplit(x, "[;=,]")[[1]]
  if (length(tokens) == 10)
    return(NULL)
  
  tmp = matrix(tokens[-(1:10)], ncol=4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
}

# Creating a data frame
lines = txt[ substr(txt, 1, 1) != "#" ]
tmp <- lapply(lines, processLine)
offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)

# Naming the data frame columns
names(offline) = c("time", "scanMac", "posX", "posY", "posZ",
                   "orientation", "mac", "signal",
                   "channel", "type")
numVars = c("time", "posX", "posY", "posZ",
            "orientation", "signal")
offline[ numVars ] = lapply(offline[ numVars ], as.numeric)

# Restricting the data from only access points
offline = offline[offline$type == "3",]
offline = offline["type" != names(offline)]
dim(offline)

# Changing the time format
offline$rawtime = offline$time
offline$time = offline$time/1000
class(offline$time) = c("POSIXt", "POSIXct")


#Summarizing the numeric data
summary(offline[ numVars ])

summary(sapply(offline[ ,c("mac", "channel", "scanMac")], as.factor))

offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]

#orientation
length(unique(offline$orientation))
plot(ecdf(offline$orientation))

roundOrientation = function(angles) {
  refs = seq(0, by = 45, length = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}
offline$angle = roundOrientation(offline$orientation)

# Drawing boxplots:
with(offline, boxplot(orientation ~ angle, xlab = "nearest 45 degree angle",
 ylab="orientation"))

# Selecting the right mac adresses
submac = names(sort(table(offline$mac), decreasing = TRUE))[1:7]

