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

roundOrientation = function(angles) {
  refs = seq(0, by = 45, length = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}

ReadData = function(txt) {
  
  # Creating a data frame
  lines = txt[ substr(txt, 1, 1) != "#" ]
  tmp <- lapply(lines, processLine)
  offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)
  
  # Naming the data frame columns
  names(offline) = c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac", "signal", "channel", "type")
  numVars = c("time", "posX", "posY", "posZ", "orientation", "signal")
  offline[ numVars ] = lapply(offline[ numVars ], as.numeric)
  
  # Restricting the data from only access points
  offline = offline[offline$type == "3",]
  offline = offline["type" != names(offline)]
  
  # Changing the time format
  offline$rawtime = offline$time
  offline$time = offline$time/1000
  class(offline$time) = c("POSIXt", "POSIXct")
  
  # Taking some repetitive columns out
  offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]
  
  # Rounding up orientations
  offline$angle = roundOrientation(offline$orientation)
  
  # Selecting the right Access Points
  submacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
  offline = offline[offline$mac %in% submacs,]

  # Getting rid of channels
  offline = offline["channel" != names(offline)]
  
}

offline = ReadData(txt)

#orientation
length(unique(offline$orientation))
plot(ecdf(offline$orientation))

# Drawing boxplots:
with(offline, boxplot(orientation ~ angle, xlab = "nearest 45 degree angle",
 ylab="orientation"))

# Getting rid of worthless positions
locDF = with(offline, by(offline, list(posX, posY), function(x) x))
locDF = locDF[!sapply(locDF, is.null)]

# Displaying the locations
locCounts = sapply(locDF, nrow)
locCounts = sapply(locDF, function(df) c(df[1, c("posX", "posY")], count = nrow(df)))
locCounts = t(locCounts)
plot(locCounts, type = "n", xlab = "", ylab = "")
text(locCounts, labels = locCounts[,3], cex = 0.8, srt = 45)

# Density plots and box plots
library(lattice)
bwplot(signal ~ factor(angle) | mac, data = offline, subset = posX == 2 & posY == 12 &
         mac != "00:0f:a3:39:dd:cd", layout = c(2,3))
densityplot(~ signal | mac + factor(angle), data = offline, subset = posX == 24 & posY == 4 &
              mac != "00:0f:a3:39:dd:cd", bw = 0.5, plot.points = FALSE)

# Signal summary
offline$posXY = paste(offline$posX, offline$posY, sep = "-")
byLocAngleAP = with(offline, by(offline, list(posXY, angle, mac), function(x) x))
signalSummary = lapply(byLocAngleAP, function(oneLoc) {
  summary = oneLoc[1,]
  summary$medSignal = median(oneLoc$signal)
  summary$avgSignal = mean(oneLoc$signal)
  summary$num = length(oneLoc$signal)
  summary$sdSignal = sd(oneLoc$signal)
  summary$iqr = IQR(oneLoc$signal)
  summary
})

offlineSummary = do.call("rbind", signalSummary)

breaks = seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
       data = offlineSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")

# Plotting a signal surface
submacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
oneAPAngle = subset(offlineSummary, mac == submacs[5] & angle == 0)
library(fields)
smoothSS = Tps(oneAPAngle[,c("posX", "posY")], oneAPAngle$avgSignal)
vizSmooth = predictSurface(smoothSS)
plot.surface(vizSmooth, type = "C")