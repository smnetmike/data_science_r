rm(list = ls())
library(fields)
library(lattice)

setwd("C:/Users/netmike/Documents/Training/data_science/data_science_in_R/geolocation")


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

ReadData = function(csvfile, subMacs = NULL) {
  
  txt = readLines(csvfile)
  
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
  if (is.null(subMacs)) {
    subMacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
  }
  offline = offline[offline$mac %in% subMacs,]

  # Getting rid of channels
  offline = offline["channel" != names(offline)]
  
}

offline = ReadData("offline.final.trace.txt")

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
bwplot(signal ~ factor(angle) | mac, data = offline, subset = posX == 2 & posY == 12 &
         mac != "00:0f:a3:39:dd:cd", layout = c(2,3))
densityplot(~ signal | mac + factor(angle), data = offline, subset = posX == 24 & posY == 4 &
              mac != "00:0f:a3:39:dd:cd", bw = 0.5, plot.points = FALSE)

### Signal summary
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

# Box plot of SD of signal vs average signal
breaks = seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
       data = offlineSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")

# smoothing signal
with(offlineSummary, smoothScatter((avgSignal - medSignal) ~ num), 
     xlab = "number of observations", ylab = "mean - median")
abline(h = 0, col = "#989ea3", lwd = 2)

# Other smoothing
lo.obj = with(offlineSummary, loess(diff ~ num, data.frame(diff = (avgSignal - medSignal), 
                                                  num = num)))
lo.obj.pr = predict(lo.obj, data.frame(num = (70:120)))
lines(x = (70:120), y = lo.obj.pr, col = "#4daf4a", lwd = 2)

### Storing mac adresses
submacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]

# Plotting a signal surface
oneAPAngle = subset(offlineSummary, mac == submacs[1] & angle == 0)
smoothSS = Tps(oneAPAngle[,c("posX", "posY")], oneAPAngle$avgSignal)
vizSmooth = predictSurface(smoothSS)
plot.surface(vizSmooth, type = "C")
points(oneAPAngle[,c("posX", "posY")], pch = 19, cex = 0.5)

# Group of surface plots
surfaceSS <- function(data, m, a) {
  cat("mac=",m,", angle=", a)
  oneAPAngle = subset(data, mac == m & angle == a)
  smoothSS = Tps(oneAPAngle[,c("posX", "posY")], oneAPAngle$avgSignal)
  vizSmooth = predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C")
  points(oneAPAngle[,c("posX", "posY")], pch = 19, cex = 0.5)
}

parCur = par(mfrow = c(2,2), mar = rep(1, 4))
mapply(surfaceSS, m = submacs[ rep(c(1, 2), each = 2) ],
       a = rep(c(0, 135), 2),
       data = list(data = offlineSummary))
par(parCur)

### Computing the distance
offlineSummary = subset(offlineSummary, mac != submacs[2])
AP = matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,
                1, 14, 33.5, 9.3, 33.5, 2.8),
             ncol = 2, byrow = TRUE, dimnames = list(submacs[-2], c("x","y")))
diffs = offlineSummary[, c("posX", "posY")] - AP[offlineSummary$mac, ]
offlineSummary$dist = sqrt(diffs[,1]^2 + diffs[,2]^2)

xyplot(signal ~ dist | factor(mac) + factor(angle), data = offlineSummary, pch = 19, 
       cex = 0.3, xlab = "distance")

### Using the online data
macs = unique(offlineSummary$mac)
online = ReadData("online.final.trace.txt", subMacs = macs)
online$posXY = paste(online$posX, online$posY, "-")

# Some properties
length(unique(online$posXY))
tabonlineXYA = table(online$posXY, online$angle)
tabonlineXYA = tabonlineXYA[1:6,]

### Creating the online summary
keepVars = c("posXY", "posX", "posY", "orientation", "angle")
byLoc = with(online, by(online, list(posXY), function(x) {
   ans = x[1, keepVars]
   avgSS = tapply(online$signal, online$mac, mean)
   y = matrix(avgSS, nrow =1, ncol = 6, dimnames = list(ans$posXY, names(avgSS)))
   cbind(ans, y)
} ))
onlineSummary = do.call("rbind", byLoc)

### Select training data for given angle
reshapeSS = function(data, varSignal = "signal",
                     keepVars = c("posXY", "posX","posY"), sampleAngle = FALSE) {
  byLocation =
    with(data, by(data, list(posXY),
                  function(x) {
                    if (sampleAngle) x = x[x$angle == sample(refs, size = 1), ]
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6,
                               dimnames = list(ans$posXY,
                                               names(avgSS)))
                    cbind(ans, y)
                  }))
  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}
selectTrain <- function(angleNewObs, signals, m) {
  
  nearestAngle = roundOrientation(angleNewObs)
  if (m %% 2 == 1) {
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  } else {
    m = m + 1
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    if (sign(angleNewObs - nearestAngle) > -1)
      angles = angles[ -1 ]
    else
      angles = angles[ -m ]
  }
  
  angles = angles + nearestAngle
  angles[angles < 0] = angles[ angles < 0 ] + 360
  angles[angles > 360] = angles[ angles > 360 ] - 360
  
  offlineSubset =
    offlineSummary[ signals$angle %in% angles, ]
  trainSS = reshapeSS(offlineSubset, varSignal = "avgSignal")
  
  trainSS
}
train130 = selectTrain(130, offlineSummary, m = 3)

### Find nearest neighbour
findNN = function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1,
                function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(trainSubset[closest, 1:3 ])
}

#
undebug(findNN)
nn = findNN(onlineSummary[1,6:11], train130)

### Predicting locations
predXY = function(newSignals, newAngles, trainData,
                  numAngles = 1, k = 3) {
  closeXY = list(length = nrow(newSignals))
  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] =
      findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }
  
  estXY = lapply(closeXY, function(x) sapply(x[ ,2:3], function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

debug(predXY)
options(error = recover, warn = 1)
estXYk3 = predXY(newSignals = onlineSummary[ , 6:11],
                 newAngles = onlineSummary[ , 4],
                 offlineSummary, numAngles = 3, k = 3)
estXYk1 = predXY(newSignals = onlineSummary[ , 6:11],
                 newAngles = onlineSummary[ , 4],
                 offlineSummary, numAngles = 3, k = 1)

# Calculating errors
calcError =
  function(estXY, actualXY)
    sum( rowSums( (estXY - actualXY)^2) )
actualXY = onlineSummary[, c("posX", "posY")]
sapply(list(estXYk1, estXYk3), calcError, actualXY)

# Cross-Validation
keepVars = c("posXY", "posX","posY", "orientation", "angle")
onlineCVSummary = reshapeSS(offline, keepVars = keepVars,
                            sampleAngle = TRUE)

v = 11
permuteLocs = sample(unique(offlineSummary$posXY))
permuteLocs = matrix(permuteLocs, ncol = v,
                     nrow = floor(length(permuteLocs)/v))
K = 20
err = rep(0, K)
for (j in 1:v) {
  onlineFold = subset(onlineCVSummary,
                      posXY %in% permuteLocs[ , j])
  offlineFold = subset(offlineSummary,
                       posXY %in% permuteLocs[ , -j])
  actualFold = onlineFold[ , c("posX", "posY")]
  for (k in 1:K) {
    estFold = predXY(newSignals = onlineFold[ , 6:11],
                     newAngles = onlineFold[ , 4],
                     offlineFold, numAngles = 3, k = k)
    err[k] = err[k] + calcError(estFold, actualFold)
  }
}