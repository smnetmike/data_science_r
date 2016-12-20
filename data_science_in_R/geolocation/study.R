rm(list = ls())

setwd("C:/Users/netmike/Documents/Training/data_science/data_science_in_R/geolocation")

txt = readLines("offline.final.trace.txt")

tokens = strsplit(txt[4], "[;=,]")[[1]]

tokens[c(2, 4, 6:8, 10)]
matrix(tokens[-(1:10)], ncol=4, byrow = TRUE)


processLine = function (x) {
  tokens = strsplit(x, "[;=,]")[[1]]
  if (length(tokens) == 10)
    return(NULL)
  
  tmp = matrix(tokens[-(1:10)], ncol=4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
}

lines = txt[ substr(txt, 1, 1) != "#" ]
tmp <- lapply(lines, processLine)
offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)

dim(offline)