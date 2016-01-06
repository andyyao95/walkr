L1.single <- function(treat,
                      data,
                      breaks,
                      weights){
  
  tdata <- data[as.logical(treat)]
  cdata <- data[as.logical(1 - treat)]
  cweights <- weights[as.logical(1 - treat)]
  
  # treatment frequency for each basket
  tfreq <- vector()
  tfreq[1] <- length(which(tdata <= breaks[1]))
  for(i in 1:(length(breaks) - 1)){
    tfreq[i + 1] <- length(intersect(which(tdata > breaks[i]),
                                 which(tdata <= breaks[i + 1])))
  }
  tfreq[length(breaks)] <-
    length(tdata[which(tdata > breaks[length(breaks)])])
  tfreq <- tfreq/length(tdata)
  
  # control frequency for each basket
  cfreq <- vector()
  cfreq[1] <- sum(cweights[which(cdata <= breaks[1])])
  for(i in 1:(length(breaks) - 1)){
    cfreq[i + 1] <- sum(cweights[intersect(which(cdata > breaks[i]),
                                 which(cdata <= breaks[i + 1]))])
  }
  cfreq[length(breaks)] <-
    sum(cweights[cdata[which(cdata > breaks[length(breaks)])]])
  cfreq <- cfreq/length(tdata)
  
  # L1 imbalance
  return(sum(abs(tfreq - cfreq))/length(breaks))
}