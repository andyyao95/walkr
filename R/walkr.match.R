
walkr.match <- function(data,
                        treat  = "treat",
                        match.var,
                        Y.var  = "outcome",
                        method = "MEAN",
                        size = 10,
                        chains = 1,
                        thin   = 1,
                        burn   = .5,
                        breaks = NULL){
  
  # also need to do multiple match.var cases
  
  stopifnot(is.data.frame(data))
  stopifnot(method %in% c("MEAN", "MSQ"))
  
  # treatment vs control groups  
  tgroup <- data[which(data[[treat]] == 1), ]
  cgroup <- data[which(data[[treat]] == 0), ]
  
  # sizes of treatment and control groups
  nt <- nrow(tgroup); nc <- nrow(cgroup)
  
  # preliminary weights from the treatment group (all equal)
  wt <- rep(1/nt, nt)
  
  
  #drop the control group subjects that lie outside the convex hull
  library(WhatIf)
  my.result <- whatif(data = matrix(tgroup[[match.var]], ncol = 1),
                      cfact = matrix(cgroup[[match.var]], ncol = 1))
  
  

  
  # run Walkr to sample weights
  if(method == "MEAN"){
    
    A <- t(matrix(cgroup[[match.var]][my.result$in.hull]))
    b <- mean(tgroup[[match.var]])
    cweight <- rep(0, nc)
    cweight[my.result$in.hull] <- nrow(tgroup)*rowMeans(walkr(A = A, b = b, points = size, 
                                                              chains = chains, thin = thin, 
                                                              burn = burn))

   
    
  } else if(method == "MSQ"){
    
    
    A <- rbind(matrix(cgroup[[match.var]][my.result$in.hull], ncol = nc),
               matrix(cgroup[[match.var]][my.result$in.hull]^2, ncol = nc),
                 rep(1, nc))

    b <- c(mean(tgroup[[match.var]]), 
           mean(tgroup[[match.var]]^2),
           1)
    
      
    cweight[my.result$in.hull] <- nrow(tgroup)*rowMeans(walkr(A = A, b = b, points = size, 
                                                              chains = chains, thin = thin, 
                                                              burn = burn))
    
  }
  
  # weight of treated subjects are all 1
  tweight <- rep(1, nt)
  
  # calculating L1 score
  L1 <- L1.meas(group = c(rep(1, nt), rep(0, nc)),
                # data is a hack, need L1 var of some sort
                data = data[, -c(1, 12)],
                breaks = breaks,
                weights = c(tweight, cweight))$L1
  
  #Setting up return objects
  
  # ??
  if(length(match.var) == 1){
    dat <- data.frame(X = c(tgroup[[match.var]], cgroup[[match.var]]),
                      Weight = c(tweight, cweight),
                      Group = c(rep("Treat", nrow(tgroup)),
                                rep("Control", nrow(cgroup))))
    
    histplot <- ggplot(dat, aes(x = X, weight = Weight, fill = Group)) + 
      geom_histogram(position = "dodge", bins = 30)  +
      xlab(match.var) +
      ylab("Weight")
  }
  
  TE <- weighted.mean(tgroup[[Y.var]], tweight) -
    weighted.mean(cgroup[[Y.var]], cweight)
  
  return(list(data = data,
              tgroup = tgroup, cgroup = cgroup,
              tweight = tweight, cweight = cweight,
              L1 = L1,
              histplot = histplot,
              TE = TE))
}