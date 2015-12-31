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
  
  # run Walkr to sample weights
  if(method == "MEAN"){
    weight <- walkr(rbind(matrix(cgroup[[match.var]], ncol = nc), rep(1, nc)),
                    c(t(wt)%*%tgroup[[match.var]], 1),
                    points = size,
                    chains = chains,
                    thin = thin,
                    burn = burn)
    weight <- rowMeans(weight)
  } else if(method == "MSQ"){
    
    A <- rbind(matrix(cgroup[[match.var]], ncol = nc),
               matrix(cgroup[[match.var]]^2, ncol = nc),
               rep(1, nc))

    b <- c(t(wt)%*%tgroup[[match.var]],
           t(wt)%*%(tgroup[[match.var]]^2),
           1)
    
    weight <- walkr(A,
                    b,
                    points = size,
                    chains = chains,
                    thin = thin,
                    burn = burn)
  }
  
  # take the mean and scale up
  cweight <- rowMeans(weight) * nt; tweight <- rep(1, nt)
  
  L1 <- L1.meas(group = c(rep(1, nt), rep(0, nc)),
                # data is a hack, need L1 var of some sort
                data = data[, -c(1, 12)],
                breaks = breaks,
                weights = c(tweight, cweight))$L1
  
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