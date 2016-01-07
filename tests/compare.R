library(WhatIf)
  
##################### Data Generation ##################

run <- 10

data <- list()
for(i in 1:run){
  data[[i]] <- generate.data()
}
    
RAW <- rep(NA, run) 
NEW <- rep(NA, run)
OLD <- rep(NA, run)

result <- data.frame(RAW = RAW, NEW = NEW, OLD = OLD)

hist_RAW <- list()
hist_NEW <- list()
hist_OLD <- list()

#################### RAW  #############################
    

for(i in 1:run) {
    
  this_data <- data[[i]]
  tgroup <- this_data[which(this_data$treat == 1),][1:5,]
  cgroup <- this_data[which(this_data$treat == 0),]
  nt <- 5
  nc <- nrow(cgroup)
  
  tweight <- rep(1, nt); cweight <- rep(nt/nc, nc)
  dat <- data.frame(X = c(tgroup$age, cgroup$age),
                    Weight = c(tweight, cweight),
                    Group = c(rep("Treat", nt), rep("Control", nc)))
    
  hist_RAW[i] <- ggplot(dat, aes(x = X, weight = Weight, fill = Group)) +
    geom_histogram(position = "dodge", bins = 30)  +
    xlab("age") +
    ylab("Weight")
    
  result$RAW[i] <- L1.single(treat = c(tgroup$treat, cgroup$treat),
                             data = c(tgroup$age, cgroup$age),
                             weights = c(tweight, cweight),
                             breaks = hist(c(tgroup$age, cgroup$age), plot = F)$breaks)
    
  }

####################### NEW #################################



for(i in 1:run) {
  
  this_data <- data[[i]]
  tgroup <- this_data[which(this_data$treat == 1),][1:5,]
  cgroup <- this_data[which(this_data$treat == 0),]
  nt <- nrow(tgroup)
  nc <- nrow(cgroup)
  
  A <- t(matrix(cgroup$age))
  # A <- t(matrix(cgroup$age[my.result$in.hull]))
  b <- mean(tgroup$age)
  
  c <- generate.c(tgroup$age, cgroup$age, sample_size = 3)
  c_lowerDim <- c[1:(length(c) - 2)]
  
  cweight <- nt *rowMeans(walkr(A = A, b = b, points = 10, c = c_lowerDim))
  
  tweight <- rep(1, nrow(tgroup))
  
  dat <- data.frame(X = c(tgroup$age, cgroup$age),
                    Weight = c(tweight, cweight),
                    Group = c(rep("Treat", nt), rep("Control", nc)))
  
  hist_NEW[i] <- ggplot(dat, aes(x = X, weight = Weight, fill = Group)) +
                      geom_histogram(position = "dodge", bins = 30)  +
                        xlab("age") +
                        ylab("Weight")
  
  
  result$NEW[i] <- L1.single(treat = c(tgroup$treat, cgroup$treat), 
                            data = c(tgroup$age, cgroup$age), 
                            breaks = hist(c(tgroup$age, cgroup$age), plot = FALSE)$breaks, 
                            weights = c(tweight, cweight))
}



############################## OLD  #############################################

remove.packages("walkr")
install_github("andyyao95/walkr", ref = "1b3ee93")
library(walkr)

for(i in 1:run) {
  
  this_data <- data[[i]]
  tgroup <- this_data[which(this_data$treat == 1),][1:5,]
  cgroup <- this_data[which(this_data$treat == 0),]
  nt <- nrow(tgroup)
  nc <- nrow(cgroup)
  A <- t(matrix(cgroup$age))
  b <- mean(tgroup$age)
  
  cweight <- nrow(tgroup)*rowMeans(walkr(A = A, b = b, points = 10))
  tweight <- rep(1, nrow(tgroup))
  
  dat <- data.frame(X = c(tgroup$age, cgroup$age),
                    Weight = c(tweight, cweight),
                    Group = c(rep("Treat", nt), rep("Control", nc)))
  
  hist_OLD[i] <- ggplot(dat, aes(x = X, weight = Weight, fill = Group)) +
    geom_histogram(position = "dodge", bins = 30)  +
    xlab("age") +
    ylab("Weight")
  
  
  result$OLD[i] <- L1.single(treat = c(tgroup$treat, cgroup$treat), 
                             data = c(tgroup$age, cgroup$age), 
                             breaks = hist(c(tgroup$age, cgroup$age), plot = FALSE)$breaks, 
                             weights = c(tweight, cweight))
}




data <- generate.data()
tgroup <- data[which(data$treat == 1),][1:5]; 
cgroup <- data[which(data$treat == 0),]; 
library(WhatIf)
my.result <- whatif(data = matrix(tgroup[,2], ncol = 1),
                    cfact = matrix(cgroup[,2], ncol = 1))
A <- t(matrix(cgroup$age))
# A <- t(matrix(cgroup$age[my.result$in.hull]))
b <- mean(tgroup$age)

c <- generate.c(tgroup$age, cgroup$age, sample_size = 3)
c_lowerDim <- c[1:(length(c)-2)]

cweight <- nrow(tgroup)*rowMeans(walkr(A = A, b = b, points = 10, c = c_lowerDim))

tweight <- rep(1, nrow(tgroup))

dat <- data.frame(X = c(tgroup[["age"]], cgroup[["age"]]),
                  Weight = c(tweight, cweight),
                  Group = c(rep("Treat", nrow(tgroup)),
                            rep("Control", nrow(cgroup))))

histplot <- ggplot(dat, aes(x = X, weight = Weight, fill = Group)) + 
  geom_histogram(position = "dodge", bins = 30)  +
  xlab("age") +
  ylab("Weight")


L1.single(treat = c(tgroup$treat, cgroup$treat), 
          data = c(tgroup$age, cgroup$age), 
          breaks = hist(c(tgroup$age, cgroup$age), 
                        plot = FALSE)$breaks, weights = c(tweight, cweight))


generate.c <- function (t_age, c_age, sample_size = 50, run = 100) {
  
  probs <- rep(0, length(c_age))
  library(WhatIf)
  
  for(i in 1:run) {
    sample_t <- sample(t_age, size = sample_size, replace = F)
    my.result <- whatif(data = matrix(sample_t, ncol = 1),
                        cfact = matrix(c_age, ncol = 1))$in.hull
    probs[which(my.result==T)] <- probs[which(my.result==T)] + 100
    
  }
  
  #normalize, scale up and return
  return (probs*length(t_age)/sum(probs))
}

