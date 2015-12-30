dati1 <- generateData(DW)

####################### CEM: AGE ONLY ##############################

require(ggplot2)
require(reshape2)

# get all the names
name <- names(dati1)

# take out age
name <- name[-which(name == "age")]

# match only the age
cem.mat <- cem("treat", dati1, drop = name)

# indices of matched treated
cem.tr <- which(cem.mat$groups == "1" & cem.mat$matched == TRUE)

# indices of matched control
cem.ct <- which(cem.mat$groups == "0" & cem.mat$matched == TRUE)

cem.idx <- unique(c(cem.tr, cem.ct))
L1.meas(dati1$treat[cem.idx], dati1[cem.idx, -c(1, 12)], breaks = mybr)$L1

age.tr <- dati1$age[cem.tr]
age.ct <- dati1$age[cem.ct]

# join ages from matched treatment group with those from matched
# control group in a data frame
df <- as.data.frame(t(rbind(age.ct, age.tr[seq(age.ct)])))

# name the data frame
names(df) <- c("control", "treat")

# a joined histogram
ggplot(melt(df, na.rm = TRUE, id.vars = NULL),
       aes(value, fill = variable)) +
  geom_histogram(position = "dodge", bins = 30) +
  xlab("Age - CEM") +
  ylab("Count")

mean(dati1$outcome[cem.tr]) - mean(dati1$outcome[cem.ct])


################ MANUAL: WALKR AGE ONLY ###########################
library(walkr)

treat <- dati1[which(dati1$treat == 1), ]
control <- dati1[which(dati1$treat == 0), ]

nt <- nrow(treat); nc <- nrow(control)

# an even weight
wt <- rep(1/nt, nt)

# match on age
weight <- walkr(rbind(matrix(control$age, ncol = nc), rep(1, nc)),
                c(t(wt)%*%treat$age, 1),
                points = 100,
                thin = 2)

# take the mean and scale up
age.ct <- control$age; weight.ct <- rowMeans(weight) * nt

age.tr <- treat$age; weight.tr <- rep(1, nt)

dat <- data.frame(x = c(age.tr, age.ct),
                  freq = c(weight.tr, weight.ct),
                  grp = c(rep("Treat", length(age.tr)),
                          rep("Control", length(age.ct))))

ggplot(dat,aes(x = x,weight = freq,fill = grp)) + 
  geom_histogram(position = "dodge", bins = 30)  +
  xlab("Age - Walkr 1") +
  ylab("Count")

weighted.mean(treat$outcome, wt) - weighted.mean(control$outcome, weight.ct)


################# MANUAL: WALKR AGE + MSQ AGE #######################
library(walkr)

treat <- dati1[which(dati1$treat == 1), ]
control <- dati1[which(dati1$treat == 0), ]

nt <- nrow(treat)
nc <- nrow(control)

wt <- rep(1/nt, nt)

A <- rbind(matrix(control$age, ncol = nc),
           matrix(control$age^2, ncol = nc),
           rep(1, nc))

b <- c(t(wt)%*%treat$age,
       t(wt)%*%(treat$age^2),
       1)

# match on age, and MSQ of age
weight <- walkr(A, b, points = 20)

# take the mean and scale up
age.ct <- control$age; weight.ct <- rowMeans(weight) * nt

age.tr <- treat$age; weight.tr <- rep(1, nt)


dat <- data.frame(x = c(age.tr, age.ct),
                  freq = c(weight.tr, weight.ct),
                  grp = c(rep("Treat", length(age.tr)), 
                          rep("Control", length(age.ct))))

ggplot(dat, aes(x = x, weight = freq, fill = grp)) + 
  geom_histogram(position = "dodge", bins = 30) + 
  xlab("Age - Walkr 2") +
  ylab("Count")

t(treat$outcome)%*%wt - t(control$outcome)%*%weight.ct/nt
