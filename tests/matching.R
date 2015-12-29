require(cem)
require(MatchIt)
require(Matching)
data(DW)
tsubjects <- which(DW$treated == 1)
csubjects <- which(DW$treated == 0)

# remove the treated variable (will be simulated later)
dati <- DW[, -c(1, 9)]

# Y - var
outcome <- DW$re78

# num of treated
nt <- length(tsubjects)

# num of control
nc <- length(csubjects)

# num of total
n <- nt + nc

# treated: a vector of boolean
treated <- logical(n)
treated[tsubjects] <- TRUE

# force dati (with treated removed) to be numeric for regression purpose
dati.num <- dati
for(i in 1:dim(dati)[2]) dati.num[, i] <- as.numeric(dati[, i])

propensity  <- glm(treated~ I(age^2) + I(education^2) + black +
                     hispanic + married + nodegree + I(re74^2) + I(re75^2) +
                     u74 + u75, family = binomial(link = "logit"), data = dati)
M <- cbind(rep(1, n),
           propensity$linear.pred,
           I(log(dati$age)^2),
           I(log(dati$education)^2),
           I(log(dati$re74+0.01)^2),
           I(log(dati$re75+0.01)^2))

# misspecified weights for: Intercept, linear.pred, age, educ, I(re74^2), and I(re75^2)
propensity.coeffs <- as.matrix(c(1.00, 0.5, 0.01, -0.3, -0.01, 0.01))

# mu: log odds
mu = M %*% propensity.coeffs

# Tr.pred: actual probability, through the logit transformation
Tr.pred <- exp(mu)/(1+exp(mu))
TreatmentEffect <- 1000
TreatmentReal <- matrix(nrow = n, ncol = 1)

for(i in 1:n){
  TreatmentReal[i] = sample(0:1, 1, prob = c(1 - Tr.pred[i], Tr.pred[i]))
}

# outcome assigns a treatment effect of 1000; error ~ N(0, 100)  
outcome <- I(TreatmentEffect*TreatmentReal) +
  .1*exp(.7*log(dati$re74+0.01) + .7*log(dati$re75+0.01)) + rnorm(n, 0, 10)

treated <- TreatmentReal
tsubjects <- which(TreatmentReal==1)
csubjects <- which(TreatmentReal==0)
nt <- length(tsubjects)
nc <- length(csubjects)

# use the simulated treat vector (as opposed to the original treatment indicator in LeLonde!)
dati1 <- data.frame(treat = treated, dati)

############################## CEM AGE ONLY ####################################

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

age.tr <- dati1$age[cem.tr]
age.ct <- dati1$age[cem.ct]

# join ages from matched treatment group with those from matched
# control group in a data frame
df <- as.data.frame(t(rbind(age.ct, age.tr[seq(age.ct)])))

# name the data frame
names(df) <- c("age.ct", "age.tr")

# a joined histogram
ggplot(melt(df, na.rm = TRUE, id.vars = NULL),
       aes(value, fill = variable)) +
  geom_histogram(position = "dodge", bins = 30)



########################### WALKR AGE ONLY ################################
library(walkr)

dati1 <- cbind(dati1, outcome)

treat <- dati1[which(dati1$treat == 1), ]
control <- dati1[which(dati1$treat == 0), ]

nt <- nrow(treat)
nc <- nrow(control)

wt <- rep(1/nt, nt)

# match on age
weight <- walkr(rbind(matrix(control$age, ncol = nc), rep(1, nc)),
                c(t(wt)%*%treat$age, 1),
                points = 20)

# take the mean and scale up
weight.ct <- rowMeans(weight) * nt
age.ct <- control$age

age.tr <- treat$age
weight.tr <- rep(1, nt)

# weighted.hist <- function (a, a_weight, b, b_weight, bins = 50) {
#   
#   b <- b + 0.5
#   
#   # building fill vector
#   my_col <- c(rep("blue", length(b)), rep("red", length(a)))
#   # plotting
#   qplot(x = c(a, b), weight = c(a_weight, b_weight), fill = my_col) +
#     geom_histogram(position = "dodge", bins = bins)
# }
# 
# weighted.hist(age.tr, weight.tr, age.ct, weight.ct)


dat <- data.frame(x = c(age.tr, age.ct),
                  freq = c(weight.tr, weight.ct),
                  grp = c(rep("age.tr", length(age.tr)),
                          rep("age.ct", length(age.ct))))

ggplot(dat,aes(x = x,weight = freq,fill = grp)) + 
  geom_histogram(position = "dodge", bins = 30)

TE <- t(treat$outcome)%*%wt - t(control$outcome)%*%weight.ct/nt


########################### WALKR AGE + MSQ AGE ################################
library(walkr)

dati1 <- cbind(dati1, outcome)

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
weight <- walkr(A, b, points = 10)

# take the mean and scale up
age.ct <- control$age
weight.ct <- rowMeans(weight) * nt

age.tr <- treat$age
weight.tr <- rep(1, nt)

# weighted.hist <- function (a, a_weight, b, b_weight, bins = 50) {
#   
#   b <- b + 0.5
#   
#   # building fill vector
#   my_col <- c(rep("blue", length(b)), rep("red", length(a)))
#   # plotting
#   qplot(x = c(a, b), weight = c(a_weight, b_weight), fill = my_col) +
#     geom_histogram(position = "dodge", bins = bins)
# }
# 
# weighted.hist(age.tr, weight.tr, age.ct, weight.ct)


dat <- data.frame(x = c(age.tr, age.ct),
                  freq = c(weight.tr, weight.ct),
                  grp = c(rep("age.tr", length(age.tr)), 
                          rep("age.ct", length(age.ct))))

ggplot(dat, aes(x = x, weight = freq, fill = grp)) + 
  geom_histogram(position = "dodge", bins = 30)

TE <- t(treat$outcome)%*%wt - t(control$outcome)%*%weight.ct/nt
