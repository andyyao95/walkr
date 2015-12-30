generate.data <- function(df){
  
  fdasfaf <- DW$treated
  
  tsubjects <- which(DW$treated == 1)
  csubjects <- which(DW$treated == 0)
  
  # remove the treated variable (will be simulated later)
  dati <- DW[, -c(1, 9)]
  
  # Y - var
  outcome <- DW$re78
  
  # num of treated and control
  nt <- length(tsubjects); nc <- length(csubjects)
  
  # num of total
  n <- nt + nc
  
  # treated: a vector of boolean
  treated <- logical(n); treated[tsubjects] <- TRUE
  
  # convert dati (with treated removed) to numeric for regression purpose
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
  Tr.pred <- exp(mu)/(1 + exp(mu))
  TreatmentEffect <- 1000
  
  # sample from Tr.pred
  TreatmentReal <- matrix(sample(0:1,
                                 n,
                                 prob = c(1 - Tr.pred[i], Tr.pred[i]),
                                 replace = TRUE),
                          ncol = 1)
  
  # outcome assigns a treatment effect of 1000; error ~ N(0, 100)  
  outcome <- I(TreatmentEffect*TreatmentReal) +
    0.1*exp(0.7*log(dati$re74 + 0.01) + 0.7*log(dati$re75 + 0.01)) + rnorm(n, 0, 10)
  
  treated <- TreatmentReal
  tsubjects <- which(TreatmentReal == 1)
  csubjects <- which(TreatmentReal == 0)
  nt <- length(tsubjects); nc <- length(csubjects)
  
  # use the simulated treat vector
  # (as opposed to the original treatment indicator in LeLonde!)
  dati1 <- data.frame(treat = treated, dati)
  return(cbind(dati1, outcome))
}