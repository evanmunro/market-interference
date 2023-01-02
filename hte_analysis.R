#me is micro-enterprises
library(readstata13)
library(grf)
library(pracma)
raw = read.dta13("investments_data.dta")
arpmap = data.frame(names=colnames(raw), labels = attr(raw, "var.labels"))

data <- raw[raw$wave %in% c(2, 3, 4),]


data$ha[is.na(data$ha)] <- 0
data$crop_sales[is.na(data$crop_sales)] <- 0 
data$ani_sales_a[is.na(data$ani_sales_a)] <- 0 
data$homeprod[is.na(data$homeprod)] <- 0 
data$home_prod_w3 <- data$crop_sales/12 + data$ani_sales_a/6 + data$homeprod #+ data["ani_prod_sales_a"]
#code 0 if missing production
#data["home_prod_w3"][is.na(data["home_prod_w3"])] = 0 
#data <- data[data["home_prod_w3"]!=0, ]
#drop data with missing animal ownership 
#outcome = data$ha
missing = is.na(data$t2_c1_op) + is.na(data$ha) + is.na(data$hhsize_ae2) + is.na(data$no497)
data = data[!missing, ]
W = data$t2_c1_op==1
# per person agricultural income
Y = data$home_prod_w3/data$hhsize_ae2
Z = data$ha
Y[Y > quantile(Y, 0.99)] = quantile(Y, 0.99)
Z[Z > quantile(Z, 0.99)] = quantile(Z, 0.99)

#Sanity check of the ATEs
mean(Y[W==1]) - mean(Y[W==0])

mean(Z[W==1]) - mean(Z[W==0])

#up2_mcwagepm is a municipal wage measure
#ethnicity_hh is a indigenous language metric
#educ1_hh is head of household education metric
#edu1_sp is spouse education metric
#dage are demographic binary variables 
#no497 is if not farm house 
#big497 is if large farm in 97 
#ha97 is how many metrics in 97 
#t2_c1_op is ITT 

#other options: bathroom97, dage55p97, hhsize97, land97, landown97, hairrigation97, haseasonality97, ani97, small497
#vani97 is value 
#landlessmall497 
#tortilla97, procampo97, probecat97, solidaridad97, INI97, empleotemp97, 
#nbani97 is bundle in cows equivalent animals 
#tortilla is tortilla support , proxcampo is agricultural subsidiies , empleo temporal support, probecat support,
#INI is indigenous porverty support
X97 = data[, c("bathroom97", "dage55p_97", "landown97", "vani97", "dage0_7_97", "dage8_17_97", "dage18_54_97", "no497", "big497", "ha97", 
               "nbani197", "electricity97",
               "homeown97", "dirtfloor97", "tortilla97", "procampo97", "probecat97", "solidaridad97", "INI97", "empleotemp97")]

X = data[, c("dage0_7_97", "dage8_17_97", "dage18_54_97", "no497", "big497", "ha97", "nbani197", "electricity97", 
                  "homeown97", "dirtfloor97", "educ1_hh", "educ1_sp", "female_hh", 
                  "ethnicity_hh", "up2_mcwagepm")]

set.seed(12)

#Task 1: Compute optimal rule on a single training split and save the predicted CATEs
#for plotting in Julia

train <- sample(1:n, n / 2)

hattauY = causal_forest(X[train,], Y[train], W[train], num.trees=200)$predictions
hattauZ = causal_forest(X[train,], Z[train], W[train], num.trees=200)$predictions

write.csv(data.frame("Z" = hattauZ, "Y"= hattauY), file="hattau.csv")

forest.W <- regression_forest(X[train,], W[train], num.trees=200)
W.hat <- predict(forest.W)$predictions
eval.c <- function(c) {
  compute_gain(Z[train], W[train], hattauY > c*hattauZ, mean(W))
}
cstar = fzero(eval.c, 0.0)$x
print(cstar)


#TASK 2: Bootstrap the lift computation
tauhatin <- function(Y, X, W, idx) {
  eval.forest <- causal_forest(X[idx, ], Y[idx], W[idx], num.trees=200)
  list(predict(eval.forest, X[idx,])[,1], predict(eval.forest, X[-idx,])[,1])
}

pscorein <- function(W, X, idx){
  forest.W <- regression_forest(X[idx,], W[idx], num.trees=200)
  W.hat.train <- predict(forest.W)$predictions
  W.hat.test <- predict(forest.W, X[-idx,])[,1]
  list(W.hat.train, W.hat.test)
}

REPS = 2

set.seed(12)
results = replicate(n = 100, compute_lift(Y, X, W, Z))
print("LIFT:")
print(mean(results[1,], na.rm=TRUE))
print(sd(results[1,], na.rm=TRUE))
print("STABILITY:")
print(mean(results[2,], na.rm=TRUE))
print(sd(results[2,], na.rm=TRUE))

print("NAIVE LIFT:")
print(mean(results[4,], na.rm=TRUE))
print(sd(results[4,], na.rm=TRUE))
print("NAIVE STABILITY:")
print(mean(results[5,], na.rm=TRUE))
print(sd(results[5,], na.rm=TRUE))
print("Cutoff:")
print(mean(results[3,], na.rm=TRUE))
print(sd(results[3,], na.rm=TRUE))

compute_gain <- function(Y, W, S, q) {
  scoreLeft = (W*Y/q - (1 - W)*Y/(1 - q))*S
  scoreRight = (W*Y/q - (1 - W)*Y/(1 - q))
  mean(scoreLeft) - mean(scoreRight)*q
}


compute_lift <- function(Y, X, W, Z, q) {
  n = length(W)
  rsamp <- sample(1:n, n, replace=TRUE)
  Y = Y[rsamp]
  X = X[rsamp,]
  W = W[rsamp]
  Z = Z[rsamp]
  train <- sample(1:n, n / 2)
  yhat <- tauhatin(Y, X, W, train)
  zhat <- tauhatin(Z, X, W, train)
  W.hat <- pscorein(W, X, train)
  W.hat.train <- W.hat[[1]]
  W.hat.test <- W.hat[[2]]
  yhat.train <- yhat[[1]]
  yhat.test <- yhat[[2]]
  zhat.train <- zhat[[1]]
  zhat.test <- zhat[[2]]
  
  compute.c <- function(c) {
    compute_gain(Z[train], W[train], yhat.train > zhat.train*c, q)
  }
  cstar = fzero(compute.c, 0.0)$x
  lift = compute_gain(Y[-train], W[-train], yhat.test > zhat.test*cstar, q)
  stability =  compute_gain(Z[-train], W[-train], yhat.test > zhat.test*cstar, q)
  naive_lift = compute_gain(Y[-train], W[-train], yhat.test > 0, q)
  naive_stability =  compute_gain(Z[-train], W[-train], yhat.test > 0, q)
  c(lift, stability, cstar, naive_lift - lift, abs(naive_stability) - abs(stability))
}
  



