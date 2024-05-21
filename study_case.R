library("lpSolve")
set.seed(48643)
# Parameters 
n = 8
S = 1000
M = 10.93
nu = c(1.03, 0.84, 1.15, 2.01, 1.28, 2.40, 1.22, 0.87)
W = diag(n) + rbind(rep(0,n) ,cbind( - 1 *diag(n-1), rep(0, n-1)))

############################################################################
############################################################################
# Expected Value problem: EV################################################

modelEV = list()
modelEV$modelsense <- "min"
modelEV$obj = c(rep(0,n), rep(1,n))
modelEV$A = rbind(c(rep(1,n), rep(0,n)), cbind(diag(n), W))
modelEV$rhs = c(M, nu)
modelEV$sense = c("<=", rep(">=", n))

resultEV <- lp(modelEV$modelsense, modelEV$obj, modelEV$A, modelEV$sense, modelEV$rhs)

print('Solution:')
print(resultEV$objval)
# 0
print(resultEV$solution)

####################################################################################
#####################################################################################
# Expected result EVV  ################################################

# y_i(\omega) values are calculated by taking the maximum of two lower bounds on each y_i(\omega) (d_i formulas)
# Faster approach for calculating the expected value

rep = 100
S = 5000

x_EV = resultEV$solution[1:n]
EEVs = c()
for (it in 1:rep) { # Multiple Replications Procedure
   omegaEEV = matrix(rexp(S*n, rate = 1/nu),n)
   objValsEEV <- c()
   for (s in 1:S) {
     d = max(omegaEEV[1,s] - x_EV[1], 0)
     for (i in 2:n) {
       d = c(d, max(omegaEEV[i,s] - x_EV[i] + d[i-1], 0) )
     }
     objValsEEV = c(objValsEEV, sum(d))
   }
   EEVs = c(EEVs, mean(objValsEEV))
 }

cat("Mean EEV Value: ", mean(EEVs) ,'\n')
cat("CI on EEV: ", mean(EEVs) + sd(EEVs)*qt(0.975,rep-1)*c(-1,1) ,'\n')
cat("Relative error: ", sd(EEVs)/mean(EEVs)*100, "%" ,'\n')
cat('EV: ', resultEV$objval,'\n')

##########################################################################################
#Recourse model ################################################
S = 2000
x_S = diag(n)
for (i in 1:(S-1)) {
  x_S = rbind(x_S, diag(n))
}
WS = diag(n*S) + rbind(rep(0,n*S), cbind( - 1 *diag(n*S-1), rep(0, n*S-1)))
for (i in 1: S-1) {
  WS[i*n + 1, i*n] = 0 
}

omegaRM = matrix(rexp(S*n, rate = 1/nu),n)

model3 = list()
model3$modelsense = "min"
model3$obj = c(rep(0,n), rep(1/S, n*S))
model3$A = rbind(c(rep(1,n), rep(0,n*S)), cbind(x_S, WS))
model3$rhs = c(M,c(omegaRM)) 
model3$sense = c("<=", rep(">=", n*S))

result3 <- lp(model3$modelsense, model3$obj, model3$A, model3$sense, model3$rhs)

cat('TS Solution:',result3$objval,'\n')
cat('x values:', result3$solution[1:n],'\n')
cat('Sum of x values:', sum(result3$solution[1:n]),'\n')

# Wait and see model ################################################
S = 5000
omegaWS = matrix(rexp(S*n, rate = 1/nu),n)

objValsWS <- c()
solsWS <- c()
for (s in 1:S) {
  modelWS = list()
  modelWS$modelsense <- "min"
  modelWS$obj = c(rep(0,n), rep(1,n))
  modelWS$A = rbind(c(rep(1,n), rep(0,n)), cbind(diag(n), W))
  modelWS$rhs = c(M, omegaWS[,s])
  modelWS$sense = c("<=", rep(">=", n))
  resultWS <- lp(modelWS$modelsense, modelWS$obj, modelWS$A, modelWS$sense, modelWS$rhs)
  objValsWS <- c(objValsWS, resultWS$objval)
  solsWS <- c(resultWS$solution)
}


print('Solution:')
print(mean(objValsWS))

### question b: Sensitivity #############
objValsB <- c()
solsB <- c()
for (M_b in seq(9.5, 12.5, 0.5)) {
  model3$rhs = c(M_b,c(omegaRM))
  newResult <- lp(model3$modelsense, model3$obj, model3$A, model3$sense, model3$rhs)
  objValsB <- c(objValsB, newResult$objval)
  solsB <- rbind(solsB, newResult$solution[1:n])
}
cat('M: ', seq(9.5, 12.5, 0.5), '\n')
cat('Obj: ', objValsB, '\n')
cat('x: \n', solsB, '\n')

### question c ################################################
S = 2000
count = 0
dum = result3$solution[(n+1):((n+1)*S)]
for (i in seq(8, (n*S), by = 8)) {
  if (dum[i] > 3) {
    count = count + 1
  } 
}
cat('count', count, 'Prob', count/S,'\n')

################################################################################################
#########################################################################################
# question d################################################
# Proportional Model ####
x_prop = c(1.04, 0.85, 1.16, 2.03, 1.29, 2.42, 1.23, 0.88)
omegaDprop = matrix(rexp(S*n, rate = 1/nu),n)
objValsDprop = c()
for (s in 1:S) {
  d = max(omegaDprop[1,s] - x_prop[1], 0)
  for (i in 2:n) {
    d = c(d, max(omegaDprop[i,s] - x_prop[i] + d[i-1], 0) )
  }
  objValsDprop = c(objValsDprop, sum(d))
}
cat("Average objective value for proportionality approach:", mean(objValsDprop),'\n')

# Expected value of TS Solution ####

x_TS = result3$solution[1:n]
omegaDts = matrix(rexp(S*n, rate = 1/nu),n)
objValsDts = c()
for (s in 1:S) {
  d = max(omegaDts[1,s] - x_TS[1], 0)
  for (i in 2:n) {
    d = c(d, max(omegaDts[i,s] - x_TS[i] + d[i-1], 0) )
  }
  objValsDts = c(objValsDts, sum(d))
}
cat("Average objective value for TS approach:", mean(objValsDts),'\n')