library("lpSolve")
# Parameters 
n = 8
S = 1000
M = 10.93
nu = c(1.03, 0.84, 1.15, 2.01, 1.28, 2.40, 1.22, 0.87)
W = diag(n) + rbind(rep(0,n) ,cbind( - 1 *diag(n-1), rep(0, n-1)))
omega = matrix(0, n, S)
for (i in 1:n) {
  omega[i,] =  rexp(S, rate = 1/nu[i])
  
}
############################################################################
############################################################################
# Expected Value problem: EV

model1 = list()
model1$modelsense <- "min"
model1$obj = cbind(rep(0,n), rep(1,n))
model1$A = rbind(c(rep(1,n), rep(0,n)), cbind(diag(n), W))
model1$rhs = c(M, nu)
model1$sense = c("<=", rep(">=", n))

result1 <- lp(model1$modelsense, model1$obj, model1$A, model1$sense, model1$rhs)


print('Solution:')
print(result1$objval)
# 0
print(result1$solution)
# 1.03 0.84 1.15 2.01 1.28 2.40 1.22 0.87 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
####################################################################################
#####################################################################################
# Expected result problem EVV 
x_bar = c(1.03, 0.84, 1.15, 2.01, 1.28, 2.40, 1.22, 0.87)
x_bar_S = rep(x_bar, S)
WS = diag(n*S) + rbind(rep(0,n*S) ,cbind( - 1 *diag(n*S-1), rep(0, n*S-1)))

model2 = list()
model2$modelsense <- "min"
model2$obj = rep(1/S, n*S)
model2$A = WS
model2$rhs = c(omega) - x_bar_S
model2$sense = rep(">=", n*S)

result2 <- lp(model2$modelsense, model2$obj, model2$A, model2$sense, model2$rhs)


print('Solution:')
print(result2$objval)
print(result2$solution)
#######################################################################################33
##########################################################################################
#Recourse model 
x_S = diag(n)
for (i in 1: (S-1)) {
  x_S = rbind(x_S, diag(n))
}


model3 = list()
model3$modelsense <- "min"
model3$obj = cbind(rep(0,n), rep(1/S, n*S))
model3$A =rbind(c(rep(1,n), rep(0,n*S)), cbind(x_S, WS))
model3$rhs = c(M,c(omega)) 
model3$sense = c ("<=", rep(">=", n*S))

result3 <- lp(model3$modelsense, model3$obj, model3$A, model3$sense, model3$rhs)


print('Solution:')
print(result3$objval)
print(result3$solution[1:n])
################################################################################################3
##########################################################################################3333
# Wait and see model
model4 = list()
model4$modelsense <- "min"
model4$obj = cbind(rep(0,n), rep(1/S, n*S))
model4$A =rbind(c(rep(1,n), rep(0,n*S)), cbind(diag(n*S), WS))
model4$rhs = c(omega) 
model4$sense = c ("<=", rep(">=", n*S))

result4 <- lp(model4$modelsense, model4$obj, model4$A, model4$sense, model4$rhs)


print('Solution:')
print(result4$objval)
print(result4$solution[1:n])
