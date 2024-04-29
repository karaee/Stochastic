################################################
###########     Farmer's problem      ##########
################################################
# In this R-script you can find a sample code for solving
# the farmer's problem introduced in the textbook by Louveaux and Birge.
# We will solve the deterministic problem and the stochastic problem
# The latter using the large-scale deterministic equivalent problem and
# using a basic L-shaped method.
# To implement these algorithms we need an LP solver. Here we use lpSolve (a package from R).
# It is not as good as Cplex or GUROBI, but the package is easy to install and does not require a license

# NOTE: This code is work-in-progress. If you find any errors, or have suggestions for improvement,
#.      Please let me know.


# Install the package lpSolve
################################################
############    lpSolve Example      ############
################################################
# Below we show how to use lpSolve to solve a simple LP problem:
# max x1 + x2 + 2*x3 : x1 + 2x2 + 3x3 <= 4 and x1 + x2 >= 1 

# Load the lpSolve package
library('lpSolve')

# create a "model" (=list containing all model characteristics)
model <- list()


# Specify all coefficients/characteristic of the model
#### max cx: Ax <= b
#### It is possible to specify whether you are maximizing/minimizing
#### whether inequalities are >=, <=
#### and types of the variables. Continuous variables are always ">= 0".

model$A          <- matrix(c(1,2,3,1,1,0), nrow=2, ncol=3, byrow=T)
model$obj        <- c(1,1,2)
model$modelsense <- "max"
model$rhs        <- c(4,1)
model$sense      <- c("<=", '>=')


result <- lp(model$modelsense, model$obj, model$A, model$sense, model$rhs)



# Print the solution:
print('Solution:')
print(result$objval)
print(result$solution)


#####################################################
########## Deterministic Farmer's Problem  #########
#####################################################
### Here we solve the deterministic farmer's problem
### Notice that we already use SP notation, where x
### are first-stage variables, (y,w) second-stage variables
### with unit cost vectors c and q.
### A, b correspond to first-stage constraints, W and h and the
### technology matrix T (=Tech) to second-stage constraints

# For my convenience, I make all second-stage inequalities '>='
# Check for yourself that all matrices of coefficients etc. are correct

library('lpSolve')

A    <- matrix(c(1,1,1),nrow = 1)
Tech <- rbind(diag(c(2.5,3,20)),matrix(rep(0,3),nrow=1))
W    <- matrix(c(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,-1,-1),c(0,0,0,0,-1,0)),nrow = 4,byrow=T)
b    <- 500 
c    <- c(150,230,260)
q    <- c(238,-170,210,-150,-36,-10)
h    <- c(200,240,0,-6000)

FP <- list()


FP$A          <- cbind(rbind(A,Tech),rbind(matrix(rep(0,dim(A)[1]*dim(W)[2]),nrow=dim(A)[1]),W))
FP$obj        <- c(c,q)
FP$modelsense <- "min"
FP$rhs        <- c(b,h)
FP$sense      <- c('<=','>=','>=','>=','>=')
FP$vtype      <- 'C'



result <- lp(FP$modelsense, FP$obj, FP$A, FP$sense, FP$rhs)


print('Solution:')
print(result$objval)
print(result$solution)
#print(result)


###################################################################
############# Stochastic Farmer's Problem #########################
###################################################################
# Next we solve a stochastic farmer's problem where we have 3 possible
# values for the yield. We solve this problem using the large-scale
# deterministic equivalent problem. We set up the code in such a way that
# it can also be used to solve the LSDE of an SAA (sample average approximation).

library('lpSolve')

# This function creates zero-matrices of size m x n
zeros <- function(m,n){
	Z <- matrix(rep(0,n*m),nrow=m)
	return(Z)
}



# Here, there are three options for specifying the distribution of the "yield" xi (or rather %dev of
# avg yield).
# The first is with three scenarios xi = 0.8, 1, 1.2; all with probability 1/3.
# The second is by drawing a sample of S scenarios from the distribution of xi (As an example we 
# assume that xi is uniformly distributed on [0.8,1.2])
# The third is the same as the second, but uses Latin Hypercube Sampling instead of Monte Carlo sampling

# S  = number of scenarios
# p  = "vector of probabilities"
# xi = "vector of realizations"

# Option 1
S <- 3
p <- rep(1/S,S)
xi <- c(0.8,1,1.2)


# Option 2; My computer can handle S = 1000;  S=10000 is already (too) difficult
# Uncomment to use this option:
#S <- 1000
#p <- rep(1/S,S)
#xi <- runif(S, min = 0.8, max = 1.2)    # If you use sampling; try running your algorithm/code twice
                                        # and check if the solution values, and objective change
                                        # and by how much.

# Option 3; My computer can handle S = 1000;  S=10000 is already (too) difficult
# Uncomment to use this option:
#library('DiceDesign')
#S <- 1000
#p <- rep(1/S,S)
#xi <- lhsDesign(S,1,randomized=TRUE)$design*0.4+0.8    #lhsDesign creates "random numbers" 
                                                           # (between 0 and 1). You need to apply the inverse
                                                           # transform (F^{-1}) to these numbers. In our
                                                           # example this means *0.4 and +0.8 to get a 
                                                           # uniform distribution on [0.8, 1.2]     




# In this example there is only uncertainty in the technology matrix T. In fact, using T from the 
# deterministic setting and a realization xi[i]. The technolgy matrix in scenario i equals xi[i]*T
# For other problems/case study, the dependence of h,T,q on xi/omega may be "more difficult"

# This is the "average/nominal" technology matrix
Tech <- rbind(diag(c(2.5,3,20)),matrix(rep(0,3),nrow=1))

A    <- matrix(c(1,1,1),nrow = 1)
W    <- matrix(c(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,-1,-1),c(0,0,0,0,-1,0)),nrow = 4,byrow=T)
b    <- 500 
c    <- c(150,230,260)
q    <- c(238,-170,210,-150,-36,-10)
h    <- c(200,240,0,-6000)

# Dimension of the problem: m = "number of constraints"; n = "number of variables; 1st and 2nd stage
m1 <- dim(A)[1]
m2 <- dim(W)[1]
n1 <- dim(A)[2]
n2 <- dim(W)[2]



# What remains is to create the matrix of coefficients for this large-scale problem
# It has the form (if there are 4 scenarios):
#. ( A   0  0  0  0 )
#. ( T1  W  0  0  0 )
#. ( T2  0  W  0  0 )
#. ( T3  0  0  W  0 )
#. ( T4  0  0  0  W )


# Verify that the following lines of code gives the above matrix; If you forgot what 
# the kronecker product of two matrices is: Google it! or enter ?kronecker
# Also notice that kronecker(xi,Tech) works here because T1 = xi[1]*Tech

Z  <- zeros(m1,n2*S)
LSDE.A <- cbind(A,Z)
LSDE.A <- rbind(LSDE.A,cbind(kronecker(xi,Tech),kronecker(diag(rep(1,S)),W)))


# Construct the other parameters/coefficients.
# Notice the use of the "rep" function, because for all scenarios, right-hand sides and objectives 
# are the same/similar

FP <- list()

FP$A          <- LSDE.A
FP$obj        <- c(c,p*rep(q,S))
FP$modelsense <- "min"
FP$rhs        <- c(b,rep(h,S))
FP$sense      <- c('<=',rep(c('>=','>=','>=','>='),S))



result <- lp(FP$modelsense, FP$obj, FP$A, FP$sense, FP$rhs)

print('Solution:')
print(result$objval)
print(result$solution[1:3])  #only prints first-stage solutions



#############################################################################
################## L-shaped algorithm #######################################
#############################################################################
# Below is a basic L-shaped algorithm for solving the same 
# Stochastic Farmer's problem

# Initialization (= the same as for LSDE)
library('lpSolve')


zeros <- function(m,n){
	Z <- matrix(rep(0,n*m),nrow=m)
	return(Z)
}


# Option 1
S <- 3
p <- rep(1/S,S)
xi <- c(0.8,1,1.2)


# Option 2; My computer can handle S = 10000 and even 100000 (It might take a while);
# This is already more than you typically need in practice
# Uncomment to use this option:
#S <- 1000
#p <- rep(1/S,S)
#xi <- runif(S, min = 0.8, max = 1.2)         

# Option 3; My computer can handle S = 10000 and even 100000 (It might take a while);
# This is already more than you typically need in practice
# Uncomment to use this option:
#library('DiceDesign')
#S <- 1000
#p <- rep(1/S,S)
#xi <- lhsDesign(S,1,randomized=TRUE)$design*0.4+0.8        







Tech <- rbind(diag(c(2.5,3,20)),matrix(rep(0,3),nrow=1))

A    <- matrix(c(1,1,1),nrow = 1)
W    <- matrix(c(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,-1,-1),c(0,0,0,0,-1,0)),nrow = 4,byrow=T)
b    <- 500 
c    <- c(150,230,260)
q    <- c(238,-170,210,-150,-36,-10)
h    <- c(200,240,0,-6000)

m1 <- dim(A)[1]
m2 <- dim(W)[1]
n1 <- dim(A)[2]
n2 <- dim(W)[2]

# Create master problem of the form min cx + theta^+ - theta^-: Ax <= b
# Remember: all continuous variables are non-negative!
# Later we will add optimality cuts of the form theta^+-theta^- >= a*x+b
# We set a lower bound on theta = theta^+ - thetat^-  -> upper bound on theta^-

Master <- list()

Master$A <- cbind(A,zeros(m1,2))
Master$A <- rbind(Master$A,c(rep(0,n1+1),1))
Master$obj        <- c(c,1,-1)
Master$modelsense <- "min"
Master$rhs        <- c(b,10**10) 
Master$sense      <- c('<=','<=')



# Solve the master problem

Master.result <- lp(Master$modelsense, Master$obj, Master$A, Master$sense, Master$rhs)

print('Solution:')
print(Master.result$objval)
print(Master.result$solution)   # Not surprising the first time we solve the master problem x = 0


# let current.x denote the current first-stage solution:
current.x <- Master.result$solution[1:n1]


# Next initalize the subproblem(s). We will use the same subproblem iteratively
# since only the rhs of these subproblems change. They are given by:
# v(h-T[i]x) = min qy: Wy >= h-T[i]x
# Since lpSolve does not provide dual solutions, we solve the dual of v:
# Dual problem: max lambda^T(h-T[i]x): W^T lambda <= q^T

Sub <- list()
Sub$A          <- t(W)
Sub$obj        <- h - xi[1]*(Tech %*% current.x)
Sub$modelsense <- "max"
Sub$rhs        <- q
Sub$sense      <- rep('<=',n2)


# In the following repeat loop, we iteratively solve the subproblems, add optimality cuts
# and solve the master problem until the optimality criterion is met.
# We keep track of the number of iterations. I've set the max equal to 1000. You can change that 
# of course.

it = 0
repeat{

	current.x <- Master.result$solution[1:n1]
	
	# At each iteration we need to determine Q(x) and a subgradient u of Q at x,
	# where x is the current first-stage solution.
	# The objective value of the subproblem equals v(h-T[i]x)
	# -> Q(x) = sum p[i]v(h-T[i]x)
	# The dual solution lambda[i,x] is a subgradient of v
	# Then, u := -sum p[i]*lambda[i,x]*T[i] is the subgradient of Q at x:
	# -> Q(x) >= Q(x.current) + u*(x-x.current)      <=>
	# Optimality cut: theta^+-theta^- - u*x >= Q(x.current) - u*x.current
	
	Q <- 0
	u <- 0
	
	# We have to calculate T[i]x several times, which is equal to xi[i]*Tech %*% x
	# Hence, calculate Tx = Tech %*% x here
	
	# Observe that u can for this example also be calculated using
	# u := (-sum p[i]*lambda[i,x]*xi[i])*Tech
	
	Tx <- Tech %*% current.x 

	for(i in 1:S){
		Sub$obj <- h - xi[i]*Tx
		Sub.result <- lp(Sub$modelsense, Sub$obj, Sub$A, Sub$sense, Sub$rhs)
		Q <- Q + p[i]*Sub.result$objval
    	u <- u - p[i]*xi[i]*Sub.result$solution 
	}
	
	# Calculate u:
	u <- u %*% Tech
	
	# Stopping criterion. Here I have selected "epsilon" = 10**-5. 	
	if(c %*% current.x+Q < Master.result$objval + 10**-5){
		break
	}
	
	# Add optimality cut and resolve the master problem
	Master$A     <- rbind(Master$A,c(-u,1,-1))
	Master$rhs   <- c(Master$rhs,Q-u %*% current.x)
	Master$sense <- c(Master$sense,'>=') 
	Master.result <- lp(Master$modelsense, Master$obj, Master$A, Master$sense, Master$rhs)

	print('Iteration:')
	print(it)
	print('Objective value current solution')
	print(c %*% current.x+Q)
	print('Solution:')
	print(Master.result$objval)
	print(Master.result$solution)
	flush.console()
	
	it <- it + 1
	if(it > 1000){
		break
    }
}	

# print the solution 
print('Solution:')
print(Master.result$objval)
print(Master.result$solution[1:3]) 





