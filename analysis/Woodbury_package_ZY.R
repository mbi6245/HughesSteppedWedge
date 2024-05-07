# Apply Woodbury formular to inverse V, using the WoodburyMatrix package

library("WoodburyMatrix")
library("sandwich")
library("clubSandwich")
library("lme4")


# Import data and fit lmer model, creat factors for clustering

load("ctdata.RData")
ctdata$fcluster = factor(ctdata$hdist)
ctdata$ftime = factor(ctdata$time)
ctdata$clustime = interaction(ctdata$fcluster,ctdata$ftime)
rslt = lmer(ct ~ ftime + treat + (1 | fcluster) + (1 | clustime), data=ctdata)

# Extract terms we need from lmer fitted object, rslt

beta = matrix(fixef(rslt),ncol=1)
np = dim(beta)[1]
gamma = matrix(c(ranef(rslt)$clustime,ranef(rslt)$fcluster),ncol=1)
nq = dim(gamma)[1]
X = model.matrix(rslt,type="fixed")
Z = model.matrix(rslt,type="random")
Y = rslt@resp$y
eta = predict(rslt,type="link")
ginv_eta = predict(rslt,type="response")

# Note that delta = Identity matrix for identity link
delta = diag(nobs(rslt)) # nobs is number of observations = 5448 in the ctdata
deltainv = delta # in this case

#
P = deltainv%*%(Y-ginv_eta) + eta # In this case P = Y, since ginv_eta + eta = 0
e = matrix(P - X%*%beta,ncol=1) # make it a matrix with 1 column
XtVX = vcov(rslt) # model based variance
theta = as.data.frame(VarCorr(rslt)) 
R = diag(c(rep(theta$vcov[1],ngrps(rslt)["clustime"]),rep(theta$vcov[2],ngrps(rslt)["fcluster"]))) #R(theta)
G = ngrps(rslt)["fcluster"]
sum=matrix(0,np,np)

start_time_WB_P <- proc.time()
for (g in 1:G){
  grp = attr(rslt,"frame")[,"fcluster"] == g
  ng = sum(grp)
  Sigma = theta$vcov[theta$grp=="Residual"]*diag(ng)
  
  WB_A <- diag(1/diag(deltainv[grp,grp]%*%Sigma%*%deltainv[grp,grp])) # this is diagonal, which is the first term of WB
  
  WB_U <- Z[grp,] # refer U,C,V in the definition of Woodbury formula: W_UCV = Z[grp,]%*%R%*%t(Z[grp,])
  
  WB_C <- solve(R)
  
  WB_V <- t(Z[grp,])
  
  #Vinv = solve(V) # Use Woodbury to compute Vinv
  
  # Create a WoodburyMatrix object
  W <- WoodburyMatrix(A = WB_A, B = WB_C, U = WB_U, V = WB_V)
  
  # Compute the inverse
  Vinv <- solve(W)
  
  H = X[grp,]%*%XtVX%*%t(X[grp,])%*%Vinv
  Q = t(X[grp,])%*%Vinv%*%X[grp,]%*%XtVX
  F = diag(ng) # for now
  A = diag(np) # for now
  
  sum = sum + A%*%t(X[grp,])%*%Vinv%*%t(F)%*%e[grp,]%*%t(e[grp,])%*%F%*%Vinv%*%X[grp,]%*%A
  
}
end_time_WB_P <- proc.time() - start_time_WB_P


c = 1 # for now
robustVar = c*XtVX%*%sum%*%XtVX
diag(robustVar) # Our result

end_time_WB_P

# > diag(rvar0) # classic Liang-Zeger
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0001379472 0.0002188822 0.0002572821 0.0003607838 0.0003631747 0.0001417230
