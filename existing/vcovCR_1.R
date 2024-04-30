library("sandwich")
library("clubSandwich")
library("lme4")
#############
ctdata$fcluster = factor(ctdata$hdist)
ctdata$ftime = factor(ctdata$time)
ctdata$clustime = interaction(ctdata$fcluster,ctdata$ftime)
rslt = lmer(ct ~ ftime + treat + (1 | fcluster) + (1 | clustime), data=ctdata)
diag(vcov(rslt))  # model-based
rvar0 = vcovCR(rslt, cluster=ctdata$fcluster, type="CR0")
diag(rvar0) # classic Liang-Zeger
#####################
# write function to reproduce classic Liang-Zeger result
#####################
# Here are some ways of extracting information from the model object
beta=matrix(fixef(rslt),ncol=1)
np=dim(beta)[1]
gamma = matrix(c(ranef(rslt)$clustime,ranef(rslt)$fcluster),ncol=1)
nq=dim(gamma)[1]
X = model.matrix(rslt,type="fixed")
Z = model.matrix(rslt,type="random")
Y = rslt@resp$y # example of "slots" in R
cluster = clubSandwich:::get_outer_group(rslt) # maybe there is a better way to do this
eta = predict(rslt,type="link")
ginv_eta = predict(rslt,type="response")
# Note that delta = Identity matrix for identity link
# more generally, we can use family(rslt)$link to find the link
delta = diag(nobs(rslt))
deltainv = delta # in this case
theta = as.data.frame(VarCorr(rslt))
sigma = theta$vcov[theta$grp=="Residual"]
R = diag(c(rep(theta$vcov[1],ngrps(rslt)["clustime"]),rep(theta$vcov[2],ngrps(rslt)["fcluster"])))
G = ngrps(rslt)["fcluster"]
mbv = vcov(rslt)

diag(vcov(rslt))  # model-based
diag(rvar0) # classic Liang-Zeger
diag(robustVar) # my code