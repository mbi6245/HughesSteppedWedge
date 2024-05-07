#get family of object, which gives link
family(obj)

#easy way to build R matrix
sigma2 = sigma(obj)^2
lambda = getME(obj,"Lambda")
R = lambda%*%t(lambda)*sigma2