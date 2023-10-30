source("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OtherMethods/store.R")
num_obs = 50
sigma = 1
dist_noise = "normal"
i = 1
setwd("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data")
img = read.csv(paste0("img_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
X = read.csv(paste0("cov_mat_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
X = as.matrix(X)
voxels = read.csv(paste0("voxels_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))

Y = list()
for (i in 1:dim(img)[2]) {
  Yi = array(dim = c(32,32,8))
  Yi[as.matrix(voxels)] = img[,i]
  Y[[i]] = Yi
}


Para = c(2,2,2)
K = 8
out = mystor(Y, X, Para, K)
STOREB = out$B.hat
STOREbeta = matrix(nrow = 8192, ncol = 3)
for (j in 1:3){
  STOREbeta[,j] = as.vector(STOREB[,,,j])
}
setwd("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data")
MSE_STORE = c()
for (num_obs in c(50)) {
  for (sigma in c(0.5,1, 2, 4, 6)) {
    for (dist_noise in c("normal")) {
      for (i in 1:5) {
        
        img = read.csv(paste0("img_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
        X = read.csv(paste0("cov_mat_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
        X = as.matrix(X)
        voxels = read.csv(paste0("voxels_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
        true_beta = read.csv(paste0("true_beta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
        true_beta = as.matrix(true_beta)
        
        Y = list()
        for (ni in 1:dim(img)[2]) {
          Yi = array(dim = c(32,32,8))
          Yi[as.matrix(voxels)] = img[,ni]
          Y[[ni]] = Yi
        }
        
        
        Para = c(2,2,2)
        K = 8
        out = mystor(Y, X, Para, K)
        STOREB = out$B.hat
        STOREbeta = matrix(nrow = 8192, ncol = 3)
        for (j in 1:3){
          STOREbeta[,j] = as.vector(STOREB[,,,j])
        }
        write.csv(STOREbeta, paste0("STOREbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"), row.names = FALSE)
        
        MSE = mean((true_beta - STOREbeta)^2)
        MSE_STORE = c(MSE_STORE, MSE)
        
      }
    }
  }
}
data_STORE = data.frame(num_obs = rep(c(50, 100), each=25), 
                        sigma = rep(rep(c(0.5,1,2,4,6), each=5), 2),MSE = MSE_STORE)
write.csv(data_STORE, "./data_STORE.csv", row.names=F)






img = read.csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/tuneTest/sigma2/img_num_obs_50.csv")
X = read.csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/tuneTest/sigma2/cov_mat_num_obs_50.csv")
X = as.matrix(X)
voxels = read.csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/tuneTest/voxels20208.csv")
#true_beta = read.csv(paste0("true_beta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
#true_beta = as.matrix(true_beta)

Y = list()
for (ni in 1:dim(img)[2]) {
  Yi = array(dim = c(20,20,8))
  Yi[as.matrix(voxels)] = img[,ni]
  Y[[ni]] = Yi
}


Para = c(2,2,2)
K = 8
out = mystor(Y, X, Para, K)
STOREB = out$B.hat
STOREbeta = matrix(nrow = 3200, ncol = 3)
for (j in 1:3){
  STOREbeta[,j] = as.vector(STOREB[,,,j])
}
write.csv(STOREbeta, "/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/tuneTest/sigma2/STOREbeta.csv", row.names = FALSE)

