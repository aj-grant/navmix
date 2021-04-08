library(MASS)
library(fossil)
library(mclust)
library(expm)
library(parallel)
library(skmeans)
library(truncnorm)

cl = makeCluster(6)
clusterEvalQ(cl, library(MASS))
clusterEvalQ(cl, library(fossil))
clusterEvalQ(cl, library(mclust))
clusterEvalQ(cl, library(skmeans))
clusterEvalQ(cl, library(truncnorm))
clusterEvalQ(cl, library(expm))

beta1_sim = function(n, low, upp, s_noise){
  K = length(n) - 1
  b = vector(length = sum(n))
  a = 1
  for (j in 1:K){
    b[a:(a+n[j]-1)] = runif(n[j], low[j], upp[j])
    a = a + n[j]
  }
  b[a:(a+n[(K+1)]-1)] = rnorm(n[(K+1)], 0, s_noise)
  return(b)
}
beta2_sim = function(n, d, th_bd, s_clust, s_noise, beta1){
  K = length(n) - 1
  b = vector(length = sum(n))
  a = 1
  for (j in 1:K){
    b[a:(a+n[j]-1)] = tan(rtruncnorm(n[j], a = d[j] - th_bd, b = d[j] + th_bd, mean = d[j], sd = s_clust)) * beta1[a:(a+n[j]-1)]
    a = a + n[j]
  }
  b[a:(a+n[(K+1)]-1)] = rnorm(n[(K+1)], 0, s_noise)
  return(b)
}

s_noise = 0.2
s_clust0 = list('m2' = 0.2, 'm9' = 0.2)
th_bd0 = list('m2' = 0.2, 'm9' = 0.2)

n0 = list('k1' = c(100, 20), 'k2' = c(50, 50, 20), 'k4' = c(40, 20, 20, 20, 20), 'k8' = c(30, 10, 10, 10, 10, 10, 10, 10, 20))

N0 = list('1sam' = 20000,
          '2sam' = list(m2 = c(20000, 40000),
                        m9 = c(20000, 23000, 26000, 29000, 32000, 35000, 38000,
                               41000, 44000)))

d_m2 = list('k1' = matrix(pi / c(4), nrow = 1),
            'k2' = matrix(pi / c(4, -4), nrow = 2),
            'k4' = matrix(pi / c(4, -4, 4, -4), nrow = 4))

d_m9 = list('k1' = matrix(pi / rep(4, 8), nrow = 1),
            'k2' = matrix(pi / c(rep(c(4, -4), 2), rep(c(8, 8), 4), rep(c(8, -8), 2)), nrow = 2),
            'k4' = matrix(pi / c(rep(c(4, 4, -4, -4), 2), rep(c(8, 4, 8, -4), 2), rep(c(8, Inf, 8, Inf), 2), rep(c(8, Inf, -8, Inf), 2)), nrow = 4)
            )


M = 1000

clusterExport(cl, c('beta1_sim', 'beta2_sim', 's_clust0', 'th_bd0', 's_noise', 'n0', 'N0', 'd_m2', 'd_m9', 'sstat', 'M'))

#1 sample sims
clusterSetRNGStream(cl, 20201201)
D_k1_m2_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k1'
  K = length(n) - 1
  d = d_m2$'k1'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

clusterSetRNGStream(cl, 20201202)
D_k2_m2_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k2'
  K = length(n) - 1
  d = d_m2$'k2'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

clusterSetRNGStream(cl, 20201203)
D_k4_m2_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k4'
  K = length(n) - 1
  d = d_m2$'k4'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

clusterSetRNGStream(cl, 20201204)
D_k1_m9_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k1'
  K = length(n) - 1
  d = d_m9$'k1'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

clusterSetRNGStream(cl, 20201205)
D_k2_m9_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k2'
  K = length(n) - 1
  d = d_m9$'k2'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

clusterSetRNGStream(cl, 20201206)
D_k4_m9_s1 = parLapply(cl, 1:M, function(i){
  n = n0$'k4'
  K = length(n) - 1
  d = d_m9$'k4'
  m = dim(d)[2] + 1
  N = N0$'1sam'
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data
  maf = runif(sum(n), 0.01, 0.5)
  G = sapply(1:sum(n), function(j){rbinom(N, 2, maf[j])})
  gu = runif(m, -0.8, 0.8)
  U = rnorm(N, 0, 1)
  e = sapply(1:m, function(k){
    gu[k] * U + sqrt(1 - gu[k]^2) * rnorm(N, 0, 1)
  })
  X = G %*% B + e
  #Compute summary statistics
  bhat = se = matrix(nrow = sum(n), ncol = m)
  for (j in 1:sum(n)){
    for (k in 1:m) {
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  corX = cor(X)
  return(list('bhat' = bhat, 'se' = se, 'corX' = corX))
})

########################################################################################################################
#Results
clusterExport(cl, c('D_k1_m2_s1', 'D_k2_m2_s1', 'D_k4_m2_s1', 'D_k1_m9_s1', 'D_k4_m9_s1', 'D_k2_m9_s1',
                    'gxclust_K', 'row_norm', 'col_norm', 'C_vMF', 'f_vMF', 'gxclust', 'gxclust_K'))

########################################################################################################################
#Cluster 1 sample datasets
clusterSetRNGStream(cl, 20201213)
R_k1_m2_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m2_s1[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201214)
R_k2_m2_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m2_s1[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201215)
R_k4_m2_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m2_s1[[j]]
  z0 = c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201216)
R_k1_m9_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m9_s1[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201217)
R_k2_m9_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m9_s1[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201218)
R_k4_m9_s1 = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m9_s1[[j]]
  z0 = c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

save('R_k1_m2_s1', 'R_k2_m2_s1', 'R_k4_m2_s1', 'R_k1_m9_s1', 'R_k4_m9_s1', 'R_k2_m9_s1', file = 'sims_res_s1.R')

########################################################################################################################
#Cluster 1 sample datasets with estimated correlation
clusterSetRNGStream(cl, 20201213)
R_k1_m2_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m2_s1[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201214)
R_k2_m2_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m2_s1[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201215)
R_k4_m2_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m2_s1[[j]]
  z0 =c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201216)
R_k1_m9_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m9_s1[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201217)
R_k2_m9_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m9_s1[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201218)
R_k4_m9_s1_cor = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m9_s1[[j]]
  z0 = c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = matrix(nrow = nrow(D$bhat), ncol = ncol(D$bhat))
  for (i in 1:n){
    S1 = diag(D$se[i, ])
    S = S1 %*% D$corX %*% S1
    B_std[i, ] = solve(sqrtm(S), D$bhat[i, ])
  }
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

########################################################################################################################
#2 sample sims
set.seed(20201207)
D_k1_m2_s2 = lapply(1:M, function(i){
  n = n0$'k1'
  K = length(n) - 1
  d = d_m2$'k1'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m2
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

set.seed(20201208)
D_k2_m2_s2 = lapply(1:M, function(i){
  n = n0$'k2'
  K = length(n) - 1
  d = d_m2$'k2'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m2
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

set.seed(20201209)
D_k4_m2_s2 = lapply(1:M, function(i){
  n = n0$'k4'
  K = length(n) - 1
  d = d_m2$'k4'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m2
  s_clust = s_clust0$'m2'
  th_bd = th_bd0$'m2'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

set.seed(20201210)
D_k1_m9_s2 = lapply(1:M, function(i){
  n = n0$'k1'
  K = length(n) - 1
  d = d_m9$'k1'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m9
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

set.seed(20201211)
D_k2_m9_s2 = lapply(1:M, function(i){
  n = n0$'k2'
  K = length(n) - 1
  d = d_m9$'k2'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m9
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

set.seed(20201212)
D_k4_m9_s2 = lapply(1:M, function(i){
  n = n0$'k4'
  K = length(n) - 1
  d = d_m9$'k4'
  m = dim(d)[2] + 1
  N = N0$'2sam'$m9
  s_clust = s_clust0$'m9'
  th_bd = th_bd0$'m9'
  #Generate genetic variant-trait associations
  beta1 = beta1_sim(n, c(0.05, 0.05, -0.4, -0.4), c(0.4, 0.4, -0.05, -0.05), s_noise)
  B = matrix(nrow = sum(n), ncol = m)
  B[, 1] = beta1
  for (k in 2:m){
    B[, k] = beta2_sim(n, d[, (k-1)], th_bd, s_clust, s_noise, beta1)
  }
  #Generate indiviudal level data & summary statistics
  maf = runif(sum(n), 0.01, 0.5)
  bhat = se = matrix(nrow = sum(n), ncol = m)
  gu = runif(m, -0.8, 0.8)
  for (k in 1:m){
    G = sapply(1:sum(n), function(j){rbinom(N[k], 2, maf[j])})
    U = rnorm(N[k], 0, 1)
    e = sapply(1:m, function(j){
      gu[j] * U + sqrt(1 - gu[j]^2) * rnorm(N[k], 0, 1)
    })
    X = G %*% B + e
    for (j in 1:sum(n)){
      mod = sstat(X[, k], G[, j])
      bhat[j, k] = mod$bhat
      se[j, k] = mod$se
    }
  }
  return(list('bhat' = bhat, 'se' = se))
})

clusterExport(cl, c('D_k1_m2_s2', 'D_k2_m2_s2', 'D_k4_m2_s2', 'D_k1_m9_s2', 'D_k4_m9_s2', 'D_k2_m9_s2'))

########################################################################################################################
#Cluster 2 sample datasets
clusterSetRNGStream(cl, 20201213)
R_k1_m2_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m2_s2[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201214)
R_k2_m2_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m2_s2[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201215)
R_k4_m2_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m2_s2[[j]]
  z0 = c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201216)
R_k1_m9_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k1_m9_s2[[j]]
  z0 = c(rep(1, 100))
  Kmax = 4
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201217)
R_k2_m9_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k2_m9_s2[[j]]
  z0 = c(rep(1, 50), rep(2, 50))
  Kmax = 5
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})

clusterSetRNGStream(cl, 20201218)
R_k4_m9_s2 = parSapply(cl, 1:M, function(j){
  #
  D = D_k4_m9_s2[[j]]
  z0 = c(rep(1, 40), rep(2, 20), rep(3, 20), rep(4, 20))
  Kmax = 7
  #
  n = dim(as.matrix(D$bhat))[1]
  B_std = D$bhat / D$se
  z_prop = gxclust_K(B_std, K = Kmax, pj_ini = 0.05)
  z_mclust = Mclust(B_std, G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_mclust_prop = Mclust(row_norm(B_std), G = 1:Kmax, verbose = FALSE, initialization = list(noise = sample(n, 5)))
  z_prop_mem = as.numeric(z_prop$clust$z)
  z_mclust_mem = as.numeric(z_mclust$classification)
  z_mclust_prop_mem = as.numeric(z_mclust_prop$classification)
  r_prop = rand.index(z0, z_prop_mem[1:100])
  r_mclust = rand.index(z0, z_mclust_mem[1:100])
  r_mclust_prop = rand.index(z0, z_mclust_prop_mem[1:100])
  return(list(rind = c(r_prop, r_mclust, r_mclust_prop), noclust = c(ncol(z_prop$clust$g)-1, max(z_mclust_mem),
                                                                     max(z_mclust_prop_mem)),
              nonoise = c(sum(z_prop_mem == max(z_prop_mem)), sum(z_mclust_mem == 0), sum(z_mclust_prop_mem == 0)),
              z_prop$clust$g, z_mclust$z, z_mclust_prop$z))
})
