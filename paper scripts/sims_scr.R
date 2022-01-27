library(fossil)
library(mclust)
library(expm)
library(tictoc)
library(parallel)
library(navmix)

cl = makeCluster(6)
clusterEvalQ(cl, library(mclust))
clusterEvalQ(cl, library(navmix))
clusterEvalQ(cl, library(expm))

N0 = list('m2' = rep(20000, 2), 'm9' = rep(20000, 9))
n0 = list('k1' = c(80, 20), 'k2' = c(40, 40, 20), 'k4' = c(30, 20, 20, 10, 20))
g0 = list('g0' = 0, 'g4' = 0.4, 'g8' = 0.8)
d_m2 = list('k1' = c(1, 1),
            'k2' = rbind(c(1, 1), c(1, -1)),
            'k4' = rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)))
d_m9 = list('k1' = c(1, 1, 1, rep(0.5, 6)),
            'k2' = rbind(c(1, 1, 1, rep(0.5, 6)), c(1, -1, -1, rep(0.5, 4), -0.5, -0.5)),
            'k4' = rbind(c(1, 1, 1, rep(0.5, 6)), c(rep(1, 5), rep(0, 4)), c(-1, -1, -1, rep(0.5, 4), -0.5, -0.5), c(rep(-1, 5), rep(0, 4))))

sstat = function(Y, X, intercept = TRUE){
  n = length(Y)
  if (intercept == TRUE){
    xx = cbind(rep(1, n), X)
  }
  else {xx = X}
  mod = lm.fit(xx, Y)
  bhat= c(mod$coefficients[2])
  s = t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)
  se = sqrt((c(s) * solve(t(xx) %*% xx))[2,2])
  return(list("bhat" = bhat, "se" = se))
}

M = 1000

clusterExport(cl, c('N0', 'n0', 'g0', 'd_m2', 'd_m9', 'sstat', 'M'))

########################################################################################################################
#gamma = 0
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210723)
D_g0_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g0'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m2$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210724)
D_g0_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g0'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m2$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210725)
D_g0_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g0'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m2$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

#m = 9
clusterSetRNGStream(cl, 20210726)
D_g0_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g0'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m9$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210727)
D_g0_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g0'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m9$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210728)
D_g0_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g0'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m9$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

########################################################################################################################
#gamma = 0.4
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210723)
D_g4_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g4'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m2$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210724)
D_g4_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g4'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m2$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210725)
D_g4_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g4'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m2$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

#m = 9
clusterSetRNGStream(cl, 20210726)
D_g4_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g4'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m9$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210727)
D_g4_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g4'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m9$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210728)
D_g4_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g4'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m9$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

########################################################################################################################
#gamma = 0.8
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210723)
D_g8_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g8'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m2$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210724)
D_g8_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g8'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m2$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210725)
D_g8_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m2'
  g = g0$'g8'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m2$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

#m = 9
clusterSetRNGStream(cl, 20210726)
D_g8_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g8'
  n = n0$'k1'
  K = length(n) - 1
  d = matrix(d_m9$'k1', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210727)
D_g8_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g8'
  n = n0$'k2'
  K = length(n) - 1
  d = matrix(d_m9$'k2', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

clusterSetRNGStream(cl, 20210728)
D_g8_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  N = N0$'m9'
  g = g0$'g8'
  n = n0$'k4'
  K = length(n) - 1
  d = matrix(d_m9$'k4', nrow = K)
  #
  n_all = sum(n)
  n_sub = sum(n[1:K])
  m = length(d) / K
  
  b = matrix(rep(0, n_sub * K), ncol = K)
  for (k in 1:K){
    p2 = runif(1, 0.05, 0.2)
    b2 = cbind(runif(n[k], 0.03, 0.06), rnorm(n[k], 0.1, 0.02))
    b[(sum(n[0:(k-1)])+1):sum(n[0:k]), k] = sapply(1:n[k], function(j){sample(b2[j, ], 1, prob = c(1-p2, p2))})
  }
  a = sapply(1:m, function(l){runif(n[K+1], -0.1, 0.1)})
  N_max = max(N)
  
  U = rnorm(N_max, 0, 1)
  maf = runif(n_all, 0.01, 0.5)
  G = sapply(1:n_all, function(j){rbinom(N_max, 2, maf[j])})
  L = sapply(1:K, function(k){G[, 1:n_sub] %*% b[, k]})
  X = sapply(1:m, function(l){L %*% d[, l] + G[, (n_sub+1):n_all] %*% a[, l] + g * U + sqrt(1 - g^2) * rnorm(N_max, 0, 1)})
  corX = cor(X)
  
  bhat = se = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    for (l in 1:m){
      Xmod = sstat(X[1:N[l], l], G[1:N[l], j])
      bhat[j, l] = Xmod$bhat
      se[j, l] = Xmod$se
    }
  }
  
  b_std = bhat / se
  b_prop = navmix::row_norm(b_std)
  b_cor = matrix(nrow = n_all, ncol = m)
  for (j in 1:n_all){
    S1 = diag(se[j, ])
    S = S1 %*% corX %*% S1
    b_cor[j, ] = solve(sqrtm(S), bhat[j, ])
  }
  
  return(list(b_std = b_std, b_prop = b_prop, b_cor = b_cor))
})

########################################################################################################################
#Results
########################################################################################################################

clusterExport(cl, c('D_g0_m2_k1', 'D_g0_m2_k2', 'D_g0_m2_k4', 'D_g0_m9_k1', 'D_g0_m9_k2', 'D_g0_m9_k4',
                    'D_g4_m2_k1', 'D_g4_m2_k2', 'D_g4_m2_k4', 'D_g4_m9_k1', 'D_g4_m9_k2', 'D_g4_m9_k4',
                    'D_g8_m2_k1', 'D_g8_m2_k2', 'D_g8_m2_k4', 'D_g8_m9_k1', 'D_g8_m9_k2', 'D_g8_m9_k4'))

########################################################################################################################
#gamma = 0
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g0_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210824)
R_g0_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210825)
R_g0_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g0_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210827)
R_g0_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210828)
R_g0_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

########################################################################################################################
#gamma = 4
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g4_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k1
  n = n0$k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210824)
R_g4_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k2
  n = n0$k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210825)
R_g4_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k4
  n = n0$k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g4_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k1
  n = n0$k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210827)
R_g4_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k2
  n = n0$k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210828)
R_g4_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k4
  n = n0$k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

########################################################################################################################
#gamma = 0.8
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g8_m2_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k1
  n = n0$k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210824)
R_g8_m2_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k2
  n = n0$k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210825)
R_g8_m2_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k4
  n = n0$k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g8_m9_k1 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k1
  n = n0$k1
  K = 1
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210827)
R_g8_m9_k2 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k2
  n = n0$k2
  K = 2
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

clusterSetRNGStream(cl, 20210828)
R_g8_m9_k4 = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k4
  n = n0$k4
  K = 4
  #
  b_std = D[[i]]$b_std
  b_prop = D[[i]]$b_prop
  b_cor = D[[i]]$b_cor
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z))
})

########################################################################################################################
#Results: variants with > 1 sig association only
########################################################################################################################

########################################################################################################################
#gamma = 0
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g0_m2_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210824)
R_g0_m2_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210825)
R_g0_m2_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m2_k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g0_m9_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210827)
R_g0_m9_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210828)
R_g0_m9_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g0_m9_k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

########################################################################################################################
#gamma = 0.4
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g4_m2_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k1
  n = n0$k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210824)
R_g4_m2_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k2
  n = n0$k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210825)
R_g4_m2_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m2_k4
  n = n0$k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g4_m9_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k1
  n = n0$k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210827)
R_g4_m9_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k2
  n = n0$k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210828)
R_g4_m9_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g4_m9_k4
  n = n0$k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

########################################################################################################################
#gamma = 0.8
########################################################################################################################
#m = 2
clusterSetRNGStream(cl, 20210823)
R_g8_m2_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k1
  n = n0$k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210824)
R_g8_m2_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k2
  n = n0$k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210825)
R_g8_m2_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m2_k4
  n = n0$k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

#m = 9
clusterSetRNGStream(cl, 20210826)
R_g8_m9_k1_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k1
  n = n0$k1
  K = 1
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210827)
R_g8_m9_k2_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k2
  n = n0$k2
  K = 2
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

clusterSetRNGStream(cl, 20210828)
R_g8_m9_k4_sig = parLapply(cl, 1:M, function(i){
  #
  D = D_g8_m9_k4
  n = n0$k4
  K = 4
  #
  sig = which(apply(abs(D[[i]]$b_std), 1, max)>qnorm(1-5e-4))
  b_std = D[[i]]$b_std[sig, ]
  b_prop = D[[i]]$b_prop[sig, ]
  b_cor = D[[i]]$b_cor[sig, ]
  p = dim(b_std)[1]
  g_navmix_std = navmix(b_std)
  g_navmix_cor = navmix(b_cor)
  g_mclust_std = Mclust(b_std, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  g_mclust_prop = Mclust(b_prop, initialization = list(noise = sample(p, 5)), verbose = FALSE)
  
  return(list(g_navmix_std = g_navmix_std$fit$g, g_navmix_cor = g_navmix_cor$fit$g,
              g_mclust_std = g_mclust_std$z, g_mclust_prop = g_mclust_prop$z, sig = sig))
})

stopCluster(cl)
