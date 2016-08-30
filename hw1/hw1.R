library(ggplot2)
library(microbenchmark)
library(reshape2)

inversion_method = function(X, y, W) {
  A = solve(t(X) %*% W %*% X)
  return(A %*% t(X) %*% W %*% y)
}

smart_inversion_method = function(X, y, W) {
  # Uses just the diagonal of W to not multiply all the zeroes.
  w = diag(W)
  A = solve(t(X) %*% (w * X))
  return(A %*% t(X) %*% (w * y))
}

my_method = function(X, y, W) {
  # TODO: update
  # from http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
  w = diag(W)
  A = t(X) %*% (w * X)
  R = chol(A)
  L = t(R)
  z = forwardsolve(L, t(X) %*% (w * y))
  return(backsolve(R, z))
}


# my_check <- function(values) {
#   # In cases where the values should all be near-identical
#   bhat = values[[1]]
#   all(sapply(values[-1], function(x) {sum((x-bhat)^2) < 1e-9}))
# }
# 
# 
# N = 2000
# P = 1000
# X = matrix(rnorm(N * P), nrow=N)
# W = diag(nrow=N)
# y = rnorm(N)
# # Ensure correctness
# microbenchmark(inversion_method(X, y, W), smart_inversion_method(X, y, W), my_method(X, y, W), times=1, control=list(order="inorder"), check=my_check)
# 
# # Better timing
# microbenchmark(compare(inversion_method, N, P), compare(smart_inversion_method, N, P), compare(my_method, N, P), times=100)
# 
 get_time = function() {
   return(unname(proc.time()["elapsed"]))
 }

# compare = function(method, n, p) {
#   X = matrix(rnorm(N * P), nrow=N)
#   W = diag(nrow=N)
#   y = rnorm(N)
#   return(method(X, y, W))
# }

compare_fairly = function(n, p, iter=100) {
  itime = 0
  stime = 0
  mtime = 0
  for (i in 1:iter) {
    X = matrix(rnorm(n * p), nrow=n)
    W = diag(nrow=n)
    y = rnorm(n)
    
    t = get_time()
    inversion_method(X, y, W)
    itime = itime + (get_time() - t)
    
    t = get_time()
    smart_inversion_method(X, y, W)
    stime = stime + (get_time() - t)
    
    t = get_time()
    my_method(X, y, W)
    mtime = mtime + (get_time() - t)
  }
  return(data.frame(itime=c(itime/iter),
              stime=c(stime/iter),
              mtime=c(mtime/iter)))
}

compute_fairly_batch = function(ns, ps, iter=100) {
  results = data.frame()
  for (n in ns) {
    for (p in ps) {
      results = rbind(results, cbind(compare_fairly(n,p), list(n=n, p=p)))
      print(paste(n, p))
    }
  }
  return(results)
}



results = compute_fairly_batch(c(3200), c(10, 20, 40, 80, 160, 320, 640, 1280))

names(results) = c("inverse", "smart inverse", "my method", "n", "p")

melted = melt(results, id.vars = c("n", "p"))
names(melted)[3] = "method"
ggplot(melted, aes(p, value, col=method)) +
  geom_line() + geom_point() + 
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10)) + scale_x_log10() + 
  xlab("p") + ylab("average time (s)") + ggtitle("n=3200")






lm(log(`inverse`) ~ log(p), data=results)
lm(log(`smart inverse`) ~ log(p), data=results)
lm(log(`my method`) ~ log(p), data=results)





# Part (D)
library(Matrix)
N = 2000
P = 1000
S = 0.05
X = matrix(rnorm(N * P), nrow=N)
X = X * matrix(rbinom(N * P, 1, S), nrow=N)
X = subset(melt(X), value != 0)
X = sparseMatrix(
  X$Var1,
  X$Var2,
  x=X$value)
X[1:5, 1:5]

my_method_sparse = function(X, y, W) {
  myDims = dim(X)
  X = subset(melt(X), value != 0)
  X = sparseMatrix(
    X$Var1,
    X$Var2,
    x=X$value,
    dims=myDims)
  # from http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
  # Need t(X) %*% W %*% X, so break up into [t(X) %*% t(W^(1/2))] %*% [W^(1/2) %*% X]
  w = diag(W)^0.5
  # does t(A) %*% A faster
  A = crossprod(w * X)
  R = chol(A)
  L = t(R)
  z = forwardsolve(L, t(X) %*% (w * y))
  return(backsolve(R, z))
}

compare_fairly_sparse = function(n, p, s, iter=100) {
  itime = 0
  stime = 0
  mtime = 0
  sparsetime = 0
  for (i in 1:iter) {
    X = matrix(rnorm(n * p), nrow=n)
    W = diag(nrow=n)
    y = rnorm(n)
    X = X * matrix(rbinom(n * p, 1, s), nrow=n)
    
    t = get_time()
    inversion_method(X, y, W)
    itime = itime + (get_time() - t)
    
    t = get_time()
    smart_inversion_method(X, y, W)
    stime = stime + (get_time() - t)
    
    t = get_time()
    my_method(X, y, W)
    mtime = mtime + (get_time() - t)
    
    t = get_time()
    my_method_sparse(X, y, W)
    sparsetime = sparsetime + (get_time() - t)
  }
  return(data.frame(itime=c(itime/iter),
                    stime=c(stime/iter),
                    mtime=c(mtime/iter),
                    sparsetime=c(sparsetime/iter)))
}

compute_fairly_sparse_batch = function(ns, ps, ss, iter=100) {
  results = data.frame()
  for (n in ns) {
    for (p in ps) {
      for (s in ss) {
        results = rbind(results, cbind(compare_fairly_sparse(n,p, s), list(n=n, p=p, s=s)))
        print(paste(n, p, s))
      }
    }
  }
  return(results)
}

results = compute_fairly_sparse_batch(c(3200), c(10, 20, 40, 80, 160, 320, 640, 1280), c(0.05), iter = 10)
results
names(results) = c("inverse", "smart inverse", "my method", "my sparse method", "n", "p", "s")

melted = melt(results, id.vars = c("n", "p", "s"))
names(melted)[4] = "method"
ggplot(melted, aes(p, value, col=method)) +
  geom_line() + geom_point() + 
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10)) + scale_x_log10() + 
  xlab("p") + ylab("average time (s)") + ggtitle("sparsity=0.05, n=3200")











rep = microbenchmark(
  solve(matrix(rnorm(5^2), ncol=5)),
  solve(matrix(rnorm(10^2), ncol=10)),
  solve(matrix(rnorm(20^2), ncol=20)),
  solve(matrix(rnorm(40^2), ncol=40)),
  solve(matrix(rnorm(80^2), ncol=80)),
  solve(matrix(rnorm(160^2), ncol=160)),
  solve(matrix(rnorm(320^2), ncol=320)),
  solve(matrix(rnorm(640^2), ncol=640)),
  times=50
)
n = c(5,10,20,40,80,160,320,640)
y = c(68.86962, 85.40070, 144.90512, 404.46966, 1854.11608,
      11909.96368, 79218.60874, 586615.34966)
y = c(50.4660,68.5045,130.7075,390.3295,1831.7670,11323.6575,
      78072.8670,583686.5095)
summary(lm(log(y) ~ log(n)))
proc.time()["user.child"]
  

  
  
  