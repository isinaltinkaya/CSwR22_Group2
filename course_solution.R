library(splines) #for splineDesign
library(fda) # for bsplinepen
library(bestglm) #for LOOCV
#use smooth.spline

#https://cswr.nrhstat.org/bivariate.html?q=spline#splines

pen_mat <- function(inner_knots, r=3) {
  knots <- sort(c(rep(range(inner_knots), r), inner_knots))
  d <- diff(inner_knots)  # The vector of knot differences; b - a 
  g_ab <- splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 #Simpson's rule
}

smoother <- function(y, inner_knots, lambda, r=3){ 
  Phi <- splineDesign(c(rep(range(inner_knots), r), inner_knots), inner_knots)
  Omega <- pen_mat(inner_knots)
  
  Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi) %*% y
  )
}

#todo: generate sine/cosine wave and add noise
#try to fit a smoother and see if we can get close to the original wave
#if we do: carry on i.e implement LOOCV for lambda
#in the future: fit something more "exotic"

wave <- seq(0, 10, 0.2)
y <- sin(wave) + rnorm(wave)
plot(wave, y)

smo <- smoother(y, wave, 0.3)
lines(wave, smo, col='blue')
lines(wave, sin(wave), col='red')
lines(wave, smo - sin(wave))


loocv_smoother <- function(inner_knots, y, lambda){
  n <-  length(lambda)#length(x) == length(y)
  loss <- rep(Inf, n)
  Sm <- list()
  Omega <- pen_mat(inner_knots) #penalty matrix
  Phi <- splineDesign(c(rep(range(inner_knots), 3), inner_knots), inner_knots)
  tphi_phi <- crossprod(Phi)
  #base functions evaluation
  #f_hat = Phi.(Phi^T.Phi + lambda.Omega)^-1.Phi^T.y
  #= S.y
  #loocv = sum(((y[i]-f_hat[i])/(1-S[i,i]))^2)
  for(i in 1:n){
    if(det(tphi_phi + lambda[i] * Omega) != 0){ #make sure we have a solution
      S <- Phi %*% solve(tphi_phi + lambda[i] * Omega) %*% t(Phi)
      #S_lambda matrix
      Sm[[i]] <- S #spares computation time for the return value
      loss[i] <- mean(((y - S * y)/(1 - diag(S)))^2, na.rm = TRUE)
    }
  }
  ind <- which.min(loss) #index of best lambda i.e loss i.e S matrix
  return(list(lambda=lambda, best_lambda=lambda[ind],
              S=Sm[[ind]], loss=loss, index=ind))
}

l_m <- loocv_smoother(wave, y, seq(0, 1, 0.1))
#l_m <- loocv_smoother(wave, y, c(0.3))
print(l_m$loss)
lines(wave, l_m$S %*% y)
