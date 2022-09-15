library(splines) #for splineDesign
library(fda) # for bsplinepen
library(bestglm) #for LOOCV
#use smooth.spline

#https://cswr.nrhstat.org/bivariate.html?q=spline#splines

pen_mat <- function(inner_knots, rep=3) {
  knots <- sort(c(rep(range(inner_knots), rep), inner_knots))
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

smoother <- function(y, inner_knots, lambda, rep=3){ 
  Phi <- splineDesign(c(rep(range(inner_knots), rep), inner_knots), inner_knots)
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

wave <- seq(0, 10, 0.1)
y <- sin(wave) + rnorm(wave)
plot(wave, y)

smo <- smoother(y, wave, 0.3)
lines(wave, smo, add=T)
lines(wave, sin(wave))
lines(wave, smo - sin(wave))
