# http://localhost:3721/web/viewer.html
smooth.spline
splines::splineDesign


# LOOCV

#LOOCV from CSwR
loocv <- function(k,y){
  f_h <- run_mean(y,k)
  mean(( (y-f_hat) / (1-1/k) )^2, na.rm=TRUE)
}



