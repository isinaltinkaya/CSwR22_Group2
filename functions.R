################################################################################
# DEVELOPMENT HISTORY OF FUNCTIONS
################################################################################
# 
# LOAD LIBRARIES

library(microbenchmark)
library(splines)
library(compiler)
library(readr)
library(tidyr)


################################################################################

################################################################################
# LOAD DATA


CSwR_dir<-"../CSwR/"


Nuuk <- read_table(paste0(CSwR_dir,"data/nuuk.dat.txt"), 
                   col_names = c("Year", 1:12), 
                   na = "-999", 
                   skip = 1) %>% 
  gather(key = "Month", value = "Temp_Nuuk", -Year, convert = TRUE) %>% 
  mutate(Temp_Nuuk = Temp_Nuuk / 10)


Nuuk<-Nuuk[Nuuk$Year>1866,]
Nuuk_year <-  group_by(Nuuk, Year) %>% 
  summarise(
    Median = median(Temp_Nuuk),
    High = max(Temp_Nuuk), 
    Low = min(Temp_Nuuk),
    Temperature = mean(Temp_Nuuk)
  )




################################################################################

###

# Computing the vector of knot differences (b - a) 

## we try to use a method that is alternative to diff
## and bytecompiled



diff_v2 <- function(v) {
  v[2:length(v)]-v[1:(length(v)-1L)]
}

#use byte compiling for faster diff()
diff_v3 <- compiler::cmpfun(function(v) {
  v[2:length(v)]-v[1:(length(v)-1L)]
})

diff_v4 <- compiler::cmpfun(function(v) {
  l<-length(v)
  v[2:l]-v[1:(l-1L)]
})

x = rnorm(1000)

all(diff(x)==diff_v2(x))
all(diff(x)==diff_v3(x))
all(diff(x)==diff_v4(x))


microbenchmark(times = 10000, unit = "ms",diff(x),diff_v2(x),diff_v3(x),diff_v4(x))
#thus decided to use diff_v3

xe2 = rnorm(100)
xe3 = rnorm(1000)
xe4 = rnorm(10000)
microbenchmark(times=100, unit="ms", diff(xe2),diff(xe3),diff(xe4),
               diff_v3(xe2),diff_v3(xe3),diff_v3(xe4))

rm(x,xe2,xe3,xe4)
###
################################################################################
# FUNCTIONS

#renamed diff_v3 as get_diff for clarity
get_diff <- compiler::cmpfun(function(v) {
  v[2:length(v)]-v[1:(length(v)-1L)]
})


## 


pen_mat <- function(inner_knots) {
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
  d <- diff(inner_knots)  # The vector of knot differences; b - a 
  g_ab <- splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}



pen_mat_v2<-function(inner_knots){
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
  d <- get_diff(inner_knots)
  g_ab <- splineDesign(knots, inner_knots, derivs=2)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}
x = rnorm(1000)

all(pen_mat(x) == pen_mat_v2(x))

microbenchmark(times = 100, unit = "ms",pen_mat(x),pen_mat_v2(x))
#v2 is faster

#
microbenchmark(times=100, 
               sort(x,method="auto"), # default
               sort(x,method="quick"),
               sort(x,method="radix"),
               sort(x,method="shell"))

# quick outperforms others
# "quick" uses Singleton (1969)'s implementation of Hoare's Quicksort method
# we will use quick

pen_mat_v3<-function(inner_knots, o=4){
  #extend by order-1
  r<-o-1
  
  knots <- sort(c(rep(range(inner_knots), r), inner_knots),
                method="quick")
  d <- get_diff(inner_knots)
  g_ab <- splineDesign(knots, inner_knots, derivs=2,ord=o)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}


microbenchmark(times = 100, unit = "ms",pen_mat(x),pen_mat_v2(x),pen_mat_v3(x))
#v3 is the fastest

#we will use v3, renamed as penalty matrix for readability
penalty_matrix<-function(inner_knots,o=4){
  
  #extend by order-1
  r<-o-1
  
  knots <- sort(c(rep(range(inner_knots), r), inner_knots),
                method="quick")
  d <- get_diff(inner_knots)
  g_ab <- splineDesign(knots, inner_knots, derivs=2,ord=o)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}


# Smoothing function
#
# x       inner knots
# y       data
# lambda  tuning parameter
smoother_v1 <- function(x,y, lambda,o=4){ 
  
  #polynomial order=4 by default
  o<-4
  r<-o-1
  
  inner_knots<-x
  
  Phi <- splineDesign(c(rep(range(inner_knots), r), inner_knots), inner_knots)
  Omega <- penalty_matrix(inner_knots,o=o)
  
  Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi) %*% y
  )
}


loocv_v1 <- function(lambda, y, x){
  
  inner_knots<-x
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
  Phi <- splineDesign(knots,inner_knots)
  
  
  #Phi <- splineDesign(c(rep(range(inner_knots), 3), inner_knots),inner_knots)
  
  Omega <- penalty_matrix(inner_knots)
  
  S <- Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi))
  
  Sii <- diag(S)
  #mean(((y - S %*% y) / (1 - Sii))^2, na.rm = TRUE)
  sum(((y - S %*% y) / (1 - Sii))^2, na.rm = TRUE)
}



choose_lambda_loocv_v1 <- function(y,x,ls=seq(50,250,2)){
  LAMBDAS <- ls
  LOOCV<-sapply(ls, loocv_v1, y=y,x=x)
  ls[which.min(LOOCV)]
  #  qplot(ls, LOOCV) + geom_vline(xintercept = lambda_opt, color = "red")
}


choose_lambda_loocv_v1(y=Nuuk_year$Temperature,x=Nuuk_year$Year)

p_Nuuk <- ggplot(Nuuk_year, aes(x = Year, y = Temperature)) + 
  geom_point()
p_Nuuk + geom_line(aes(
  y = smoother_v1(x=Year,y=Temperature,lambda=choose_lambda_loocv_v1(y=Temperature,x=Year)
  )),
  color = "blue")





smoother_loocv_v1<-function(x,y,...){
  
  #polynomial order=4 by default
  o<-4
  r<-o-1
  
  inner_knots<-x
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),
                method="quick")
  Phi <- splineDesign(knots, inner_knots)
  
  #Phi <- splineDesign(c(rep(range(inner_knots), r), inner_knots), inner_knots)
  Omega <- penalty_matrix(inner_knots,o=o)
  
  lambda<-choose_lambda_loocv_v1(y=y,x=x,ls=seq(50,250,2))
  Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi) %*% y
  )
}



p_Nuuk <- ggplot(Nuuk_year, aes(x = Year, y = Temperature)) + 
  geom_point()
p_Nuuk + geom_line(aes(
  y = smoother_loocv_v1(x=Year,y=Temperature)),
  color = "blue")



