
###

# Computing the vector of knot differences (b - a) 

## we try to use a method that is alternative to diff
## and bytecompiled

library(compiler)

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

###

# FUNCTIONS

#renamed diff_v3 as get_diff for clarity
get_diff <- compiler::cmpfun(function(v) {
  v[2:length(v)]-v[1:(length(v)-1L)]
})


## 

################################################################################
# LOAD LIBRARIES

library(microbenchmark)
library(splines)

data(NuukData)
data("Nuuk_year")


################################################################################
# LOAD DATA


CSwR_dir<-"../CSwR"


Nuuk <- read_table(paste0(CSwR_dir,"/data/nuuk.dat.txt"), 
                   col_names = c("Year", 1:12), 
                   na = "-999", 
                   skip = 1) %>% 
  gather(key = "Month", value = "Temp_Nuuk", -Year, convert = TRUE) %>% 
  mutate(Temp_Nuuk = Temp_Nuuk / 10) %>%
  filter("Year">1866)

Qaqortoq <- read_table(paste0(CSwR_dir,"/data/qaqortoq.dat"),
                       col_names = c("Year", 1:12), 
                       na = "-999", 
                       skip = 1) %>% 
  gather(key = "Month", value = "Temp_Qaqortoq", -Year, convert = TRUE) %>% 
  mutate(Temp_Qaqortoq = Temp_Qaqortoq / 10) %>% 
  filter("Year" > 1866)

p_Nuuk <- ggplot(Nuuk, aes(x = Year, y = Temperature)) + 
  geom_point()

################################################################################


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

pen_mat_v3<-function(inner_knots){
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
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

microbenchmark(times = 100, unit = "ms",pen_mat(x),pen_mat_v2(x),pen_mat_v3(x))
#v3 is the fastest

#we will use v3, renamed as penalty matrix for readability
penalty_matrix<-function(inner_knots){
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
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




# Smoothing function
#
# x       inner knots
# lambda  tuning parameter
smoother_v1 <- function(x,y, lambda){ 
  
  inner_knots<-x
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
  
  Phi <- splineDesign(knots,inner_knots)
  Omega <- penalty_matrix(inner_knots)
  

  Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi) %*% y
  )
}

x<-seq(0,20,0.1)
y <- sin(x) + rnorm(x)
plot(x,y)

data<-data.frame(x=x,y=y)

ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=smoother_v1(x,y,0.5)),color="blue")


loocv_v1 <- function(lambda, y){
  inner_knots<-x
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
  
  Phi <- splineDesign(knots,inner_knots)
  Omega <- penalty_matrix(inner_knots)
  
  
  S <- Phi %*% solve(crossprod(Phi) + lambda * Omega, t(Phi))
  Sii <- diag(S)
  sum(((y - S %*% y) / (1 - Sii))^2, na.rm = TRUE) 
}

lambda <- seq(1,200,1)
LOOCV<-sapply(lambda, loocv_v1, y=data$y)
LOOCV[which.min(LOOCV)]

qplot(lambda,LOOCV)+
  geom_vline(xintercept = LOOCV[which.min(LOOCV)], color="red")


get_best_lambda_loocv_v1 <- function(y, n=200){
  lambda <- seq(1,n,1)
  LOOCV<-sapply(lambda, loocv_v1, y=data$y)
  LOOCV[which.min(LOOCV)]
}


get_best_lambda_loocv_v1(data$y)




smoother_loocv_v1<-function(x,y){
  inner_knots<-x
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots),method="quick")
  
  Phi <- splineDesign(knots,inner_knots)
  Omega <- penalty_matrix(inner_knots)
  
  lambda<-get_best_lambda_loocv_v1(data$y)
  
  Phi %*% solve(
    crossprod(Phi) + lambda * Omega, 
    t(Phi) %*% y
  )
}














### Nuuk data 

Nuuk <- read_table(paste0(CSwR_dir,"/data/nuuk.dat.txt"), 
                    col_names = c("Year", paste0("V",c(1:12))), 
                    na = "-999", 
                    skip = 1)


Nuuk <- read_table(paste0(CSwR_dir,"/data/nuuk.dat.txt"), 
                   col_names = c("Year",1:12), 
                   na = "-999", 
                   skip = 1)


data<-Nuuk%>%
  gather(key = "Month", value = "Temperature", -Year, convert = TRUE) %>% 
  mutate(Temperature = Temperature / 10)

data<-data[data$Year > 1866,]

  
data<-group_by(data,Year)%>%
  summarise(
    Temperature=mean(Temperature))

data<-data[! is.na(data$Temperature),]
#data%>%summarise(Temperature=mean(data$Temperature))
#data <- tibble(x = data$Year, y = mean(data$Temperature))


ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=smoother_loocv_v1(x,y)),color="blue")



ggplot(data, aes(x = Year, y = Temperature)) + 
  geom_point()+
  geom_line(aes(y=smoother_loocv_v1(data$Year,data$Temperature)),color="blue")

