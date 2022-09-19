
x<-seq(0,20,0.1)
y <- sin(x) + rnorm(x)
plot(x,y)




x<-seq(0,20,0.1)
y <- sin(x) + rnorm(x)

data<-data.frame(x,y)

ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=loocv_smoother(x,y,seq(0,1,0.1))),color="blue")


wave <- seq(0, 10, 0.2)
y <- sin(wave) + rnorm(wave)
data<-data.frame(x=wave,y=y)

#plot(x,y)

#wave <- seq(0, 10, 0.2)
#y <- sin(wave) + rnorm(wave)

smoother_loocv_v1(data$x,data$y)


#data<-data.frame(x=wave,y=y)

ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=smoother_loocv_v1(x,y)),color="blue")


ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=loocv_smoother(x,y, seq(0, 1, 0.1))),color="blue")


ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=smoother_loocv_v1(x,y)),color="blue")

ggplot(data, aes(x = x, y = y)) + 
  geom_point()+
  geom_line(aes(y=loocv_smoother(x,y,seq(0,1,0.1))),color="blue")

plot(wave, y)
lines(wave, l_m$S %*% y)
lines(wave, smoother_loocv_v1(wave,y))


ggplot(data, aes(x = Year, y = Temperature)) + 
  geom_point()+
  geom_line(aes(y=smoother_loocv_v1(data$Year,data$Temperature)),color="blue")

ggplot(Nuuk_year, aes(x = year, y = mean)) + 
  geom_point() + geom_smooth(se = FALSE) + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 10), color = "red", se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"), color = "purple", se = FALSE)+
  geom_line(aes(y=smoother_loocv_v1(Nuuk_year$year,Nuuk_year$mean)),color="blue")



x<-seq(0,20,0.1)
y <- sin(x) + rnorm(x)
data<-data.frame(x=x,y=y)

###