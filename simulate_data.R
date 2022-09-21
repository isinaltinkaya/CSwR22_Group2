
simulate_data<-function(from=0, to=20, step=0.1, signal=sin, noise=rnorm){
  x<-seq(from,to,step)
  data_y<-signal(x)
  y<-data_y + noise(x)
  structure(list(df=data.frame(x=x,y=y), real=data_y),class="simulation_data")
}


#Usage:

#Simulate sin waves
data<-simulate_data(signal=sin)

p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()

p_data +
  geom_line(aes(y=smoother_v1(data$df,data$df$x,data$df$y)$y),
            color="blue")+
  geom_line(aes(y=data$real),
            color="green"
)



#Simulate cos waves
data<-simulate_data(signal=cos)

p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()

p_data +
  geom_line(aes(y=smoother_v1(data$df,data$df$x,data$df$y)$y),
            color="blue")+
  geom_line(aes(y=data$real),
            color="green"
  )

#Simulate logarithmic sequence
data<-simulate_data(signal=log)

p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()

p_data +
  geom_line(aes(y=smoother_v1(data$df,data$df$x,data$df$y)$y),
            color="blue")+
  geom_line(aes(y=data$real),
            color="green"
  )

#Simulate exponential sequence
data<-simulate_data(signal=function(x) 2^x)

p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()

p_data +
  geom_line(aes(y=smoother_v1(data$df,data$df$x,data$df$y)$y),
            color="blue")+
  geom_line(aes(y=data$real),
            color="green"
  )

#Simulate random values
data<-simulate_data(signal=runif)

p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()

p_data +
  geom_line(aes(y=smoother_v1(data$df,data$df$x,data$df$y)$y),
            color="blue")+
  geom_line(aes(y=data$real),
            color="green"
  )

