################################################################################
#
# BENCHMARKING
#

source("shared.R")
source("functions.R")
source("simulate_data.R")



################################################################################
# LOAD DATA

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

data<-structure(list(df=data.frame(x=Nuuk_year$Year,y=Nuuk_year$Temperature)),class="simulation_data")


p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()


p_data +
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y)$y),
            color="yellow")+
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y,useSVD = TRUE)$y),
            color="orange")+
  geom_line(aes(y=smooth.spline(data$df$x,data$df$y,all.knots=TRUE)$y),
            color="red")+
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y,useSVD = TRUE, subsample=TRUE)$y),
            color="darkred")



range(smooth.spline(data$df$x,data$df$y,all.knots=TRUE)$y - smoother_v1(data$df$x,data$df$y,useSVD = TRUE, subsample=TRUE)$y)

#smoother_v2(data,useSVD=FALSE, subsample = FALSE)

mb<-microbenchmark(smoother_v1(data$df$x,data$df$y),
               smoother_v2(data,useSVD=FALSE, subsample = TRUE),
               smoother_v2(data,useSVD=TRUE, subsample = TRUE),
               smoother_v2(data,useSVD=TRUE,subsample = FALSE),
               smooth.spline(data$df$x,data$df$y),
               times=20)
autoplot(mb)


################################################################################


data<-simulate_data()


p_data<-ggplot(data$df, aes(x=x,y=y)) + 
  geom_point()


smoother_v1(data$df$x,data$df$y,useSVD=FALSE, subsample=TRUE)
smoother_v1(data$df$x,data$df$y,useSVD=TRUE, subsample=FALSE)
smoother_v1(data$df$x,data$df$y,useSVD=TRUE, subsample=TRUE)

smoother_v2(data$df,useSVD=FALSE, subsample=TRUE)
smoother_v2(data$df,useSVD=TRUE, subsample=FALSE)
smoother_v2(data$df,useSVD=TRUE, subsample=TRUE)


smooth.spline(data$df$x,data$df$y,all.knots=TRUE)

p_data +
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y)$y),
            color="yellow")+
  geom_line(aes(y=data$real),
            color="green")+
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y,useSVD = TRUE)$y),
            color="orange")+
  geom_line(aes(y=smooth.spline(data$df$x,data$df$y,all.knots=TRUE)$y),
            color="red")

p_data +
  geom_line(aes(y=smoother_v1(data$df$x,data$df$y,useSVD = TRUE, subsample=TRUE)$y),
            color="darkred")

mb<-microbenchmark(smoother_v1(data$df$x,data$df$y),
                   smoother_v2(data,useSVD=FALSE, subsample = TRUE),
                   smoother_v2(data,useSVD=TRUE, subsample = TRUE),
                   smoother_v2(data,useSVD=TRUE,subsample = FALSE),
                   smooth.spline(data$df$x,data$df$y),
                   times=20)

autoplot(mb)
