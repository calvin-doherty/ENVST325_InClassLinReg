# read in greenhouse gas data from reservoirs
ghg <- read.csv("/cloud/project/Deemer_GHG_Data.csv")

install.packages(c("dplyr", "ggplot2", "olsrr", "PerformanceAnalytics"))
library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)

ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)

# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)

# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe

summary(mod.full)

res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)

plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 

# check full model
full.step$model

# plot AIC over time
plot(full.step )

# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")

# look at prediction with 95% confidence interval of the mean

predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")

### Almond Orchard Evapotranspiration Data

ETdat <- read.csv("/cloud/project/ETdata.csv")

unique(ETdat$crop)

library(lubridate)
library(ggplot2)
library(forecast)
library(dplyr)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")
newAlmond <- forecast(model4)
newAlmond 

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

### Homework Question 1

ghg$CO2_Transformed <- 1/(ghg$co2+1000)

CO2_model <- lm(CO2_Transformed ~ airTemp+
                  HydroV+
                  Region+
                  log.age+
                  log.precip+
                  log.DIP, data=ghg) #uses the data argument to specify dataframe

summary(CO2_model)

### Homework Question 2

## average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

## average fields for each month for pistachios
pistachio <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% # only use pistachio fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# pistachio ET time series
pistachio_ts <- ts(pistachio$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
pistachio_dec <- decompose(pistachio_ts)
# plot decomposition
plot(pistachio_dec)

## average fields for each month for fallow/idle fields
fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>% # only use fallow/idle fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# fallow ET time series
fallow_ts <- ts(fallow$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose fallow ET time series
fallow_dec <- decompose(fallow_ts)
# plot decomposition
plot(fallow_dec)

## average fields for each month for corn
corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>% # only use corn fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# corn ET time series
corn_ts <- ts(corn$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose corn ET time series
corn_dec <- decompose(corn_ts)
# plot decomposition
plot(corn_dec)

## average fields for each month for table grapes
grapes <- ETdat %>% # ET data
  filter(crop == "Grapes (Table/Raisin)") %>% # only use grape fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# grapes ET time series
grapes_ts <- ts(grapes$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose grapes ET time series
grapes_dec <- decompose(grapes_ts)
# plot decomposition
plot(grapes_dec)

# visualize the data

all_crops <- bind_rows(
  almond %>% mutate(crop = "Almond"),
  grapes %>% mutate(crop = "Grapes"),
  pistachio %>% mutate(crop = "Pistachio"),
  fallow %>% mutate(crop = "Fallow"),
  corn %>% mutate(crop = "Corn")
)

ggplot(all_crops, aes(x=ymd(date),y=ET.in))+
  geom_point(aes(col=crop))+
  geom_line(aes(col=crop))+
  labs(x="year", y="Monthly evapotranspiration (in)")

### Homework Question 3

## Fallow
fallow_y <- na.omit(fallow_ts)

model4_fallow <- arima(fallow_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

# calculate fits
AR_fit4_fallow <- fallow_y - residuals(model4_fallow)

#plot data
plot(fallow_y)

# plot fits
points(AR_fit4_fallow, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3"),
       bty="n")
newFallow <- forecast(model4_fallow)
newFallow 

#make dataframe for plotting
newFallowF <- data.frame(newFallow)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newFallowF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newFallowF$dateF[24])+  # Plotting original data
  geom_line(data = newFallowF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newFallowF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

## Pistachio

pistachio_y <- na.omit(pistachio_ts)

model4_pistachio <- arima(pistachio_y , # data 
                       order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

# calculate fits
AR_fit4_pistachio <- pistachio_y - residuals(model4_pistachio)

#plot data
plot(pistachio_y)

# plot fits
points(AR_fit4_pistachio, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3"),
       bty="n")
newPistachio <- forecast(model4_pistachio)
newPistachio 

#make dataframe for plotting
newPistachioF <- data.frame(newPistachio)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newPistachioF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = pistachio, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachio$date[1]),newPistachioF$dateF[24])+  # Plotting original data
  geom_line(data = newPistachioF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newPistachioF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")


