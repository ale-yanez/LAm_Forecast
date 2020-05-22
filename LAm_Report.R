# Funciones y Directorios ####
rm(list=ls())

# library(stringr)
# library(knitr)
# library(dplyr)
 library(ggplot2)
 library(reshape)
# library(ggthemes)
# library(R2admb)
# library(gridExtra)
 library(ggpubr)

# source('~/Documents/Rwork/Functions/Funciones/functions.R')
# source('~/Documents/Rwork/Functions/Funciones/read.report.R')
# source('~/Documents/Rwork/Functions/Funciones/Fn_PBRs.R')
 source('~/Documents/Rwork/Functions/Funciones/Fn_DiagramaFase.R')
 source('~/Documents/Rwork/Functions/multiplot.R')
 source('~/Documents/Rwork/Functions/read.admb.R')

dir.1<-'~/Documents/GitHub/LAm_Forecast/' # Folder donde se almacenan las 3 carpetas y el .R
dir.2<-'~/Documents/GitHub/LAm_Forecast/1_2019_full/'
dir.3<-'~/Documents/GitHub/LAm_Forecast/2_2019_Oct/'
dir.4<-'~/Documents/GitHub/LAm_Forecast/3_2018_full/'
dir.5<-'~/Documents/GitHub/LAm_Forecast/4_2017_full/'
dir.6<-'~/Documents/GitHub/LAm_Forecast/5_2016_full/'
dir.7<-'~/Documents/GitHub/LAm_Forecast/6_2015_full/'



# Modifico nombres .rep para leerlos
setwd(dir.2)
dir()
system('mv LAM.rep LAM_19_full.rep')
system('mv LAM.std LAM_19_full.std')
system('mv LAM.par LAM_19_full.par')
system('mv LAM.cor LAM_19_full.cor')

setwd(dir.3)
dir()
system('mv LAM.rep LAM_19_oct.rep')
system('mv LAM.std LAM_19_oct.std')
system('mv LAM.par LAM_19_oct.par')
system('mv LAM.cor LAM_19_oct.cor')

setwd(dir.4)
dir()
system('mv LAM.rep LAM_18_full.rep')
system('mv LAM.std LAM_18_full.std')
system('mv LAM.par LAM_18_full.par')
system('mv LAM.cor LAM_18_full.cor')

setwd(dir.5)
dir()
system('mv LAM.rep LAM_17_full.rep')
system('mv LAM.std LAM_17_full.std')
system('mv LAM.par LAM_17_full.par')
system('mv LAM.cor LAM_17_full.cor')

setwd(dir.6)
dir()
system('mv LAM.rep LAM_16_full.rep')
system('mv LAM.std LAM_16_full.std')
system('mv LAM.par LAM_16_full.par')
system('mv LAM.cor LAM_16_full.cor')

setwd(dir.7)
dir()
system('mv LAM.rep LAM_15_full.rep')
system('mv LAM.std LAM_15_full.std')
system('mv LAM.par LAM_15_full.par')
system('mv LAM.cor LAM_15_full.cor')


# Lee Reportes de 6 carpetas ####
setwd(dir.2)
out1 <- read.admb("LAM_19_full")
names(out1)
std1 <- read.table('LAM_19_full.std', header=T, sep="", na="NA", fill=T)
names(std1)
# --------------------------------------------------------------- 

setwd(dir.3)
out2 <- read.admb("LAM_19_oct")
names(out2)
std2 <- read.table('LAM_19_oct.std', header=T, sep="", na="NA", fill=T)
names(std2)
# --------------------------------------------------------------- 

setwd(dir.4)
out3 <- read.admb("LAM_18_full")
names(out3)
std3 <- read.table('LAM_18_full.std', header=T, sep="", na="NA", fill=T)
names(std3)
# --------------------------------------------------------------- 

setwd(dir.5)
out4 <- read.admb("LAM_17_full")
names(out4)
std4 <- read.table('LAM_17_full.std', header=T, sep="", na="NA", fill=T)
names(std4)
# --------------------------------------------------------------- 

setwd(dir.6)
out5 <- read.admb("LAM_16_full")
names(out5)
std5 <- read.table('LAM_16_full.std', header=T, sep="", na="NA", fill=T)
names(std5)
# --------------------------------------------------------------- 

setwd(dir.7)
out6 <- read.admb("LAM_15_full")
names(out6)
std6 <- read.table('LAM_15_full.std', header=T, sep="", na="NA", fill=T)
names(std6)
# --------------------------------------------------------------- 

# Para graficar ... ####

yrs <- out1$YRS
nyrs <- length(yrs)
tallas <- seq(10,52,1)
class(tallas)
M <- 0.3
Brms <- out1$BDoLP*0.4 # Respecto a 2019 full
Frms <- out1$Fpbr[3] # Respecto a 2019 full
B0 <- out1$BDoLP

#Observado
obsD <- out1$Desemb[1,]
obsC <- out1$CPUE[1,] ; obsC[obsC <= 0.01]   <-NA
obsS <- out1$BCRU[1,] ; obsS[obsS <= 1]  <-NA

#predichos y estimados 
predD1         <- out1$Desemb[2,]
predD2         <- out2$Desemb[2,]
predD3         <- out3$Desemb[2,]
predD4         <- out4$Desemb[2,]
predD5         <- out5$Desemb[2,]
predD6         <- out6$Desemb[2,]

Rec_est1      <- subset(std1,name=='Restim')$value
Rec_est2      <- subset(std2,name=='Restim')$value
Rec_est3      <- subset(std3,name=='Restim')$value
Rec_est4      <- subset(std4,name=='Restim')$value
Rec_est5      <- subset(std5,name=='Restim')$value
Rec_est6      <- subset(std6,name=='Restim')$value

desvRec1      <- subset(std1,name=='log_dev_Ro')$value
desvRec2      <- subset(std2,name=='log_dev_Ro')$value
desvRec3      <- subset(std3,name=='log_dev_Ro')$value
desvRec4      <- subset(std4,name=='log_dev_Ro')$value
desvRec5      <- subset(std5,name=='log_dev_Ro')$value
desvRec6      <- subset(std6,name=='log_dev_Ro')$value

BT_est1       <- subset(std1,name=='BT')$value
BT_est2       <- subset(std2,name=='BT')$value
BT_est3       <- subset(std3,name=='BT')$value
BT_est4       <- subset(std4,name=='BT')$value
BT_est5       <- subset(std5,name=='BT')$value
BT_est6       <- subset(std6,name=='BT')$value

BD_est1       <- subset(std1,name=='BD')$value
BD_est2       <- subset(std2,name=='BD')$value
BD_est3       <- subset(std3,name=='BD')$value
BD_est4       <- subset(std4,name=='BD')$value
BD_est5       <- subset(std5,name=='BD')$value
BD_est6       <- subset(std6,name=='BD')$value

F_est1        <- exp(subset(std1,name=='log_Fh')$value)
F_est2        <- exp(subset(std2,name=='log_Fh')$value)
F_est3        <- exp(subset(std3,name=='log_Fh')$value)
F_est4        <- exp(subset(std4,name=='log_Fh')$value)
F_est5        <- exp(subset(std5,name=='log_Fh')$value)
F_est6        <- exp(subset(std6,name=='log_Fh')$value)


# std 
stdpredD      <- subset(std1,name=="pred_Desemb")$std
stdpredC      <- subset(std1,name=="pred_CPUE")$std
stdpredS      <- subset(std1,name=="pred_Bcru")$std

stdRec1       <- subset(std1,name=='Restim')$std
stdRec2       <- subset(std2,name=='Restim')$std
stdRec3       <- subset(std3,name=='Restim')$std
stdRec4       <- subset(std4,name=='Restim')$std
stdRec5       <- subset(std5,name=='Restim')$std
stdRec6       <- subset(std6,name=='Restim')$std

stddesvRec1   <- subset(std1,name=='log_dev_Ro')$std
stddesvRec2   <- subset(std2,name=='log_dev_Ro')$std
stddesvRec3   <- subset(std3,name=='log_dev_Ro')$std
stddesvRec4   <- subset(std4,name=='log_dev_Ro')$std
stddesvRec5   <- subset(std5,name=='log_dev_Ro')$std
stddesvRec6   <- subset(std6,name=='log_dev_Ro')$std

stdBT1        <- subset(std1,name=='BT')$std
stdBT2        <- subset(std2,name=='BT')$std
stdBT3        <- subset(std3,name=='BT')$std
stdBT4        <- subset(std4,name=='BT')$std
stdBT5        <- subset(std5,name=='BT')$std
stdBT6        <- subset(std6,name=='BT')$std

stdBD1        <- subset(std1,name=='BD')$std
stdBD2        <- subset(std2,name=='BD')$std
stdBD3        <- subset(std3,name=='BD')$std
stdBD4        <- subset(std4,name=='BD')$std
stdBD5        <- subset(std5,name=='BD')$std
stdBD6        <- subset(std6,name=='BD')$std

stdF1         <- subset(std1,name=='log_Fh')$std
stdF2         <- subset(std2,name=='log_Fh')$std
stdF3         <- subset(std3,name=='log_Fh')$std
stdF4         <- subset(std4,name=='log_Fh')$std
stdF5         <- subset(std5,name=='log_Fh')$std
stdF6         <- subset(std6,name=='log_Fh')$std

# Confidence Intervals
cvdes         <- rep(0.1,nyrs)
cvcpue        <- rep(0.15,nyrs)
cvsurv        <- rep(0.30,nyrs)

rec1_lwr      <-Rec_est1-1.96*stdRec1
rec1_upr      <-Rec_est1+1.96*stdRec1
# rec2_lwr      <-Rec_est2-1.96*stdRec2
# rec2_upr      <-Rec_est2+1.96*stdRec2
# rec3_lwr      <-Rec_est3-1.96*stdRec3
# rec3_upr      <-Rec_est3+1.96*stdRec3

desvrec1_lwr  <- desvRec1-1.96*stddesvRec1
desvrec1_upr  <- desvRec1+1.96*stddesvRec1
# desvrec2_lwr  <- desvRec2-1.96*stddesvRec2
# desvrec2_upr  <- desvRec2+1.96*stddesvRec2
# desvrec3_lwr  <- desvRec3-1.96*stddesvRec3
# desvrec3_upr  <- desvRec3+1.96*stddesvRec3

BT1_lwr       <-BT_est1-1.96*stdBT1
BT1_upr       <-BT_est1+1.96*stdBT1
# BT2_lwr       <-BT_est2-1.96*stdBT2
# BT2_upr       <-BT_est2+1.96*stdBT2
# BT3_lwr       <-BT_est3-1.96*stdBT3
# BT3_upr       <-BT_est3+1.96*stdBT3

BD1_lwr       <-BD_est1-1.96*stdBD1
BD1_upr       <-BD_est1+1.96*stdBD1
# BD2_lwr       <-BD_est2-1.96*stdBD2
# BD2_upr       <-BD_est2+1.96*stdBD2
# BD3_lwr       <-BD_est3-1.96*stdBD3
# BD3_upr       <-BD_est3+1.96*stdBD3

F1_lwr        <-exp(log(F_est1)-1.96*stdF1)
F1_upr        <-exp(log(F_est1)+1.96*stdF1)
# F2_lwr        <-exp(log(F_est2)-1.96*stdF2)
# F2_upr        <-exp(log(F_est2)+1.96*stdF2)
# F3_lwr        <-exp(log(F_est3)-1.96*stdF3)
# F3_upr        <-exp(log(F_est3)+1.96*stdF3)

obsD95i   <- obsD*exp(-1.96*cvdes); obsD95s <- obsD*exp(1.96*cvdes)
obsC95i   <- obsC*exp(-1.96*cvcpue); obsC95s <-obsC*exp(1.96*cvcpue)
obsS95i   <- obsS*exp(-1.96*cvsurv); obsS95s <-obsS*exp(1.96*cvsurv)



#Plots Comparativos ####
setwd(dir.1)
dir()


# Reclutamiento ####

p8 <- ggplot(data = NULL, aes(x = yrs)) + 
  geom_line(aes(y = Rec_est1, colour = '2019_full', linetype = '2019_full')) +
  geom_line(aes(y = c(Rec_est2), colour = '2019_oct', linetype = '2019_oct')) +
  geom_line(aes(y = c(Rec_est3, NA), colour = '2018_full', linetype = '2018_full')) +
  geom_line(aes(y = c(Rec_est4, NA, NA), colour = '2017_full', linetype = '2017_full')) +
  geom_line(aes(y = c(Rec_est5, NA, NA, NA), colour = '2016_full', linetype = '2016_full')) +
  geom_line(aes(y = c(Rec_est6, NA, NA, NA, NA), colour = '2015_full', linetype = '2015_full')) +
  
  
  geom_ribbon(data=NULL, aes(ymin=rec1_lwr, ymax=rec1_upr), fill = 'grey37', alpha = 0.2) + 
  #geom_ribbon(data=NULL, aes(ymin=c(rec2_lwr), ymax=c(rec2_upr)), fill = 'grey70', alpha = 0.4) + 
  #geom_ribbon(data=NULL, aes(ymin=c(rec3_lwr, NA), ymax=c(rec3_upr, NA)), fill = 'grey90', alpha = 0.4) + 
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'gold1', 'green2', 'darkorchid2', 'coral2'),
                     limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full'),
                     breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'dashed', 'longdash', 'solid', 'dashed', 'longdash'),
                        limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full'),
                        breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full'))

p8 <- p8 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Reclutas x 10^6') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2019, by = 2),1))

p8


# Desvíos Reclutamiento ####

p9 <- ggplot(data = NULL, aes(x = yrs)) + 
  geom_line(aes(y = desvRec1, colour = '2019_full', linetype = '2019_full')) +
  geom_line(aes(y = c(desvRec2), colour = '2019_oct', linetype = '2019_oct')) +
  geom_line(aes(y = c(desvRec3, NA), colour = '2018_full', linetype = '2018_full')) +
  geom_line(aes(y = c(desvRec4, NA, NA), colour = '2017_full', linetype = '2017_full')) +
  geom_line(aes(y = c(desvRec5, NA, NA, NA), colour = '2016_full', linetype = '2016_full')) +
  geom_line(aes(y = c(desvRec6, NA, NA, NA, NA), colour = '2015_full', linetype = '2015_full')) +
  
  geom_ribbon(data=NULL, aes(ymin=desvrec1_lwr, ymax=desvrec1_upr), fill = 'grey37', alpha = 0.2) + 
  geom_line(aes(y = c(rep(0,35)), colour = '', linetype = '')) +
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'gold1', 'green2', 'darkorchid2', 'coral2','black'),
                     limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', ''),
                     breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', '')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'dashed', 'longdash', 'solid', 'dashed', 'longdash','dotted'),
                        limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full',''),
                        breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', ''))

p9 <- p9 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Desvíos Reclutamientos') + xlab('Años') + 
  scale_x_continuous(breaks=round(seq(min(yrs), 2019, by = 2),1))

p9

#install.packages("ggarrange", dependencies=TRUE, repos='http://cran.rstudio.com/')


plot <- ggarrange(p8, #+ rremove("x.text"), 
                  p9, 
          #labels = c("A", "B"),
          ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend = "bottom")

ggexport(plot, filename = "Figs/Reclutas.pdf")



# Biomasa Total ####

p10 <- ggplot(data = NULL, aes(x = yrs)) + 
  geom_line(aes(y = BT_est1, colour = '2019_full', linetype = '2019_full')) +
  geom_line(aes(y = c(BT_est2), colour = '2019_oct', linetype = '2019_oct')) +
  geom_line(aes(y = c(BT_est3, NA), colour = '2018_full', linetype = '2018_full')) +
  geom_line(aes(y = c(BT_est4, NA, NA), colour = '2017_full', linetype = '2017_full')) +
  geom_line(aes(y = c(BT_est5, NA, NA, NA), colour = '2016_full', linetype = '2016_full')) +
  geom_line(aes(y = c(BT_est6, NA, NA, NA, NA), colour = '2015_full', linetype = '2015_full')) +
  
  geom_ribbon(data=NULL, aes(ymin=BT1_lwr, ymax=BT1_upr), fill = 'grey37', alpha = 0.2) + 
  #geom_ribbon(data=NULL, aes(ymin=c(BT2_lwr), ymax=c(BT2_upr)),fill = 'grey70', alpha = 0.4) + 
  geom_line(aes(y = c(rep(Brms,35)), colour = '', linetype = '')) +
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'gold1', 'green2', 'darkorchid2', 'coral2', 'green'),
                     limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full','Brms'),
                     breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full','Brms')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'dashed', 'longdash','solid', 'dashed', 'longdash','longdash'),
                        limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full','Brms'),
                        breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full','Brms'))

p10 <- p10 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Biomasa Total (t)') + xlab('Años') + 
  scale_x_continuous(breaks=round(seq(min(yrs), 2019, by = 2),1))

p10


# Biomasa Desovante ####

p11 <- ggplot(data = NULL, aes(x = yrs)) + 
  geom_line(aes(y = BD_est1, colour = '2019_full', linetype = '2019_full')) +
  geom_line(aes(y = c(BD_est2), colour = '2019_oct', linetype = '2019_oct')) +
  geom_line(aes(y = c(BD_est3, NA), colour = '2018_full', linetype = '2018_full')) +
  geom_line(aes(y = c(BD_est4, NA, NA), colour = '2017_full', linetype = '2017_full')) +
  geom_line(aes(y = c(BD_est5, NA, NA, NA), colour = '2016_full', linetype = '2016_full')) +
  geom_line(aes(y = c(BD_est6, NA, NA, NA, NA), colour = '2015_full', linetype = '2015_full')) +
  
  
  geom_ribbon(data=NULL, aes(ymin=BD1_lwr, ymax=BD1_upr), fill = 'grey37', alpha = 0.2) + 
  #geom_ribbon(data=NULL, aes(ymin=c(BD2_lwr), ymax=c(BD2_upr)),fill = 'grey70', alpha = 0.4) + 
  geom_line(aes(y = c(rep(Brms,35)), colour = 'Brms', linetype = 'Brms')) +
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'gold1', 'green2', 'darkorchid2', 'coral2','chartreuse3'),
                     limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'Brms'),
                     breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'Brms')) +
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'dashed', 'longdash', 'solid', 'dashed', 'longdash','dotted'),
                        limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'Brms'),
                        breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'Brms'))

p11 <- p11 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Biomasa Desovante (t)') + xlab('Años') + 
  scale_x_continuous(breaks=round(seq(min(yrs), 2019, by = 2),1))

p11


plot <- ggarrange(p10, p11, #+ rremove("x.text"), 
                  #labels = c("A", "B"),
                  ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend = "bottom")


ggexport(plot, filename = "Figs/Biomasas.pdf")



# Mortalidad por Pesca ####

p12 <- ggplot(data = NULL, aes(x = yrs)) + 
  geom_line(aes(y = F_est1, colour = '2019_full', linetype = '2019_full')) +
  geom_line(aes(y = c(F_est2), colour = '2019_oct', linetype = '2019_oct')) +
  geom_line(aes(y = c(F_est3, NA), colour = '2018_full', linetype = '2018_full')) +
  geom_line(aes(y = c(F_est4, NA, NA), colour = '2017_full', linetype = '2017_full')) +
  geom_line(aes(y = c(F_est5, NA, NA, NA), colour = '2016_full', linetype = '2016_full')) +
  geom_line(aes(y = c(F_est6, NA, NA, NA, NA), colour = '2015_full', linetype = '2015_full')) +
  
  
  geom_ribbon(data=NULL, aes(ymin=F1_lwr, ymax=F1_upr), fill = 'grey37', alpha = 0.2) + 
  geom_line(aes(y = c(rep(M,35)), colour = 'M', linetype = 'M')) +
  geom_line(aes(y = c(rep(Frms,35)), colour = 'Frms', linetype = 'Frms')) +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'gold1', 'green2', 'darkorchid2', 'coral2', 'chartreuse3','black'),
                     limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'M', 'Frms'),
                     breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'M', 'Frms')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'dashed', 'longdash', 'solid', 'dashed', 'longdash','twodash','dotted'),
                        limits = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'M', 'Frms'),
                        breaks = c('2019_full', '2019_oct', '2018_full', '2017_full', '2016_full', '2015_full', 'M', 'Frms'))

p12 <- p12 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Mortalidad por Pesca (1/años)') + xlab('Años') + 
  scale_x_continuous(breaks=round(seq(min(yrs), 2019, by = 2),1))

p12
#ggsave(p12, file='Figs/Fig12TRY.pdf', width=7.5, height=5, dpi=300)
#dev.off()

ggexport(p12, filename = "Figs/F.pdf", width=7.5, height=5, dpi=300)

