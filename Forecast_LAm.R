# Funciones y Directorios ####

rm(list=ls())

 library(stringr)
 library(ggplot2)
 library(reshape)
 library(ggpubr)

 source('~/Documents/Rwork/Functions/Funciones/functions.R')
 source('~/Documents/Rwork/Functions/Funciones/read.report.R')

 dir.1 <-'~/Documents/GitHub/LAm_Forecast/'
 dir.2<-'~/Documents/GitHub/LAm_Forecast/1_2019_full/'
 dir.3<-'~/Documents/GitHub/LAm_Forecast/2_2019_Oct/'
 dir.4<-'~/Documents/GitHub/LAm_Forecast/3_2018_full/'
 dir.5<-'~/Documents/GitHub/LAm_Forecast/4_2017_full/'
 dir.6<-'~/Documents/GitHub/LAm_Forecast/5_2016_full/'
 dir.7<-'~/Documents/GitHub/LAm_Forecast/6_2015_full/'
 dir.8<-'~/Documents/GitHub/LAm_Forecast/CBA_Proy_2018/'



 # Estimación CBA #### 
 
#Tabla 1: CBA 2019 Octubre 2019
 
setwd(dir.4) # Directorio de Modelo Base !!!! 
dir()

std1        <- read.table('LAM_18_full.std', header=T, sep="", na="NA", fill=T)
r    <- seq(0.1,0.5,0.1) # niveles de riesgo (cuantiles)                                
nr   <- length(r)                                                                                   
CBA  <- matrix(ncol=nr)
 
CBAp    <-subset(std1,name=='CBA')$value[3]
CBApstd <-subset(std1,name=='CBA')$std[3]

 for(j in 1:nr){	
   CBA[j]<-qnorm(r[j],CBAp,CBApstd)
   }
 
 tCBA <- round(cbind(CBAp,CBApstd,CBA),0)
 colnames(tCBA) <- c("mean","std",seq(10,50,10))
 tCBA
 
setwd(dir.1)
write.table(tCBA, 'CBA_2018full.txt', append = FALSE, sep = " ", dec = ".",
             row.names = TRUE, col.names = TRUE)
 
 

# Seguir con proyecciones para el informe
 
 
# Proyección CBA bajo distintos F ####

#unlink(dir.3,recursive=T)
dir.create(file.path('~/Documents/GitHub/LAm_Forecast/','CBA_Proy_2018'))
setwd(dir.4);file.copy(c('2018_full.dat', 'LAM_18_full.rep','LAM'), dir.8)
#setwd(dir.5);file.copy(c('lamnor1903.rep','lamnor1903'), dir.2)
# setwd(dir.2)
# system('mv lamnor1809s6.dat lamnor1903.dat')

setwd(dir.8)
data        <- lisread('2018_full.dat')
names(data) <- str_trim(names(data), side="right")
dat         <- data


system('cat LAM_18_full.rep')
# 0.4557

system('./LAM -ind 2018_full.dat')


# escenario 1: Ro*1
Tasa_bdpr <- c(0.3, 0.4,0.45, 0.5053)
l_npbr <- length(Tasa_bdpr)
dat$npbr   <- l_npbr
dat$Tasa_bdpr <- Tasa_bdpr
writeData(paste("LAMs",1,'.dat', sep=''), dat, append=FALSE)


# escenario 2: Ro/2
dat$pry_nR <- 0.5
dat$Tasa_bdpr <- Tasa_bdpr
writeData(paste('LAMs',2,'.dat', sep=''), dat, append=FALSE)

# escenario 3: Ro*1.5
dat$pry_nR <- 1.5
dat$Tasa_bdpr <- Tasa_bdpr
writeData(paste('LAMs',3,".dat",sep=""), dat, append=FALSE)


# Corre Escenarios ####

#setwd(dir.3)
run <- rbind("./LAM -ind  $1.dat -r","cp LAM.rep $1.rep","cp LAM.std $1.std")
cat(run, file = (can <- file("run.sh","wb",encoding="UTF-8")),sep="\n")
close(can)

n<-3
casos <- rep(NA,n)
s     <- seq(1,n,1)
for(i in 1:n){
  casos[i]<-paste("./run.sh LAMs",s[i],sep="")
}
cat(casos,file=(can<-file("casos.sh","wb",encoding="UTF-8")),sep="\n")
close(can)

system("chmod 755 run.sh")
system("bash ./casos.sh")


# Lee reportes escenarios ####

#setwd(dir.3)
yrs  <- seq(2018,2028,1)
nyr  <- length(yrs)


rep1      <-reptoRlist('LAMs1.rep')   
rep2      <-reptoRlist('LAMs2.rep')   
rep3      <-reptoRlist('LAMs3.rep')   

Brms <- rep1$BDoLP*0.4

BD_Proy <- data.frame(cbind(yrs,rep1$BD_proy[1:11,]))
colnames(BD_Proy) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')

C_Proy <- data.frame(cbind(yrs,rep1$C_proy[1:11,]))
colnames(C_Proy) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')



# Graficos ####

#setwd(dir.8)

# BD Proyectada
p1 <- ggplot(BD_Proy, aes(x = yrs)) + 
  geom_line(aes(y = F45, colour = 'F45', linetype = 'F45')) +
  geom_line(aes(y = F40, colour = 'F40', linetype = 'F40')) +
  geom_line(aes(y = Fsq, colour = 'Fsq', linetype = 'Fsq')) + 
  geom_line(aes(y = c(rep(Brms,11)), colour = 'Brms', linetype = 'Brms')) +
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3', 'black'),
                     limits = c('F45', 'F40',  'Fsq', 'Brms'),
                     breaks = c('F45', 'F40',  'Fsq', 'Brms')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash', 'dotted'),
                        limits = c('F45', 'F40', 'Fsq', 'Brms'),
                        breaks = c('F45', 'F40', 'Fsq', 'Brms')) +
  
  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Biomasa Desovante (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank())

p1
ggexport(p1, filename = "Figs_Proyec/Fig1_BD.pdf", width=8, height=4.5, dpi=300)

# op<-par(no.readonly=TRUE)
# ps.options(horizontal=F,bg="white",onefile=FALSE,paper="special")
# postscript("Figs_Proyec/Fig1.eps",height=5, width=7.5) 
# p1
dev.off()


# Captura Proyectada

p2 <- ggplot(C_Proy, aes(x = yrs)) + 
  geom_line(aes(y = F45, colour = 'F45', linetype = 'F45')) +
  geom_line(aes(y = F40, colour = 'F40', linetype = 'F40')) +
  geom_line(aes(y = Fsq, colour = 'Fsq', linetype = 'Fsq')) + 
  geom_line(aes(y= 748, colour = 'catch19', linetype = 'catch19')) +
  geom_line(aes(y= 994, colour = 'CBA20', linetype = 'CBA20'))+
  
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3', 'darkorchid2', 'coral2'),
                     limits = c('F45', 'F40',  'Fsq','catch19','CBA20'),
                     breaks = c('F45', 'F40',  'Fsq','catch19','CBA20')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash', 'dashed', 'dashed'),
                        limits = c('F45', 'F40', 'Fsq', 'catch19','CBA20'),
                        breaks = c('F45', 'F40', 'Fsq', 'catch19','CBA20')) +
  
  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Captura (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank())

p2
ggexport(p2, filename = "Rec1_catch.pdf", width=8, height=4.5, dpi=300)

dev.off()



# Captura Proyectada Rec*0.5

BD_Proy2 <- data.frame(cbind(yrs,rep2$BD_proy[1:11,]))
colnames(BD_Proy2) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')

C_Proy2 <- data.frame(cbind(yrs,rep2$C_proy[1:11,]))
colnames(C_Proy2) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')


p2 <- ggplot(C_Proy2, aes(x = yrs)) + 
  geom_line(aes(y = F45, colour = 'F45', linetype = 'F45')) +
  geom_line(aes(y = F40, colour = 'F40', linetype = 'F40')) +
  geom_line(aes(y = Fsq, colour = 'Fsq', linetype = 'Fsq')) + 
  geom_line(aes(y= 748, colour = 'catch19', linetype = 'catch19')) +
  geom_line(aes(y= 994, colour = 'CBA20', linetype = 'CBA20'))+
  
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3', 'darkorchid2', 'coral2'),
                     limits = c('F45', 'F40',  'Fsq','catch19','CBA20'),
                     breaks = c('F45', 'F40',  'Fsq','catch19','CBA20')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash', 'dashed', 'dashed'),
                        limits = c('F45', 'F40', 'Fsq', 'catch19','CBA20'),
                        breaks = c('F45', 'F40', 'Fsq', 'catch19','CBA20')) +
  
  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Captura (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank())

p2
ggexport(p2, filename = "Rec2_catch.pdf", width=8, height=4.5, dpi=300)

dev.off()



# Captura Proyectada Rec*1.5

BD_Proy3 <- data.frame(cbind(yrs,rep3$BD_proy[1:11,]))
colnames(BD_Proy3) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')

C_Proy3 <- data.frame(cbind(yrs,rep3$C_proy[1:11,]))
colnames(C_Proy3) <- c('yrs', 'F30', 'F40', 'F45', 'Fsq')

# Captura Proyectada Rec*0.5
p2 <- ggplot(C_Proy3, aes(x = yrs)) + 
  geom_line(aes(y = F45, colour = 'F45', linetype = 'F45')) +
  geom_line(aes(y = F40, colour = 'F40', linetype = 'F40')) +
  geom_line(aes(y = Fsq, colour = 'Fsq', linetype = 'Fsq')) + 
  geom_line(aes(y= 748, colour = 'catch19', linetype = 'catch19')) +
  geom_line(aes(y= 994, colour = 'CBA20', linetype = 'CBA20'))+
  
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3', 'darkorchid2', 'coral2'),
                     limits = c('F45', 'F40',  'Fsq','catch19','CBA20'),
                     breaks = c('F45', 'F40',  'Fsq','catch19','CBA20')) +
  
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash', 'dashed', 'dashed'),
                        limits = c('F45', 'F40', 'Fsq', 'catch19','CBA20'),
                        breaks = c('F45', 'F40', 'Fsq', 'catch19','CBA20')) +
  
  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Captura (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank())

p2
ggexport(p2, filename = "Rec3_catch.pdf", width=8, height=4.5, dpi=300)

dev.off()



# Escenarios Reclutamientos ####

outBD2 <- data.frame(rbind(rep1$BD_proy[1:11,],rep2$BD_proy[1:11,],rep3$BD_proy[1:11,]))
colnames(outBD2) <- c('F30', 'F40','F45','Fsq')
outYP2 <- data.frame(rbind(rep1$C_proy[1:11,],rep2$C_proy[1:11,],rep3$C_proy[1:11,]))
colnames(outYP2) <- c('F30', 'F40','F45','Fsq')


Bd2 <- as.data.frame(outBD2[,2:4]) %>% mutate (year=rep(yrs,3)) %>% 
  mutate(Rec=c(rep('Rmed',each=11),rep('Rmed / 2',each=11),rep('Rmed x 1.5',each=11))) %>% 
  melt(id.vars=c('year','Rec'))

Yp2 <- as.data.frame(outYP2[,2:4]) %>% mutate (year=rep(yrs,3)) %>% 
  mutate(Rec=c(rep('Rmed',each=11),rep('Rmed / 2',each=11), rep('Rmed x 1.5',each=11))) %>% 
  melt(id.vars=c('year','Rec'))




figBD  <- ggplot(Bd2) + aes(x = yrs) + 
  geom_line(aes(x=year, y=value, color = variable, linetype = variable)) +
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3'),
                     limits = c('F45', 'F40',  'Fsq'),
                     breaks = c('F45', 'F40',  'Fsq')) +
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash'),
                        limits = c('F45', 'F40', 'Fsq'),
                        breaks = c('F45', 'F40', 'Fsq')) +

  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Biomasa desovante proyectada (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=7)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank()) +

  facet_wrap(~Rec, dir='v') 

figBD



figYp  <- ggplot(Yp2) + aes(x = yrs) + 
  geom_line(aes(x=year, y=value, color = variable, linetype = variable)) +
  
  scale_color_manual(name = '',
                     values = c('green4', 'red', 'royalblue3'),
                     limits = c('F45', 'F40',  'Fsq'),
                     breaks = c('F45', 'F40',  'Fsq')) +
  
  scale_linetype_manual(name = '',
                        values = c('solid', 'twodash',  'twodash'),
                        limits = c('F45', 'F40', 'Fsq'),
                        breaks = c('F45', 'F40', 'Fsq')) +

  xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2029, by = 1),1)) +
  ylab('Captura proyectada (t)') + # scale_y_continuous(breaks=round(seq(min(predC), 3, by = 0.5),1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=7)) +
  theme(legend.position = 'bottom')  + theme(legend.title=element_blank()) +
  
  facet_wrap(~Rec, dir='v') 

  figYp


p3 <-  ggarrange(figBD, figYp, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

ggexport(p3, filename = "Figs_Proyec/Fig3.pdf", width=7, height=5.5, dpi=300)

