# Funciones y Directorios ####
library(rstudioapi)
library(ggplot2)
library(reshape)
library(ggpubr)
library(devtools)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

devtools::source_url("https://github.com/ale-yanez/RFunctions/blob/master/read.admb.R?raw=TRUE")
devtools::source_url("https://github.com/ale-yanez/RFunctions/blob/master/read.report.R?raw=TRUE")

prj <- reptoRlist("../1_2019_2_full/salidas/1_2019_2_full.prj")
b20_par <- reptoRlist('~/Documents/ADMwork/IFOP/2020/Lama_model/Estatus_2008/norte/CBA_Proy/LAmNs1.rep')

prj$Rec_male

# Para graficar ... ####
yrs <- prj$Rec_male[,1]
M <- 0.3
lastyr <- tail(yrs, n=1)
# Brms <- out1$BDoLP*0.4
# Frms <- out1$Fpbr[3]


# Catch
p1 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(b20_par$Desemb[1,],b20_par$C_proy[2,3]), colour = 'b20_par', linetype = 'b20_par')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b20_par', 'prj'),
                     breaks = c('b20_par', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid'),
                        limits = c('b20_par', 'prj'),
                        breaks = c('b20_par', 'prj'))
p1 <- p1 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2021, by = 5),1))

p1
ggsave(p1, filename = "Catch_1921_20par.png", width=8.5, height=5.5, dpi=300)


# Catch 5 ultimos años
p2 <- ggplot(data = NULL, aes(x = c(2015:2021))) + 
  geom_line(aes(y = c(b20_par$Desemb[1,31:36],b20_par$C_proy[2,3]), colour = 'b20_par', linetype = 'b20_par')) +
  geom_line(aes(y = c(rep(NA,4),prj$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b20_par', 'prj'),
                     breaks = c('b20_par', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid'),
                        limits = c('b20_par', 'prj'),
                        breaks = c('b20_par', 'prj'))
p2 <- p2 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(2015, 2021, by = 1),1))

p2
ggsave(p2, filename = "Catch_1921_20par_2.png", width=8.5, height=5.5, dpi=300)





# 2018 -2020

prj18_20 <- reptoRlist("../3_2018_full/salidas/3_2018_full.prj")
b19_par <- reptoRlist("../2_2019_Oct/LAM_19_oct.rep")

# Para graficar ... ####
yrs <- prj18_20$Rec_male[,1]
lastyr <- tail(yrs, n=1)


# Catch
p3 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(b19_par$Desemb[1,],b19_par$C_proy[2,3]), colour = 'b19_par', linetype = 'b19_par')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj18_20$Y_total_proj[1:3,3]), colour = 'prj18_20', linetype = 'prj18_20'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b19_par', 'prj18_20'),
                     breaks = c('b19_par', 'prj18_20')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid'),
                        limits = c('b19_par', 'prj18_20'),
                        breaks = c('b19_par', 'prj18_20'))
p3 <- p3 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p3
ggsave(p1, filename = "Catch_1820_19par.png", width=8.5, height=5.5, dpi=300)


# Catch 5 ultimos años
p4 <- ggplot(data = NULL, aes(x = c(2014:2020))) + 
  geom_line(aes(y = c(b19_par$Desemb[1,30:35],b19_par$C_proy[2,3]), colour = 'b19_par', linetype = 'b19_par')) +
  geom_line(aes(y = c(rep(NA,4),prj18_20$Y_total_proj[1:3,3]), colour = 'prj18_20', linetype = 'prj18_20'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b19_par', 'prj18_20'),
                     breaks = c('b19_par', 'prj18_20')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid'),
                        limits = c('b19_par', 'prj18_20'),
                        breaks = c('b19_par', 'prj18_20'))
p4 <- p4 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(2014, 2020, by = 1),1))

p4
ggsave(p4, filename = "Catch_1820_19par_2.png", width=8.5, height=5.5, dpi=300)




# Catch 19 full a 2020
b19_full <- reptoRlist("../1_2019_2_full/LAmN.rep")

p5 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(b19_par$Desemb[1,],b19_par$C_proy[2,3]), colour = 'b19_par', linetype = 'b19_par')) +
  geom_line(aes(y = c(b19_full$Desemb[1,],b19_full$C_proy[2,3]), colour = 'b19_full', linetype = 'b19_full')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj18_20$Y_total_proj[1:3,3]), colour = 'prj18_20', linetype = 'prj18_20'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'royalblue3', 'red1'),
                     limits = c('b19_par', 'b19_full','prj18_20'),
                     breaks = c('b19_par', 'b19_full','prj18_20')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid', 'solid'),
                        limits = c('b19_par','b19_full', 'prj18_20'),
                        breaks = c('b19_par','b19_full', 'prj18_20'))
p5 <- p5 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p5
ggsave(p5, filename = "Catch_1820_19parFull.png", width=8.5, height=5.5, dpi=300)



# Catch 7 ultimos años 19 full a 2020

p6 <- ggplot(data = NULL, aes(x = c(2014:2020))) + 
  geom_line(aes(y = c(b19_par$Desemb[1,30:35],b19_par$C_proy[2,3]), colour = 'b19_par', linetype = 'b19_par')) +
  geom_line(aes(y = c(b19_full$Desemb[1,30:35],b19_full$C_proy[2,3]), colour = 'b19_full', linetype = 'b19_full')) +
  geom_line(aes(y = c(rep(NA,4),prj18_20$Y_total_proj[1:3,3]), colour = 'prj18_20', linetype = 'prj18_20'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'royalblue3', 'red1'),
                     limits = c('b19_par', 'b19_full','prj18_20'),
                     breaks = c('b19_par', 'b19_full','prj18_20')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid', 'solid'),
                        limits = c('b19_par','b19_full', 'prj18_20'),
                        breaks = c('b19_par','b19_full', 'prj18_20'))
p6 <- p6 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(2014, 2020, by = 1),1))

p6
ggsave(p6, filename = "Catch_1820_19parFull_2.png", width=8.5, height=5.5, dpi=300)





# Catch 7 ultimos años 19 full a 2020 (19=18)
prj18_20_2 <- reptoRlist("../3_2018_2_full/salidas/3_2018_2_full.prj")

p7 <- ggplot(data = NULL, aes(x = c(2014:2020))) + 
  geom_line(aes(y = c(b19_par$Desemb[1,30:35],b19_par$C_proy[2,3]), colour = 'b19_par', linetype = 'b19_par')) +
  geom_line(aes(y = c(b19_full$Desemb[1,30:35],b19_full$C_proy[2,3]), colour = 'b19_full', linetype = 'b19_full')) +
  geom_line(aes(y = c(rep(NA,4),prj18_20$Y_total_proj[1:3,3]), colour = 'prj18_20', linetype = 'prj18_20'))  +
  geom_line(aes(y = c(rep(NA,4),prj18_20_2$Y_total_proj[1:3,3]), colour = 'prj18_20_2', linetype = 'prj18_20_2'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'royalblue3', 'red1', 'red1'),
                     limits = c('b19_par', 'b19_full','prj18_20','prj18_20_2'),
                     breaks = c('b19_par', 'b19_full','prj18_20','prj18_20_2')) +
  scale_linetype_manual(name = '',
                        values = c('dotted', 'solid', 'solid', 'dotted'),
                        limits = c('b19_par','b19_full', 'prj18_20','prj18_20_2'),
                        breaks = c('b19_par','b19_full', 'prj18_20','prj18_20_2'))
p7 <- p7 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(2014, 2020, by = 1),1))

p7
ggsave(p7, filename = "Catch_1820_19parFull_2.png", width=8.5, height=5.5, dpi=300)





