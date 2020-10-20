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

# 2019 2021 opcion 1 ####
#prj_1  <- reptoRlist("./1_2019_full/salidas/1_2019_full.prj")
prj_1  <- reptoRlist("./LAm_Norte/1_2019_full_1/salidas/1_2019_full_1.prj")

prj_1$Rec_male

# Para graficar ... ##
yrs <- prj_1$Rec_male[,1]
M <- 0.3
lastyr <- tail(yrs, n=1)
# Brms <- out1$BDoLP*0.4
# Frms <- out1$Fpbr[3]


# Rec_total
p1 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_1$Reclutas[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_1$Rec_tot_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p1 <- p1 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Reclutas') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p1
ggsave(p1, filename = "Figs/2019_2021/Rec19_21_1.png", width=8.5, height=5.5, dpi=300)


# SSB
p2 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_1$SSB[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_1$SSB_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p2 <- p2 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('SSB') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p2
ggsave(p2, filename = "Figs/2019_2021/SSB19_21_1.png", width=8.5, height=5.5, dpi=300)


# Ftot
p3 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_1$F_total[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_1$F_total_proy[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p3 <- p3 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('F total') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p3
ggsave(p3, filename = "Figs/2019_2021/F19_21_1.png", width=8.5, height=5.5, dpi=300)


# Catch
p4 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_1$Y_total[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_1$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  annotate('text', x=2021.7, y=850, label='969 t') +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p4 <- p4 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p4
ggsave(p4, filename = "Figs/2019_2021/Catch19_21_opt1.png", width=8.5, height=5.5, dpi=300)


# Read proyecciones opcion Rec 2
prj_2 <- reptoRlist("./LAm_Norte/1_2019_full_1/salidas/Rec2.prj")
prj_3 <- reptoRlist("./LAm_Norte/1_2019_full_1/salidas/Rec3.prj")

# Catch
p5 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_1$Y_total[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  #geom_line(aes(y = c(prj_2$Y_total[,2],rep(NA,2)), colour = 'b2_18', linetype = 'b2_18')) +
  #geom_line(aes(y = c(prj_3$Y_total[,2],rep(NA,2)), colour = 'b3_18', linetype = 'b3_18'))# +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_1$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_2$Y_total_proj[1:3,3]), colour = 'prj2', linetype = 'prj2'))  +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_3$Y_total_proj[1:3,3]), colour = 'prj3', linetype = 'prj3'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1', 'green', 'gold'),
                     limits = c('b_19', 'prj','prj2','prj3'),
                     breaks = c('b_19', 'prj','prj2','prj3')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted','dotted','dotted'),
                        limits = c('b_19', 'prj','prj2','prj3'),
                        breaks = c('b_19', 'prj','prj2','prj3'))
p5 <- p5 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p5
ggsave(p5, filename = "Figs/2019_2021/Catch_Recs19_21.png", width=8.5, height=5.5, dpi=300)


# 2019 2021 opcion 2 ####
prj_opt2  <- reptoRlist("./LAm_Norte/1_2019_full_2/salidas/1_2019_full_2.prj")

prj_opt2$Rec_male

# Para graficar ...
yrs <- prj_opt2$Rec_male[,1]
lastyr <- tail(yrs, n=1)


# Catch
p6 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_opt2$Y_total[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_opt2$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  annotate('text', x=2021, y=850, label='913 t') +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p6 <- p6 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p6
ggsave(p6, filename = "Figs/2019_2021/Catch19_21_opt2.png", width=8.5, height=5.5, dpi=300)

# 2019 2021 opcion 3 ####

prj_opt3  <- reptoRlist("./LAm_Norte/1_2019_full_3/salidas/1_2019_full_3.prj")

prj_opt3$Rec_male

# Para graficar ...
yrs <- prj_opt3$Rec_male[,1]
lastyr <- tail(yrs, n=1)


# Catch
p7 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(prj_opt3$Y_total[,2],rep(NA,2)), colour = 'b_19', linetype = 'b_19')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj_opt3$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  annotate('text', x=2021, y=850, label='922 t') +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_19', 'prj'),
                     breaks = c('b_19', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_19', 'prj'),
                        breaks = c('b_19', 'prj'))
p7 <- p7 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p7
ggsave(p7, filename = "Figs/2019_2021/Catch19_21_opt3.png", width=8.5, height=5.5, dpi=300)


# 2018 2020 opcion 1 ####
P18_opt1  <- reptoRlist("./LAm_Norte/2_2018_full_1/salidas/2_2018_full_1.prj")

P18_opt1$Rec_male

# Para graficar ...
yrs <- P18_opt1$Rec_male[,1]
lastyr <- tail(yrs, n=1)


# Catch
p8 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(P18_opt1$Y_total[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),P18_opt1$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  annotate('text', x=2020, y=750, label='897 t') +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p8 <- p8 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p8
ggsave(p8, filename = "Figs/2018_2020/Catch18_20_opt1.png", width=8.5, height=5.5, dpi=300)


# 2018 2020 opcion 2 ####
P18_opt2  <- reptoRlist("./LAm_Norte/2_2018_full_2/salidas/2_2018_full_2.prj")

P18_opt2$Rec_male

# Para graficar ...
yrs <- P18_opt2$Rec_male[,1]
lastyr <- tail(yrs, n=1)


# Catch
p9 <- ggplot(data = NULL, aes(x = c(yrs, lastyr+1, lastyr+2))) + 
  geom_line(aes(y = c(P18_opt2$Y_total[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),P18_opt2$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  annotate('text', x=2020.5, y=750, label='925 t') +
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p9 <- p9 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p9
ggsave(p9, filename = "Figs/2018_2020/Catch18_20_opt2.png", width=8.5, height=5.5, dpi=300)
