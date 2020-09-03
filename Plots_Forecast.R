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

prj <- reptoRlist("./3_2018_full/salidas/3_2018_full.prj")

prj$Rec_male

# Para graficar ... ####
yrs <- prj$Rec_male[,1]
M <- 0.3
# Brms <- out1$BDoLP*0.4
# Frms <- out1$Fpbr[3]


# Rec_male
p1 <- ggplot(data = NULL, aes(x = c(yrs, 2019, 2020))) + 
  geom_line(aes(y = c(prj$Rec_male[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$Rec_male_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p1 <- p1 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Reclutas machos') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p1
ggsave(plot_rec, filename = "Rec_male.png", width=7, height=8, dpi=300)


# Rec_total
p2 <- ggplot(data = NULL, aes(x = c(yrs, 2019, 2020))) + 
  geom_line(aes(y = c(prj$Reclutas[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$Rec_tot_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p2 <- p2 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=8)) +
  theme(legend.position = 'bottom') + ylab('Reclutas') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p2
ggsave(p2, filename = "Figs/Rec_2018.png", width=8.5, height=5.5, dpi=300)


# SSB
p3 <- ggplot(data = NULL, aes(x = c(yrs, 2019, 2020))) + 
  geom_line(aes(y = c(prj$SSB[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$SSB_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p3 <- p3 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('SSB') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p3
ggsave(p3, filename = "Figs/SSB_2018.png", width=8.5, height=5.5, dpi=300)


# Ftot
p4 <- ggplot(data = NULL, aes(x = c(yrs, 2019, 2020))) + 
  geom_line(aes(y = c(prj$F_total[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$F_total_proy[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p4 <- p4 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('F total') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p4
ggsave(p4, filename = "Figs/Ftot_2018.png", width=8.5, height=5.5, dpi=300)


# Catch
p5 <- ggplot(data = NULL, aes(x = c(yrs, 2019, 2020))) + 
  geom_line(aes(y = c(prj$Y_total[,2],rep(NA,2)), colour = 'b_18', linetype = 'b_18')) +
  geom_line(aes(y = c(rep(NA,length(yrs)-1),prj$Y_total_proj[1:3,3]), colour = 'prj', linetype = 'prj'))  +
  #geom_ribbon(data=NULL, aes(ymin=c(prj$Rec_male[,4],rep(NA,2)), ymax=c(prj$Rec_male[,5],rep(NA,2)), fill = 'grey55', alpha = 0.4)) 
  
  scale_color_manual(name = '',
                     values = c('royalblue3', 'red1'),
                     limits = c('b_18', 'prj'),
                     breaks = c('b_18', 'prj')) +
  scale_linetype_manual(name = '',
                        values = c('solid', 'dotted'),
                        limits = c('b_18', 'prj'),
                        breaks = c('b_18', 'prj'))
p5 <- p5 + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=9)) +
  theme(legend.position = 'bottom') + ylab('Catch') + xlab('Años') + scale_x_continuous(breaks=round(seq(min(yrs), 2020, by = 5),1))

p5
ggsave(p5, filename = "Figs/Catch_2018.png", width=8.5, height=5.5, dpi=300)







