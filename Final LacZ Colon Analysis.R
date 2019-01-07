library(readr)
library(ggplot2)
library(tidyverse)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(readxl)

#========= Data Frames and orders ========

colon_df <- read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                       sheet = "R MF ")

colon_df <- filter(colon_df,colon_df$PFUs>100000)

saline_colon<- subset(colon_df, Treatment == 'Saline')
ENU_colon<- subset(colon_df, Treatment == 'ENU')

saline_colon_noOut <- subset(colon_df_noOut, Treatment == 'Saline')
ENU_colon_noOut<- subset(colon_df_noOut, Treatment == 'ENU')

order_diet <- c('Deficient','Control','Supplemented')  #order x-axis
colon_df$Treatment_order = factor(colon_df$Treatment, 
                                  levels = c('Saline','ENU'))

#===========Colon graph with Outliers==========

c_plt <- ggplot(colon_df, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
colon_plt<- c_plt + 
  stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5)) +
  facet_grid(Treatment_order~., scales = "free") +
  geom_point(data = colon_df) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

## Saving pdf of pic

pdf('Colon LacZ final - Outliers.pdf')
colon_plt
dev.off()

#===============Data summary with sd and mean=====

summary_dfcolon <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

colon_sd_mean <- summary_dfcolon(colon_df, varname="MF", groupnames=c("Diet","Treatment"))


#===================Plot w/out Outliers=============

colon_df_noOut <- colon_df[-c(46,30,52,55),]

c_plt1 <- ggplot(colon_df_noOut, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.shape = NA, fatten=1)
colon_plt_noOut<- c_plt1 + stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5)) +
  facet_grid(Treatment_order~., scales = "free") +
  geom_point(data = colon_df_noOut) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

pdf('Colon LacZ final NO Outliers.pdf')
colon_plt_noOut
dev.off()

colon_sd_noOut <- summary_dfcolon(colon_df_noOut, varname="MF", groupnames=c("Diet","Treatment"))

#========= Plot with INDIVIDUAL points w/out outliers (GML residuals) =========

c_nout <-  read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                      sheet = "MF Ind NoOut")

c_nout$Treat_order = factor(c_nout$Treatment,levels = c('Saline','ENU'))

c_nout_gml <- ggplot(c_nout, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
c_nout_gml<- c_nout_gml + 
  #stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5)) +
  facet_grid(Treat_order~., scales = "free") +
  geom_point(data = c_nout) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

## Saving pdf of pic

pdf('Colon LacZ No Outliers_GML Residuals.pdf')
c_nout_gml
dev.off()

#No error bars

pdf('Colon LacZ No Outliers_GML Residuals NO ERROR bars.pdf')
c_nout_gml
dev.off()

#========= Plot with Individual points NO outliers(GML residuals) + no sample 47 =========

# This graph is taking sample 47 since those points are outliers according to 
# the boxplot

c_nout$Treat_order = factor(c_nout$Treatment,levels = c('Saline','ENU'))

c_nout_gml2 <- ggplot(c_nout[-c(176,177),], aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
c_nout_gml2<- c_nout_gml2 + 
  stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5)) +
  facet_grid(Treat_order~., scales = "free") +
  geom_point(data = c_nout[-c(176,177),]) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

pdf('Colon LacZ No Outliers_GML Residuals- no 47.pdf')
c_nout_gml2
dev.off()

#========= Plot with POOLED means w/out outliers (GML residuals) ======

# DF without Ouliers from Residuals (MEAN by sample)

c_nout_mean <- aggregate(MF~Diet+Treatment+Sample, c_nout, mean)
write_csv(c_nout_mean, file("Colon Pooled means.csv"))

c_nout_mean$Treat_order = factor(c_nout_mean$Treatment,levels = c('Saline','ENU'))

#c_nout_mean[-c(47),]   # To take out sample 47

c_nout_mean_plot <- ggplot(c_nout_mean, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
c_nout_mean_plot<- c_nout_mean_plot + 
  stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5)) +
  facet_grid(Treat_order~., scales = "free") +
  geom_point(data = c_nout_mean) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

pdf('Colon LacZ No Outliers_GML Residuals MEAN values.pdf')
c_nout_mean_plot
dev.off()

#============ Tissue Comparison graphs (BM sample 54 removed)=========

tissues_df <- read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                       sheet = "C+BM+S")

# Subsetting by treatment

saline_tissues<- subset(tissues_df, Treatment == 'Saline')
ENU_tissues<- subset(tissues_df, Treatment == 'ENU')


tissues_df$Treat_order = factor(tissues_df$Treatment,levels = c('Saline','ENU'))

tissues_plot <- ggplot(tissues_df, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
tissues_plot<- tissues_plot + 
  stat_boxplot(geom ='errorbar', width = 0.2, size=0.3) +
  scale_x_discrete(name = " FA Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Tissues") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 1),
        axis.title.x = element_text(vjust = 1),
        strip.text.y = element_text(size=9, face = 'bold'),
        strip.text.x = element_text(size = 9, face = "bold"),
        plot.title = element_text(face= 'bold',vjust = 2, hjust = 0.5)) +
  facet_grid(Treat_order~Tissue, scales = "free") +
  geom_point(data = tissues_df) 
  #scale_color_manual(values = c("firebrick3", "springgreen3"))

pdf('Tissues comparison plot.pdf')
tissues_plot
dev.off()


#========= Trying the really nice graph with all tissues ========
title = "Mutant Frequency of Tissues"

theme = theme_set(theme_minimal())
theme = theme_update(legend.position="top", legend.title=element_blank(), 
                     panel.grid.major.x=element_blank())

nice = ggplot(saline_tissues, aes( y = MF, x = factor(Diet, levels = order_diet),
                                   fill= factor(Diet, levels = order_diet)))+ 
  ggtitle(title)

nice = nice + geom_boxplot(outlier.colour = NULL, color = "gray") + 
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", 
             fun.data = function(x){ return(c(y=median(x), 
                                              ymin=median(x), ymax=median(x))) })

theme = theme_update(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(), axis.title.x=element_blank())


theme = theme_update(axis.text.y=element_blank(), axis.ticks.y = element_blank(), 
                     axis.line.y = element_blank(), axis.title.y=element_blank())
#No Y Axis Label + Grey Axis Numbers
theme = theme_update(axis.line.y = element_blank(), axis.title.y=element_blank(), 
                     axis.text.y = element_text(colour="grey"), 
                     axis.ticks.y= element_line(colour="grey"))

nice = nice + facet_wrap(~Tissue, nrow = 1, scales="free")


#Same Scale Per Facet
#nice = nice + facet_grid(. ~ Tissue)




