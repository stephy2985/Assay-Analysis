# Open this file with "Final LacZ Colon Analysis"
# Many of the df in this files come from that file

library(readxl)
library(lme4)
library(stats)
library(ggplot2)

#======== Data ========

colon_df <- read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                       sheet = "R MF ")

colon_df <- filter(colon_df,colon_df$PFUs>100000)

colon_df_noOut <- colon_df[-c(46,30,52,55),]

#General with outliers
saline_colon<- subset(colon_df, Treatment == 'Saline')
ENU_colon<- subset(colon_df, Treatment == 'ENU')

#General withot outliers
colon_df_noOut <- colon_df[-c(46,30,52,55),]

saline_colon_noOut <- subset(colon_df_noOut, Treatment == 'Saline')
ENU_colon_noOut<- subset(colon_df_noOut, Treatment == 'ENU')


# Df with the MF of every run of every sample
colon_df_ind <- read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                           sheet = "MF Ind C")

# Df with no outliers from residuals.(INDIVIDUAL REPLICATES)
c_nout <-  read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                      sheet = "MF Ind NoOut")

saline_c_resid <- subset(c_nout, Treatment == 'Saline')
ENU_c_resid<- subset(c_nout, Treatment == 'ENU')

# DF without Ouliers from Residuals (MEAN by sample)

c_nout_mean <- aggregate(MF~Diet+Treatment+Sample, c_nout, mean)

saline_c_resid_mean <- subset(c_nout_mean, Treatment == 'Saline')
ENU_c_resid_mean <- subset(c_nout_mean, Treatment == 'ENU')

#=======Stats general with outliers========

b<- glm(colon_df$MF~ as.factor(colon_df$Diet)*as.factor(colon_df$Treatment), data=colon_df)
aov(b)
summary(aov(b))
TukeyHSD(aov(b))  #Summary of everything

TukeyHSD(aov(colon_df$MF~(colon_df_pfulimit$Treatment)))
TukeyHSD(aov(colon_df$MF~(colon_df_pfulimit$Diet)))

#=======Stats general w/o outliers========

x <- glm(colon_df_noOut$MF~ as.factor(colon_df_noOut$Diet)*as.factor(colon_df_noOut$Treatment), data=colon_df)
aov(x)
summary(aov(x))
TukeyHSD(aov(x))

TukeyHSD(aov(colon_df_noOut$MF~(colon_df_noOut$Treatment)))
TukeyHSD(aov(colon_df_noOut$MF~(colon_df_noOut$Diet)))


#=======Stats per Diet with outliers========

c <-  glm(saline_colon$MF~ as.factor(saline_colon$Diet), data=saline_colon)
aov(c)
summary(c)
TukeyHSD(aov(c))

TukeyHSD(aov(saline_colon$MF~saline_colon$Diet))
TukeyHSD(aov(ENU_colon$MF~ENU_colon$Diet))

#=======Stats per Diet w/o outliers==========

#make subsets of data from colon w/o outliers for each diet

d <-  glm(saline_colon_noOut$MF~ (saline_colon_noOut$Diet), data=saline_colon)
aov(d)
summary(d)

TukeyHSD(aov(saline_colon_noOut$MF~saline_colon_noOut$Diet))
TukeyHSD(aov(ENU_colon_noOut$MF~ENU_colon_noOut$Diet))


#=========== Stats with quasibinomial and **********Pairwise Tests***** RESULTS USED==========

# General data
q_yesOut <- glm(colon_df$MF~ as.factor(colon_df$Diet)*as.factor(colon_df$Treatment),family = quasipoisson())
summary(q_yesOut)
aov(q_yesOut)

#Saline
pairwise.t.test(saline_colon$MF, saline_colon$Diet, p.adjust.method = 'none')  #no adjusted P-value
pairwise.t.test(saline_colon$MF, saline_colon$Diet, p.adjust.method = 'holm')
pairwise.t.test(saline_colon$MF, saline_colon$Diet, p.adjust.method = 'bonf')

# General Data with no outliers
q_noOut <- glm(colon_df_noOut$MF~ as.factor(colon_df_noOut$Diet)*as.factor(colon_df_noOut$Treatment),family = quasipoisson())
summary(q_noOut)
aov(q_noOut)
#TukeyHSD(aov(q_noOut))

#Saline
pairwise.t.test(saline_colon_noOut$MF, saline_colon_noOut$Diet, p.adjust.method = 'none')  #no adjusted P-value
pairwise.t.test(saline_colon_noOut$MF, saline_colon_noOut$Diet, p.adjust.method = 'holm')
pairwise.t.test(saline_colon_noOut$MF, saline_colon_noOut$Diet, p.adjust.method = 'bonf')
#ENU
pairwise.t.test(ENU_colon_noOut$MF, ENU_colon_noOut$Diet, p.adjust.method = 'holm')

# Data without outliers from individual residuals

q_c_nout <- glm(c_nout$MF~ as.factor(c_nout$Diet)*as.factor(c_nout$Treatment),family = quasipoisson())
summary(q_c_nout)
aov(q_c_nout)
TukeyHSD(aov(q_c_nout))

#Saline
pairwise.t.test(saline_c_resid$MF, saline_c_resid$Diet, p.adjust.method = 'none')  #no adjusted P-value
pairwise.t.test(saline_c_resid$MF, saline_c_resid$Diet, p.adjust.method = 'holm')
pairwise.t.test(saline_c_resid$MF, saline_c_resid$Diet, p.adjust.method = 'bonf')
#ENU
pairwise.t.test(ENU_c_resid$MF, ENU_c_resid$Diet, p.adjust.method = 'holm')

# Data with mean values from residuals... no outliers

q_c_nout_mean <- glm(c_nout_mean$MF~ as.factor(c_nout_mean$Diet)*as.factor(c_nout_mean$Treatment),
                     family = quasipoisson())
summary(q_c_nout_mean)
aov(q_c_nout_mean)
TukeyHSD(aov(q_c_nout_mean))

#Saline
pairwise.t.test(saline_c_resid_mean$MF, saline_c_resid_mean$Diet, p.adjust.method = 'none')  #no adjusted P-value
pairwise.t.test(saline_c_resid_mean$MF, saline_c_resid_mean$Diet, p.adjust.method = 'holm')
pairwise.t.test(saline_c_resid_mean$MF, saline_c_resid_mean$Diet, p.adjust.method = 'bonf')
#ENU
pairwise.t.test(ENU_c_resid_mean$MF, ENU_c_resid_mean$Diet, p.adjust.method = 'holm')


#======= Generalize Linear Mixed Model ========= RESIDUALS to find outliers

glmer(colon_df_noOut$MF~ (1|colon_df_noOut$Diet)+(1 |colon_df_noOut$Treatment), data = colon_df_noOut, 
      family = poisson)

residuals(glmer(colon_df_noOut$MF~ (1|colon_df_noOut$Diet)+(1 |colon_df_noOut$Treatment), 
                data = colon_df_noOut, family = poisson))

glmer(colon_df$MF~ (1|colon_df$Diet)+(1 |colon_df$Treatment), 
      data = colon_df, family = poisson)  #of the data with with outliers in 

residuals(glmer(colon_df$MF~ (1|colon_df$Diet)+(1 |colon_df$Treatment), 
                data = colon_df, family = poisson)) 

# Finding the outliers of every sample per run
glmer(colon_df_ind$MF~ (1|colon_df_ind$Diet)+(1 |colon_df_ind$Treatment), 
      data = colon_df_ind, family = poisson)

residuals(glmer(colon_df_ind$MF~ (1|colon_df_ind$Diet)+(1 |colon_df_ind$Treatment), 
                data = colon_df_ind, family = poisson))

#============ Correlation tests!!! ============

library(ggpubr)

correl_df <-read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                           sheet = "Correlation")


saline_corr<- subset(correl_df, Treatment == 'Saline')
ENU_corr<- subset(correl_df, Treatment == 'ENU')

# Correlation Per treatment 

cor.test(saline_corr$MF, saline_corr$MF_bm, method = 'pearson') #Saline correlation C/BM
cor.test(ENU_corr$MF, ENU_corr$MF_bm, method = 'pearson') #ENU Correlation C/BM

# Correlation per Diet

cor.test(saline_corr[saline_corr$Diet == "Deficient",]$MF, 
         saline_corr[saline_corr$Diet == "Deficient",]$MF_bm, method = 'pearson')

cor.test(saline_corr[saline_corr$Diet == "Control",]$MF, 
         saline_corr[saline_corr$Diet == "Control",]$MF_bm, method = 'pearson')

cor.test(saline_corr[saline_corr$Diet == "Supplemented",]$MF, 
         saline_corr[saline_corr$Diet == "Supplemented",]$MF_bm, method = 'pearson')

ggscatter(saline_corr[saline_corr$Diet == "Deficient",], 
          x = "MF", y = "MF_bm", 
          add = "reg.line", conf.int = TRUE, 
          cor.method = "pearson", xlab = "Colon MF", ylab = "BM MF")
ggscatter(saline_corr[saline_corr$Diet == "Control",], 
          x = "MF", y = "MF_bm", 
          add = "reg.line", conf.int = TRUE, 
          cor.method = "pearson", xlab = "Colon MF", ylab = "BM MF")
ggscatter(saline_corr[saline_corr$Diet == "Supplemented",], 
          x = "MF", y = "MF_bm", 
          add = "reg.line", conf.int = TRUE, 
          cor.method = "pearson", xlab = "Colon MF", ylab = "BM MF")



