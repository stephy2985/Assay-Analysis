# Final R Script for LacZ analysis including Stats and Graphs

library(readxl)
library(lme4)
library(stats)
library(ggplot2)
library(readr)
library(tidyverse)
library(plyr)
library(dplyr)
library(RColorBrewer)

# Import all data

# Df with the MF of every run of every sample
colon_df_ind <- read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                           sheet = "MF Ind C")

#------- Finding the outliers of every sample per run----- Residual of ind.peplicates----
glmer(colon_df_ind$MF~ (1|colon_df_ind$Diet)+(1 |colon_df_ind$Treatment), 
      data = colon_df_ind, family = poisson)

ind_residuals <- residuals(glmer(colon_df_ind$MF~ (1|colon_df_ind$Diet)+(1 |colon_df_ind$Treatment), 
                                 data = colon_df_ind, family = poisson))

#---- DF USED-----

# Df with no outliers from residuals.(INDIVIDUAL REPLICATES)
c_nout <-  read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                      sheet = "MF Ind NoOut")

saline_c_resid <- subset(c_nout, Treatment == 'Saline')
ENU_c_resid<- subset(c_nout, Treatment == 'ENU')

      # Df no outliers from ind.residuals - Also NO sample #47
c_nout_no47 <-  read_excel("~/Documents/1st Year/Main Colon LacZ PFUs ,mutant plaques 3.xlsx", 
                      sheet = "MF Ind NoOut -47")

saline_c_res.47 <- subset(c_nout_no47, Treatment == 'Saline')
ENU_c_res.47<- subset(c_nout_no47, Treatment == 'ENU')


# DF without Ouliers from Residuals (MEAN by sample)

c_nout_mean <- aggregate(MF~Diet+Treatment+Sample, c_nout, mean)
write_csv(c_nout_mean, file("Colon Pooled means.csv"))

saline_c_resid_mean <- subset(c_nout_mean, Treatment == 'Saline')
ENU_c_resid_mean <- subset(c_nout_mean, Treatment == 'ENU')

# Df no Outliers fron Ind.residuals (Mean w/no 47) 
c_nout_mean_no47 <- aggregate(MF~Diet+Treatment+Sample, c_nout_no47, mean)
write_csv(c_nout_mean_no47, file("Colon Pooled means - no 47.csv"))

saline_c_res.47_mean <- subset(c_nout_mean_no47, Treatment == 'Saline')
ENU_c_res.47_mean <- subset(c_nout_mean_no47, Treatment == 'ENU')


#----------Plot of Means (no outliers - inclusive47)----- USED ----
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




#----------Plot of Means (no outliers - NO sample 47)-----

c_nout_mean_no47$Treat_order = factor(c_nout_mean_no47$Treatment,levels = c('Saline','ENU'))

#c_nout_mean[-c(47),]   # To take out sample 47

c_nout_no47_plot <- ggplot(c_nout_mean_no47, aes(x= factor(Diet, levels = order_diet), y = MF )) + 
  geom_boxplot(fill="paleturquoise", colour="cornsilk4",
               outlier.size = 0, fatten=1)
c_nout_no47_plot<- c_nout_no47_plot + 
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
  geom_point(data = c_nout_mean_no47) +
  scale_color_manual(values = c("firebrick3", "springgreen3"))

pdf('Final Colon LacZ No Outliers Ind.Residuals no47- MEAN values.pdf')
c_nout_no47_plot
dev.off()






#------Stats for Means (no out - no sample 47)---- USED!!!----
q_c_nout.no47_mean <- glm(c_nout_mean_no47$MF~ as.factor(c_nout_mean_no47$Diet)*as.factor(c_nout_mean_no47$Treatment),
                     family = quasipoisson())
summary(q_c_nout.no47_mean)
aov(q_c_nout.no47_mean)
TukeyHSD(aov(q_c_nout.no47_mean))

#Saline
pairwise.t.test(saline_c_res.47_mean$MF, saline_c_res.47_mean$Diet, p.adjust.method = 'none')  #no adjusted P-value
pairwise.t.test(saline_c_res.47_mean$MF, saline_c_res.47_mean$Diet, p.adjust.method = 'holm')
pairwise.t.test(saline_c_res.47_mean$MF, saline_c_res.47_mean$Diet, p.adjust.method = 'bonf')
#ENU
pairwise.t.test(ENU_c_res.47_mean$MF, ENU_c_res.47_mean$Diet, p.adjust.method = 'bonf')


#------Stats for Means (no out - Inclusive sample47)-----
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



#----PLots  Saline and ENU Separate (No 47)----


# SALINE individual plot
c_sal_plot <- ggplot(saline_c_res.47_mean, aes(x= factor(Diet, levels = order_diet), y = MF, fill= Diet)) + 
  geom_boxplot(colour= "#948d88", fatten = 1, coef = NULL) + #Take the NULL= normal lines of error
  scale_fill_manual(breaks = c("Deficient", "Control", "Supplemented"), 
                     values = c("#68D6E4", "#fc9399", "#f2ae8b"))


c_sal_plot<- c_sal_plot + 
  scale_x_discrete(name = "Diets") + 
  scale_y_continuous(name = "Mutant Frequency\n(10^-5)") +
  ggtitle("Mutant Frequency of Colonic Tissue
          In saline Treatmnet") +
  theme_bw() + 
  theme(axis.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title.y = element_text(vjust = 2),
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(vjust = 2, hjust = 0.5, face="bold"),
        legend.position="none",
        axis.text = element_text(colour = "black", size = rel(1))) +
  geom_point(data = saline_c_res.47_mean)
  