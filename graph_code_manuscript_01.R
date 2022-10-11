####Code used to generate graphs for MacPherson, R.A. et al., "Genetic and Genomic Analyses of Drosophila melanogaster Models of Rare Human Diseases"
####Relevant GEO Accession number: GSE213763

### Fly Models of CSS/NCBRS/CdLS
################## Sleep / Activity Plots #######################


setwd("path\\to\\directory")

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

### Locomotor Activity
gene <- read.csv("Activity.csv", header=TRUE)


bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("Bap111", "brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))
gene1<-gene1[which(gene1$Gene != "Bap111"),]


act<-ggplot(data=gene1, aes(x=Gene, y=(value), fill=Sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs( y="Total Activity\n(Counts)")+
  scale_y_discrete(limit = c(0,1000, 2000, 3000, 4000, 5000))+
  expand_limits(y=c(0,5100))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", 
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), 
                          plot.title = element_text(hjust=0.5,size=30), 
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

act


###Night Sleep
gene <- read.csv(file= "DSleep.csv", header=TRUE)

bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("Bap111", "brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))
###
gene1<-gene1[which(gene1$Gene != "Bap111"),]


ns<-ggplot(data=gene1, aes(x=Gene, y=(mean_sleep_per_ind), fill=Sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(y="Proportion of Time\nAsleep, Night")+
  ylim(0,1.05)+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)
ns


###Day Sleep
gene <- read.csv(file= "NSleep.csv", header=TRUE)

bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("Bap111", "brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))
gene1<-gene1[which(gene1$Gene != "Bap111"),]


ds<-ggplot(data=gene1, aes(x=Gene, y=(mean_sleep_per_ind), fill=Sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(y="Proportion of Time\nAsleep, Day")+
  ylim(0,1.05)+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

ds


### Sleep Bout Count, Night
gene <- read.csv("SBC_N.csv", header=TRUE)

bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("Bap111", "brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))
###
gene1<-gene1[which(gene1$Gene != "Bap111"),]


SBC_N<-ggplot(data=gene1, aes(x=Gene, y=bout_count, fill=Sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(y="Sleep Bout Count,\nNight")+
  scale_y_discrete(limit = c(0, 5, 10, 15, 20, 25))+
  expand_limits(y=c(0,26.5))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

SBC_N


plot_grid(act, ns, SBC_N, labels=c("A", "B", "C"), align = "hv", nrow=3)



########################## startle ############################
setwd("path\\to\\directory")

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gene <- read.table("raw.data.txt", header=TRUE)


bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("Bap111", "brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))
gene1<-gene1[which(gene1$Gene != "Bap111"),]

### Startle Plot
startle<-ggplot(data=gene1, aes(x=Gene, y=Time, fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(y="Time (s)")+
  scale_y_discrete(limit = c(0, 10, 20, 30, 40, 50))+
  expand_limits(y=c(0,53))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position="none",
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)
startle

### Tapping

setwd("G:\\My Drive\\Research\\CS_BB\\startle\\tapping")

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gene <- read.csv("number_tapping_raw_data.csv", header=TRUE)


bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
Sex<-gene$sex
my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene$Gene = factor(gene$Gene, 
                   levels = c("brm", "osa", "Snr1",
                              "SMC1", "SMC3", "vtd", "Control"))
###

### Percent tapping
tap<-ggplot(data=gene, aes(x=Gene, y=(percent), fill=Sex)) + 
  geom_col(position = position_dodge(width=0.7), width=0.7, color='black', size=0.5)+
  labs(y="Percent Tapping")+
  scale_y_discrete(limit = c(0,10, 20, 30))+
  expand_limits(y=c(0,31))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position="none",
                          axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)
tap


plot_grid(startle, tap, labels = c("A", "B"), align = "hv", nrow = 1)



################# brain morphology #####################


### gross abnormalities
setwd("path\\to\\directory")


library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
gene <- read.csv("norm_brain_lengths_averaged.csv", header=TRUE)
gene<-as.data.frame(gene)


bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene


my_palette = c(brewer.pal(5, "Set1")[c(4,3)])

alpha<-gene1[gene1$lobe=="alpha",]
alpha$gene = factor(alpha$gene, 
                    levels = c("brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))

beta<-gene1[gene1$lobe=="beta",]
beta$gene = factor(beta$gene, 
                   levels = c("brm", "osa", "Snr1",
                              "SMC1", "SMC3", "vtd", "Control"))

gene2 <- read.csv("norm_brain_lengths_ellipsoid_2021_07_16.csv", header=TRUE)
gene2<-as.data.frame(gene2)

gene3<-gene2


verellip<-gene3[gene3$lobe=="verellip",]
verellip$gene = factor(verellip$gene, 
                       levels = c("brm", "osa", "Snr1",
                                  "SMC1", "SMC3", "vtd", "Control"))

horellip<-gene3[gene3$lobe=="horellip",]
horellip$gene = factor(horellip$gene, 
                       levels = c("brm", "osa", "Snr1",
                                  "SMC1", "SMC3", "vtd", "Control"))



### Plots

alobe<-ggplot(data=alpha, aes(x=gene, y=avg_norm_length, fill=sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(x="Gene", y="Normalized Length,\nAlpha Lobe")+
  #scale_y_discrete(limit = c(0.2, 0.4, 0.6, 0.8))+
  expand_limits(y=c(0.2, 0.8))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

blobe<-ggplot(data=beta, aes(x=gene, y=avg_norm_length, fill=sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(x="Gene", y="Normalized Length,\nBeta Lobe")+
  #scale_y_discrete(limit = c(0.2, 0.4, 0.6, 0.8))+
  expand_limits(y=c(0.2, 0.8))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

hellip<-ggplot(data=horellip, aes(x=gene, y=norm_length, fill=sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(x="Gene", y="Normalized Horizontal Length,\nEllipsoid Body")+
  #scale_y_continuous()+
  expand_limits(y=c(0, 0.6))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"), 
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)


vellip<-ggplot(data=verellip, aes(x=gene, y=norm_length, fill=sex)) + 
  geom_boxplot(outlier.size=0.25)+
  labs(x="Gene", y="Normalized Vertical Length,\nEllipsoid Body")+
  #scale_y_discrete(limit = c(0.2, 0.4, 0.6, 0.8))+
  expand_limits(y=c(0, 0.6))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

plot_grid(alobe, blobe, vellip, hellip, nrow=2, ncol=2, labels = c("A", 'B', 'C', 'D'), vjust = 0.9)



################## abnormal brains #######################

setwd("path\\to\\directory")


library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gene <- read.csv("Brain_abnormal_graphing.csv", header=TRUE)
gene<-as.data.frame(gene)

bold.italic<-element_text(face="bold.italic", size=12, colour = "black") 
gene1<-gene


my_palette = c(brewer.pal(5, "Set1")[c(4,3)])
gene1$Gene = factor(gene1$Gene, 
                    levels = c("brm", "osa", "Snr1",
                               "SMC1", "SMC3", "vtd", "Control"))

alpha<-gene1[gene1$Trait=="narrow alpha lobe head/impaired outgrowth",]
beta<-gene1[gene1$Trait=="beta lobe midline defects",]
total<-gene1[gene1$Trait=="Percent Abnormal",]

### brain morphology plots %
a<-ggplot(data=alpha, aes(x=Gene, y=Percent, fill=Sex)) + 
  geom_col(position = "dodge", color='black')+
  labs(x="Gene", y="Percent of Brains")+
  scale_y_discrete(limit = c(0, 20, 40, 60,80, 100))+
  expand_limits(y=c(0,101))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

b<-ggplot(data=beta, aes(x=Gene, y=Percent, fill=Sex)) + 
  geom_col(position = "dodge", color='black')+
  labs(x="Gene", y="Percent of Brains")+
  scale_y_discrete(limit = c(0, 20, 40, 60,80, 100))+
  expand_limits(y=c(0,101))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)

t<-ggplot(data=total, aes(x=Gene, y=Percent, fill=Sex)) + 
  geom_col(position = "dodge", color='black')+
  labs(x="Gene", y="Percent Abnormal Brains")+
  scale_y_discrete(limit = c(0, 20, 40, 60,80, 100))+
  expand_limits(y=c(0,101))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = element_text(size=12, color = "black"),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = element_text(size=12),axis.text.x = bold.italic)


plot_grid(a, b, t, nrow=1, ncol=3, labels = c("A", 'B', 'C'))



########### Venn Diagram Code ############

#Import Lists of DEG

setwd("path\\to\\directory_brm")
brm <- read.csv(file= "01_pvalues_FDR_9657brm.csv", header=TRUE)
setwd("path\\to\\directory_osa")
osa <- read.csv(file= "01_pvalues_FDR_9657osa.csv", header=TRUE)
setwd("path\\to\\directory_Snr1")
Snr1 <- read.csv(file= "01_pvalues_FDR_9657Snr1.csv", header=TRUE)

setwd("path\\to\\directory_smc1")
SMC1 <- read.csv(file= "01_pvalues_FDR_9657SMC1.csv", header=TRUE)
setwd("path\\to\\directory_smc3")
SMC3 <- read.csv(file= "01_pvalues_FDR_9657SMC3.csv", header=TRUE)
setwd("path\\to\\directory_vtd")
vtd <- read.csv(file= "01_pvalues_FDR_9657vtd.csv", header=TRUE)

setwd("path\\to\\directory_css")
css <- read.csv(file= "01_pvalues_FDR_9657css.csv", header=TRUE)

setwd("path\\to\\directory_cdls")
cdls <- read.csv(file= "01_pvalues_FDR_9657cdls.csv", header=TRUE)



#Subset by FDR value
#Females
brm_F_5E2 <- subset(brm, brm$FDR_bysex_female<0.05, select = Flybase_id)
brm_F_5E3 <- subset(brm, brm$FDR_bysex_female<0.005, select = Flybase_id)
brm_F_5E4 <- subset(brm, brm$FDR_bysex_female<0.0005, select = Flybase_id)

osa_F_5E2 <- subset(osa, osa$FDR_bysex_female<0.05, select = Flybase_id)
osa_F_5E3 <- subset(osa, osa$FDR_bysex_female<0.005, select = Flybase_id)
osa_F_5E4 <- subset(osa, osa$FDR_bysex_female<0.0005, select = Flybase_id)

Snr1_F_5E2 <- subset(Snr1, Snr1$FDR_bysex_female<0.05, select = Flybase_id)
Snr1_F_5E3 <- subset(Snr1, Snr1$FDR_bysex_female<0.005, select = Flybase_id)
Snr1_F_5E4 <- subset(Snr1, Snr1$FDR_bysex_female<0.0005, select = Flybase_id)

css_F_5E2 <- subset(css, css$FDR_bysex_female<0.05, select = Flybase_id)
css_F_5E3 <- subset(css, css$FDR_bysex_female<0.005, select = Flybase_id)
css_F_5E4 <- subset(css, css$FDR_bysex_female<0.0005, select = Flybase_id)



#Males

brm_M_5E2 <- subset(brm, brm$FDR_bysex_male<0.05, select = Flybase_id)
brm_M_5E3 <- subset(brm, brm$FDR_bysex_male<0.005, select = Flybase_id)
brm_M_5E4 <- subset(brm, brm$FDR_bysex_male<0.0005, select = Flybase_id)

osa_M_5E2 <- subset(osa, osa$FDR_bysex_male<0.05, select = Flybase_id)
osa_M_5E3 <- subset(osa, osa$FDR_bysex_male<0.005, select = Flybase_id)
osa_M_5E4 <- subset(osa, osa$FDR_bysex_male<0.0005, select = Flybase_id)

Snr1_M_5E2 <- subset(Snr1, Snr1$FDR_bysex_male<0.05, select = Flybase_id)
Snr1_M_5E3 <- subset(Snr1, Snr1$FDR_bysex_male<0.005, select = Flybase_id)
Snr1_M_5E4 <- subset(Snr1, Snr1$FDR_bysex_male<0.0005, select = Flybase_id)


css_M_5E2 <- subset(css, css$FDR_bysex_male<0.05, select = Flybase_id)
css_M_5E3 <- subset(css, css$FDR_bysex_male<0.005, select = Flybase_id)
css_M_5E4 <- subset(css, css$FDR_bysex_male<0.0005, select = Flybase_id)



#Subset by FDR value

#Females
SMC1_F_5E2 <- subset(SMC1, SMC1$FDR_bysex_female<0.05, select = Flybase_id)
SMC1_F_5E3 <- subset(SMC1, SMC1$FDR_bysex_female<0.005, select = Flybase_id)
SMC1_F_5E4 <- subset(SMC1, SMC1$FDR_bysex_female<0.0005, select = Flybase_id)

SMC3_F_5E2 <- subset(SMC3, SMC3$FDR_bysex_female<0.05, select = Flybase_id)
SMC3_F_5E3 <- subset(SMC3, SMC3$FDR_bysex_female<0.005, select = Flybase_id)
SMC3_F_5E4 <- subset(SMC3, SMC3$FDR_bysex_female<0.0005, select = Flybase_id)

vtd_F_5E2 <- subset(vtd, vtd$FDR_bysex_female<0.05, select = Flybase_id)
vtd_F_5E3 <- subset(vtd, vtd$FDR_bysex_female<0.005, select = Flybase_id)
vtd_F_5E4 <- subset(vtd, vtd$FDR_bysex_female<0.0005, select = Flybase_id)

cdls_F_5E2 <- subset(cdls, cdls$FDR_bysex_female<0.05, select = Flybase_id)
cdls_F_5E3 <- subset(cdls, cdls$FDR_bysex_female<0.005, select = Flybase_id)
cdls_F_5E4 <- subset(cdls, cdls$FDR_bysex_female<0.0005, select = Flybase_id)


#Males

SMC1_M_5E2 <- subset(SMC1, SMC1$FDR_bysex_male<0.05, select = Flybase_id)
SMC1_M_5E3 <- subset(SMC1, SMC1$FDR_bysex_male<0.005, select = Flybase_id)
SMC1_M_5E4 <- subset(SMC1, SMC1$FDR_bysex_male<0.0005, select = Flybase_id)

SMC3_M_5E2 <- subset(SMC3, SMC3$FDR_bysex_male<0.05, select = Flybase_id)
SMC3_M_5E3 <- subset(SMC3, SMC3$FDR_bysex_male<0.005, select = Flybase_id)
SMC3_M_5E4 <- subset(SMC3, SMC3$FDR_bysex_male<0.0005, select = Flybase_id)

vtd_M_5E2 <- subset(vtd, vtd$FDR_bysex_male<0.05, select = Flybase_id)
vtd_M_5E3 <- subset(vtd, vtd$FDR_bysex_male<0.005, select = Flybase_id)
vtd_M_5E4 <- subset(vtd, vtd$FDR_bysex_male<0.0005, select = Flybase_id)

cdls_M_5E2 <- subset(cdls, cdls$FDR_bysex_male<0.05, select = Flybase_id)
cdls_M_5E3 <- subset(cdls, cdls$FDR_bysex_male<0.005, select = Flybase_id)
cdls_M_5E4 <- subset(cdls, cdls$FDR_bysex_male<0.0005, select = Flybase_id)



```

#Count number of intersections between datasets
```{r}

#CSS/NCRBS_Female_FDR<0.05
length(brm_F_5E2$Flybase_id)
length(osa_F_5E2$Flybase_id)
length(Snr1_F_5E2$Flybase_id)
length(intersect(brm_F_5E2$Flybase_id,osa_F_5E2$Flybase_id))
length(intersect(brm_F_5E2$Flybase_id,Snr1_F_5E2$Flybase_id))
length(intersect(osa_F_5E2$Flybase_id,Snr1_F_5E2$Flybase_id))
length(intersect(intersect(brm_F_5E2$Flybase_id, osa_F_5E2$Flybase_id), Snr1_F_5E2$Flybase_id))

#CSS/NCRBS_Male_FDR<0.05
length(brm_M_5E2$Flybase_id)
length(osa_M_5E2$Flybase_id)
length(Snr1_M_5E2$Flybase_id)
length(intersect(brm_M_5E2$Flybase_id,osa_M_5E2$Flybase_id))
length(intersect(brm_M_5E2$Flybase_id,Snr1_M_5E2$Flybase_id))
length(intersect(osa_M_5E2$Flybase_id,Snr1_M_5E2$Flybase_id))
length(intersect(intersect(brm_M_5E2$Flybase_id, osa_M_5E2$Flybase_id), Snr1_M_5E2$Flybase_id))


#CDLS/Female_FDR<0.05
length(SMC1_F_5E2$Flybase_id)
length(SMC3_F_5E2$Flybase_id)
length(vtd_F_5E2$Flybase_id)
length(intersect(SMC1_F_5E2$Flybase_id,SMC3_F_5E2$Flybase_id))
length(intersect(SMC1_F_5E2$Flybase_id,vtd_F_5E2$Flybase_id))
length(intersect(SMC3_F_5E2$Flybase_id,vtd_F_5E2$Flybase_id))
length(intersect(intersect(SMC1_F_5E2$Flybase_id, SMC3_F_5E2$Flybase_id), vtd_F_5E2$Flybase_id))


#CDLS/Male_FDR<0.05
length(SMC1_M_5E2$Flybase_id)
length(SMC3_M_5E2$Flybase_id)
length(vtd_M_5E2$Flybase_id)
length(intersect(SMC1_M_5E2$Flybase_id,SMC3_M_5E2$Flybase_id))
length(intersect(SMC1_M_5E2$Flybase_id,vtd_M_5E2$Flybase_id))
length(intersect(SMC3_M_5E2$Flybase_id,vtd_M_5E2$Flybase_id))
length(intersect(intersect(SMC1_M_5E2$Flybase_id, SMC3_M_5E2$Flybase_id), vtd_M_5E2$Flybase_id))

```
#Making Plots
library(RBGL)
library(graph)
library(devtools)
library(gridtext)
#install_github("js229/Vennerable")
library(Vennerable)
library(grid)
library(RColorBrewer)
library(gridExtra)



CSS_F_E2 <- Venn(SetNames=c("brm","osa","Snr1"),Weight = c(`100` = 202, `010` = 524, `001`=2328, `110`=69, `101`=156, `011`=386, `111`= 156))
gridExtra::grid.arrange(grid::grid.grabExpr(plot(CSS_F_E2, doWeights = TRUE, type = "circles", show=list(Faces=FALSE))),top=richtext_grob("CSS/NCBRS Females",gp=gpar(fontsize=24,font=1)))

CSS_M_E2 <- Venn(SetNames=c("brm","osa","Snr1"),Weight = c(`100` = 638, `010` = 248, `001`=1269, `110`=499, `101`=1274, `011`=249, `111`= 584))
gridExtra::grid.arrange(grid::grid.grabExpr(plot(CSS_M_E2, doWeights = TRUE, type = "circles", show=list(Faces=FALSE))),top=richtext_grob("CSS/NCBRS Males",gp=gpar(fontsize=24,font=1)))


CDLS_F_E2 <- Venn(SetNames=c("SMC1","SMC3","vtd"),Weight = c(`100` = 883, `010` = 942, `001`=207, `110`=1271, `101`=113, `011`=225, `111`= 273))
gridExtra::grid.arrange(grid::grid.grabExpr(plot(CDLS_F_E2, doWeights = TRUE, type = "circles", show=list(Faces=FALSE))),top=richtext_grob("CDLS/NCBRS Females",gp=gpar(fontsize=24,font=1)))

CDLS_M_E2 <- Venn(SetNames=c("SMC1","SMC3","vtd"),Weight = c(`100` = 1201, `010` = 210, `001`=452, `110`=155, `101`=382, `011`=139, `111`= 657))
gridExtra::grid.arrange(grid::grid.grabExpr(plot(CDLS_M_E2, doWeights = TRUE, type = "circles", show=list(Faces=FALSE))),top=richtext_grob("CDLS/NCBRS Males",gp=gpar(fontsize=24,font=1)))




css_cdls<- Venn(SetNames=c("CSS-males","CSS-females","CdLS-males", "CdLS-females"),Weight = c(`1000` = 1376, `0100` = 1136, `0010`=218, `0001`=874, 
                                                                                              `1100`=573, `1010`=1803, `1001`=231,`0110`=43, `0101`=1620, `0011`=108,
                                                                                              `1110`= 372, `1101`= 251, `1011`= 178, `0111`= 29,
                                                                                              `1111`= 156))
gridExtra::grid.arrange(grid::grid.grabExpr(plot(css_cdls, doWeights = TRUE, type = "Chow-Ruskey", show=list(Faces=FALSE))),top=richtext_grob("css vs cdls",gp=gpar(fontsize=24,font=1)))


############# co reg genes ###################


setwd("path\\to\\directory")

library(readr)
gene <- read_csv("coregs.csv")


bold.italic<-element_text(face="bold.italic", size=12, colour = "black", angle=45, hjust=1,) 
gene1<-gene


creg_star<-ggplot(data=gene1[which(gene1$Assay == "coreg_startle"),], aes(x=Gene, y=Mean, fill=Sex)) + 
  geom_col(position = position_dodge(width=0.7), width=0.7, color='black', size=0.5)+
  geom_errorbar(aes(ymin=Mean-SED, ymax=Mean+SED), width=0.2, position=position_dodge(0.7))+
  labs(y="Average Time (s)\n(Experimental - Control)")+
  #scale_y_discrete(limit = c())+
  geom_hline(yintercept = 0, color = "black")+
  expand_limits(y=c(-10,10))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = 'none',axis.text = element_text( size=12, color = "black"), 
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=5)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = element_text(size=12),axis.text.x = bold.italic)


creg_tap<-ggplot(data=gene1[which(gene1$Assay == "coreg_tapping"),], aes(x=Gene, y=Mean, fill=Sex)) + 
  geom_col(position = position_dodge(width=0.7), width=0.7, color='black', size=0.5)+
  labs(y="Percent Tapping\n(Experimental - Control)")+
  scale_y_discrete(limit = c(0, 10, 20, 30, 40))+
  expand_limits(y=c(0,41))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = 'none',axis.text = element_text( size=12, color = "black"), 
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=5)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = element_text(size=12),axis.text.x = bold.italic)


creg_ns<-ggplot(data=gene1[which(gene1$Assay == "coreg_ns"),], aes(x=Gene, y=Mean, fill=Sex)) + 
  geom_col(position = position_dodge(width=0.7), width=0.7, color='black', size=0.5)+
  geom_errorbar(aes(ymin=Mean-SED, ymax=Mean+SED), width=0.2, position=position_dodge(0.7))+
  labs(y="Proportion of Time\nAsleep, Night\n(Experimental - Control)")+
  #scale_y_discrete(limit = c())+
  geom_hline(yintercept = 0, color = "black")+
  expand_limits(y=c(-0.1,0.05))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = 'none',axis.text = element_text( size=12, color = "black"), axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=5)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = element_text(size=12),axis.text.x = bold.italic)

creg_act<-ggplot(data=gene1[which(gene1$Assay == "coreg_act"),], aes(x=Gene, y=Mean, fill=Sex)) + 
  geom_col(position = position_dodge(width=0.7), width=0.7, color='black', size=0.5)+
  geom_errorbar(aes(ymin=Mean-SED, ymax=Mean+SED), width=0.2, position=position_dodge(0.7))+
  labs(y="Total Activity\n(counts)\n(Experimental - Control)")+
  #scale_y_discrete(limit = c())+
  geom_hline(yintercept = 0, color = "black")+
    expand_limits(y=c(0,400))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = 'none',axis.text = element_text( size=12, color = "black"), axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=5)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = element_text(size=12),axis.text.x = bold.italic)


plot_grid(creg_star, creg_tap, creg_act, creg_ns, labels = c("A", "B", "C", "D"), nrow=2, ncol=2, align = "hv", hjust = -0.5, vjust=01)




############### reaction norm plots #######################

#graph code for the double RNAi lines for the CSS/NCBRS/CdLS project
#reaction norm plots
#last updated by Rebecca MacPherson 2022-04-18 with p values

################ data prep ###################

setwd("path\\to\\directory")


library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gene <- read.csv("summary_stats.csv", header=TRUE)

gene1<-gene



gene1$snr1_gene = factor(gene1$snr1_gene, 
                         levels = c("0","1"))
gene1$alt_gene = factor(gene1$alt_gene, 
                        levels = c("0","1"))


snr1_labs<-c("Snr1+", "Snr1-")
CG5877_labs<-c("CG5877+", "CG5877-")
odc1_labs<-c("odc1+", "odc1-")

my_palette_m = c(brewer.pal(4, "Paired")[c(4,3)])
my_palette_f = c(brewer.pal(10, "Paired")[c(10,9)])


################ startle ###################


###startle, male, for 55690 line, with standard error
s_cg5877_m<-ggplot(data=(filter(gene1, Test %in% "startle" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "M")),
                   aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Time (s)")+
  scale_x_discrete(labels = snr1_labs)+
  #scale_y_discrete(limit = c(0, 10, 20, 30, 40, 50))+
  #expand_limits(y=c(0,51))+
  ylim(20, 50)+ #scales the y axis to start at 30
  expand_limits(y=c(20,51))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=50, label="italic(p)< 0.0001", parse = TRUE)




###startle, female, for 55690 line
s_cg5877_f<-ggplot(data=(filter(gene1, Test %in% "startle" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "F")),
                   aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Time (s)")+
  scale_x_discrete(labels = snr1_labs)+
  ylim(20, 50)+ #scales the y axis to start at 30
  expand_limits(y=c(20,51))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=50, label="italic(p)== 0.0681", parse = TRUE)



###startle, male, for 64498 line, with standard error
s_odc1_m<-ggplot(data=(filter(gene1, Test %in% "startle" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "M")),
                 aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Time (s)")+
  scale_x_discrete(labels = snr1_labs)+
  ylim(20, 50)+ #scales the y axis to start at 30
  expand_limits(y=c(20,51))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=50, label="italic(p)< 0.0001", parse = TRUE)




###startle, female, for 64498 line
s_odc1_f<-ggplot(data=(filter(gene1, Test %in% "startle" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "F")),
                 aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Time (s)")+
  scale_x_discrete(labels = snr1_labs)+
  ylim(20, 50)+ #scales the y axis to start at 30
  expand_limits(y=c(20,51))+
  #ylim(30, 50)+ #scales the y axis to start at 30
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=50, label="italic(p)< 0.0001", parse = TRUE)



################ night sleep ###################

###ns, male, for 55690 line, with standard error
ns_cg5877_m<-ggplot(data=(filter(gene1, Test %in% "ns" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "M")),
                    aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Proportion of Time\nAsleep, Night")+
  scale_x_discrete(labels = snr1_labs)+
  #scale_y_discrete(limit = c(0, 0.25, 0.50, 0.75, 1.00))+
  expand_limits(y=c(0.75,1.01))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "right",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=1.0, label="italic(p)== 0.0086", parse = TRUE)




###ns, female, for 55690 line
ns_cg5877_f<-ggplot(data=(filter(gene1, Test %in% "ns" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "F")),
                    aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Proportion of Time\nAsleep, Night")+
  scale_x_discrete(labels = snr1_labs)+
  expand_limits(y=c(0.75,1.01))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "right",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=1.0, label="italic(p)== 0.0051", parse = TRUE)



###ns, male, for 64498 line, with standard error
ns_odc1_m<-ggplot(data=(filter(gene1, Test %in% "ns" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "M")),
                  aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Proportion of Time\nAsleep, Night")+
  scale_x_discrete(labels = snr1_labs)+
  expand_limits(y=c(0.75,1.01))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "right",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=1, label="italic(p)== 0.0001", parse = TRUE)




###ns, female, for 64498 line
ns_odc1_f<-ggplot(data=(filter(gene1, Test %in% "ns" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "F")),
                  aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Proportion of Time\nAsleep, Night")+
  scale_x_discrete(labels = snr1_labs)+
  expand_limits(y=c(0.75,1.01))+
  #ylim(30, 50)+ #scales the y axis to start at 30
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "right",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=1, label="italic(p)== 0.0038", parse = TRUE)





############### activity ###################

###act, male, for 55690 line, with standard error
act_cg5877_m<-ggplot(data=(filter(gene1, Test %in% "act" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "M")),
                     aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Total Activity\n(Counts)")+
  scale_x_discrete(labels = snr1_labs)+
  scale_y_discrete(limit = c(0, 200, 400, 600, 800))+
  expand_limits(y=c(0,850))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=800, label="italic(p)== 0.0033", parse = TRUE)




###act, female, for 55690 line
act_cg5877_f<-ggplot(data=(filter(gene1, Test %in% "act" & Line %in% c("55690","55690-32372", "32372", "tripsc") & Sex %in% "F")),
                     aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Total Activity\nCounts")+
  scale_x_discrete(labels = snr1_labs)+
  scale_y_discrete(limit = c(0, 200, 400, 600, 800))+
  expand_limits(y=c(0,850))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=CG5877_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=800, label="italic(p)< 0.0001", parse = TRUE)



###act, male, for 64498 line, with standard error
act_odc1_m<-ggplot(data=(filter(gene1, Test %in% "act" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "M")),
                   aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Total Activity\nCounts")+
  scale_x_discrete(labels = snr1_labs)+
  scale_y_discrete(limit = c(0, 200, 400, 600, 800))+
  expand_limits(y=c(0,850))+
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_m)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=800, label="italic(p)< 0.0001", parse = TRUE)




###act, female, for 64498 line
act_odc1_f<-ggplot(data=(filter(gene1, Test %in% "act" & Line %in% c("64498","64498-32372", "32372", "tripsc") & Sex %in% "F")),
                   aes(x=snr1_gene, y=Mean, color = alt_gene, group = alt_gene))+ 
  geom_line(size=1.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-Std.Error, ymax=Mean+Std.Error),  width = 0.05)+
  labs(y="Total Activity\nCounts")+
  scale_x_discrete(labels = snr1_labs)+
  scale_y_discrete(limit = c(0, 200, 400, 600, 800))+
  expand_limits(y=c(0,850))+
  #ylim(30, 50)+ #scales the y axis to start at 30
  scale_color_manual(name="",
                     breaks=c("0", "1"),
                     labels=odc1_labs,
                     values=my_palette_f)+
  theme_classic() +
  theme(legend.text = bold.italic, legend.position = "none",
        axis.text = element_text(size=15, color = "black"), axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=10)), axis.title = element_text(size=15),axis.text.x = bold.italic)+
  annotate("text", x=2.2, y=800, label="italic(p)== 0.3352", parse = TRUE)





m_cg<-plot_grid(s_cg5877_m,  act_cg5877_m, ns_cg5877_m,
                labels = c("A", "B", "C"),rel_widths = c(1.17,1.3,1.88),
                ncol=3, nrow=1)
f_cg<-plot_grid(s_cg5877_f,  act_cg5877_f, ns_cg5877_f,
                labels = c("D", "E", "F"),rel_widths = c(1.17,1.3,1.88),
                ncol=3, nrow=1)


plot_grid(m_cg, f_cg, align = "hv",
          ncol=1, nrow=2)




m_od<-plot_grid(s_odc1_m,  act_odc1_m, ns_odc1_m,
                labels = c("A", "B", "C"),rel_widths = c(1.17,1.3,1.78),
                ncol=3, nrow=1)
f_od<-plot_grid(s_odc1_f,  act_odc1_f, ns_odc1_f,
                labels = c("D", "E", "F"),rel_widths = c(1.17,1.3,1.78),
                ncol=3, nrow=1)


plot_grid(m_od, f_od, align = "hv",
          ncol=1, nrow=2)


