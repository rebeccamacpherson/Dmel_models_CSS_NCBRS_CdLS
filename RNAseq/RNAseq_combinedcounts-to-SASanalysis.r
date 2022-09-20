#RNA seq analysis
#Written by Rebecca MacPherson as part of MacPherson et al 2022 (Listed in BioRxiv at the time of submission of this work) - Pleiotropic effects of the lncRNA Uhg4 in Drosophila melanogaster
#Last Edited by Rebecca MacPherson 2022-05-21


#####Filtering and normalization of counts data prior to statistical analysis#####

setwd("path\\to\\directory")


#load libraries
library(dplyr)
library(matrixStats)
library(data.table)
library(tidyr)
library(edgeR)

combined_counts<-read.table(file="combined_counts.txt", header = TRUE)

combined_counts <- as.data.frame(combined_counts)


# Adding rownames to combined_counts and renaming as Data

Data<-combined_counts
Data<-Data[,-1]
row.names(Data)<-combined_counts[,1]

#Remove the length column from Data and call it Data_nolength
Data_nolength<-Data[,-1]


#Data filtered by Median count>2
Datamatrix<-data.matrix(Data_nolength)
Medians<-rowMedians(Datamatrix)
Medians<-data.frame(Medians)
row.names(Medians)<-row.names(Data_nolength)
DatacbindMedian<-cbind(Data_nolength,Medians)
filterMedian<-filter(DatacbindMedian,Medians>=2)
filterMedian<-subset(filterMedian,select=-c(Medians))

#Filter out genes with proportion of non zero samples greater than 0.25
#'15915' is the no of rows in filterMedian; "42" is the number of columns. 
colnotzero<-rep(0,15915)
colprop<-rep(0,15915)
for(i in 1:15915) {
  ctr<-0
  for(j in 1:42) {
    if(filterMedian[i,j] > 0) ctr=ctr+1
  }
  colnotzero[i]<-ctr
  colprop[i]<-colnotzero[i]/42.0
}
colprop<-data.frame(colprop)
filterMediancolprop<-cbind(filterMedian,colprop)
filtercolprop<-filter(filterMediancolprop,colprop>0.25)
filteredData<-subset(filtercolprop,select=-c(colprop))

#Merging based on Geneid
filteredData<-setDT(filteredData,keep.rownames = TRUE)
Geneid<-filteredData$rn
filteredData1<-cbind(filteredData,Geneid)
merged<-merge(filteredData1,cc_cg5877[c("Geneid","Length")],by="Geneid")
rownames(merged)<-merged$Geneid
merged<-subset(merged,select=-c(Geneid))


# Prep data for normalization
Datafornormalization<-relocate(merged,"Length", .before = "CSS_22_merged_L001_sorted.bam")
head(Datafornormalization)
View(Datafornormalization)
#Datafornormalization includes Geneid (rownames) and Length (1st column of dataframe)
#


#####GeTMM normalization#####
##"The raw readcount matrix (tab-delimited text file) was used, in which the first column 
#holds the geneID from Ensembl that are used as row names in data matrix (x) in R, the second column 
#of the text file (thus the first column in x) 
#holds the gene length in kb and the remaining columns contain read counts of each sample." - (Smid et al., 2018)

# calculate RPK
y<-as.data.frame(Datafornormalization)
rownames(y)<-Datafornormalization$rn
y<-y[,-1]

#genelength in kb
genelength_kb<-((y[,1])/1000)
a<-cbind(y,genelength_kb)
x<-relocate(a,"genelength_kb", .before = "CSS_22_merged_L001_sorted.bam")
x<-x[,-1]

head(x)

rpk <- (x[,2:ncol(x)]/x[,1])


# remove length col in x
x <- x[,-1]
# for normalization purposes, no grouping of samples
group <- c(rep("A",ncol(x)))


#GeTMM
rpk.norm <- DGEList(counts=rpk,group=group)
rpk.norm <- calcNormFactors(rpk.norm)
norm.counts.rpk_edger <- cpm(rpk.norm)


write.csv(norm.counts.rpk_edger,file="norm_counts_GeTMM_CSS_omnibus.csv")

#######OMNIBUS ANALYSIS################
#####Transpose for input into SAS#####
a<-as.data.frame(norm.counts.rpk_edger)


#have to retain sample information within the df, not as column or row names
cnA<-colnames(a)
y<-rbind(cnA,a)
rny<-rownames(y)
z<-cbind(rny,y)
b<-transpose(z)

colnames(b)<-rownames(z)
b<-b[-1,] #remove the first row; that info is now in column names for pivot_longer
norm_counts_transposed<-b
#
norm_counts_forSAS<-pivot_longer(norm_counts_transposed,cols=starts_with("FBgn"),names_to = "Flybase_id",values_to = "norm_counts")



#insertion of metadata information for further data analysis
sex_exp<-as.data.frame(rep(c("M","F"), times=6, each=47745))
colnames(sex_exp)<-c("sex")
sex_con<-as.data.frame(rep(c("F","M"), times=1, each = 47745))
colnames(sex_con)<-c("sex")
sex<-rbind(sex_exp, sex_con)


replicate<-rep(c("1", "2", "3"), times=14, each = 15915)
line<-rep(c("osa", "brm", "Snr1", "SMC1", "vtd", "SMC3", "control"), times=1, each = 95490)
mutant<-rep(c("knockdown", "control"), times=c(572940, 95490))

metadata_forSAS<-cbind(norm_counts_forSAS, line, sex, replicate, mutant)


write.csv(metadata_forSAS,file="norm_counts_forSAS_RNAseq_CSSomni.csv", quote=FALSE)



##### SAS analysis ###########

#Last modified by Rebecca MacPherson May 2022
#data below are in SAS notation
##########################################################################

/*Importing the data*/
FILENAME REFFILE 'path\\to\\directory.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=668429;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_css_omni_02;
ods output LSMeans=LSMeans_css_omni_02;
ods output ClassLevels= ClassLevels_css_omni_02;
ods output NObs= NObs_css_omni_02;
ods output OverallANOVA= OverallANOVA_css_omni_02;
ods output ModelANOVA= ModelANOVA_css_omni_02;
ods output SlicedANOVA= SlicedANOVA_css_omni_02;
ods output Contrasts=Contrasts_css_omni_02;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex; 
run;





##### Post SAS FDR correction #####
# Created by Rebecca MacPherson 2022-01_26
# Last edited by Rebecca MacPherson on 2022-05-17

#grab the initial ANOVA output file from the glm SAS analysis
setwd("sas\\outputs\directory")


#output file from SAS is stacked improperly and includes Type I and Type III results
# 1-  Retain just the Type III results
# 2   Create lists of p values for each gene that correspond to each model term (Line, Sex, LinexSex)

model.anova.df<-read.csv("MODELANOVA.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType =="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex' & model.anova.df$HypothesisType =="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex' & model.anova.df$HypothesisType =="3"),]


# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line, model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex", "ProbF_linexsex", "FDR_linexsex")     

#write table
setwd("path\\tp\\directory")

write.csv(id_pvalues, file="01_pvalues_FDR_allgenes.csv", quote=FALSE)

########################################################################


######## Subset analyses with 9657 genes ##########

##################subset data based on SAS run in omnibus model###########################
# filtered list only includes genes with FDR<0.05 based on the omnibus model

setwd("path\\to\\directory")

df<-read.csv(file="01_pvalues_FDR_allgenes.csv", header=TRUE)
filterL05<-as.data.frame(df[df$FDR_line<0.05, 2])
colnames(filterL05)<-c("col1")
filterLS05<-as.data.frame(df[df$FDR_linexsex<0.05, 2])
colnames(filterLS05)<-c("col1")

filterLLS05<-unique(rbind(filterL05, filterLS05))#get list of all genes with LxS OR L at FDR<0.05

setwd("path\\to\\directory")

full_nc<-read.csv(file="norm_counts_forSAS_RNAseq_CSSomni.csv", header=TRUE)

#subset the 9657 genes
filter_nc<-full_nc %>%
  filter(Flybase_id %in% filterLLS05$col1)

#subset the 9657 genes with CSS only samples
filter_css<-filter_nc %>%
  filter(line %in% c("brm", "osa", "Snr1", "control"))

#subset the 9657 genes with CSS only samples
filter_cdls<-filter_nc %>%
  filter(line %in% c("SMC1", "SMC3", "vtd", "control"))


#subset the 9657 genes with pairwise only samples
filter_cssbrm<-filter_nc %>%
  filter(line %in% c("brm", "control"))

filter_cssosa<-filter_nc %>%
  filter(line %in% c("osa", "control"))

filter_csssnr1<-filter_nc %>%
  filter(line %in% c("Snr1", "control"))

filter_cdlssmc1<-filter_nc %>%
  filter(line %in% c("SMC1", "control"))

filter_cdlssmc3<-filter_nc %>%
  filter(line %in% c("SMC3", "control"))

filter_cdlsvtd<-filter_nc %>%
  filter(line %in% c("vtd", "control"))

setwd("G:\\My Drive\\Research\\CS_BB\\RNAseq\\Analysis-merged\\SAS\\subset\\01_presas_datamanipulation")
write.csv(filter_css, file = "filter9657_cssonly_forsas.csv", quote=FALSE)
write.csv(filter_cdls, file = "filter9657_cdlsonly_forsas.csv", quote=FALSE)
write.csv(filter_nc, file = "filter9657_forsas.csv", quote=FALSE)
write.csv(filter_cssbrm, file = "filter9657_brm_forsas.csv", quote=FALSE)
write.csv(filter_cssosa, file = "filter9657_osa_forsas.csv", quote=FALSE)
write.csv(filter_csssnr1, file = "filter9657_Snr1_forsas.csv", quote=FALSE)
write.csv(filter_cdlssmc1, file = "filter9657_SMC1_forsas.csv", quote=FALSE)
write.csv(filter_cdlssmc3, file = "filter9657_SMC3_forsas.csv", quote=FALSE)
write.csv(filter_cdlsvtd, file = "filter9657_vtd_forsas.csv", quote=FALSE)



##next - run SAS 

##### SAS analysis ###########

#Last modified by Rebecca MacPherson June 2022
#data below are in SAS notation

####### css only ##########
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\directory.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=231760;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_css_01;
ods output LSMeans=LSMeans_filter9657_css_01;
ods output ClassLevels= ClassLevels_filter9657_css_01;
ods output NObs= NObs_filter9657_css_01;
ods output OverallANOVA= OverallANOVA_filter9657_css_01;
ods output ModelANOVA= ModelANOVA_filter9657_css_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_css_01;
ods output Contrasts=Contrasts_filter9657_css_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;

run;


####### cdls only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\directory.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=231760;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_cdls_01;
ods output LSMeans=LSMeans_filter9657_cdls_01;
ods output ClassLevels= ClassLevels_filter9657_cdls_01;
ods output NObs= NObs_filter9657_cdls_01;
ods output OverallANOVA= OverallANOVA_filter9657_cdls_01;
ods output ModelANOVA= ModelANOVA_filter9657_cdls_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_cdls_01;
ods output Contrasts=Contrasts_filter9657_cdls_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;

####### brm vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_brm_01;
ods output LSMeans=LSMeans_filter9657_brm_01;
ods output ClassLevels= ClassLevels_filter9657_brm_01;
ods output NObs= NObs_filter9657_brm_01;
ods output OverallANOVA= OverallANOVA_filter9657_brm_01;
ods output ModelANOVA= ModelANOVA_filter9657_brm_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_brm_01;
ods output Contrasts=Contrasts_filter9657_brm_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;


###### osa vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';

PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_osa_01;
ods output LSMeans=LSMeans_filter9657_osa_01;
ods output ClassLevels= ClassLevels_filter9657_osa_01;
ods output NObs= NObs_filter9657_osa_01;
ods output OverallANOVA= OverallANOVA_filter9657_osa_01;
ods output ModelANOVA= ModelANOVA_filter9657_osa_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_osa_01;
ods output Contrasts=Contrasts_filter9657_osa_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;


###### RNA3 (Snr1) vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_Snr1_01;
ods output LSMeans=LSMeans_filter9657_Snr1_01;
ods output ClassLevels= ClassLevels_filter9657_Snr1_01;
ods output NObs= NObs_filter9657_Snr1_01;
ods output OverallANOVA= OverallANOVA_filter9657_Snr1_01;
ods output ModelANOVA= ModelANOVA_filter9657_Snr1_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_Snr1_01;
ods output Contrasts=Contrasts_filter9657_Snr1_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;


###### RNA4 (SMC1) vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_SMC1_01;
ods output LSMeans=LSMeans_filter9657_SMC1_01;
ods output ClassLevels= ClassLevels_filter9657_SMC1_01;
ods output NObs= NObs_filter9657_SMC1_01;
ods output OverallANOVA= OverallANOVA_filter9657_SMC1_01;
ods output ModelANOVA= ModelANOVA_filter9657_SMC1_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_SMC1_01;
ods output Contrasts=Contrasts_filter9657_SMC1_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;

###### RNA5 (SMC3) vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_SMC3_01;
ods output LSMeans=LSMeans_filter9657_SMC3_01;
ods output ClassLevels= ClassLevels_filter9657_SMC3_01;
ods output NObs= NObs_filter9657_SMC3_01;
ods output OverallANOVA= OverallANOVA_filter9657_SMC3_01;
ods output ModelANOVA= ModelANOVA_filter9657_SMC3_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_SMC3_01;
ods output Contrasts=Contrasts_filter9657_SMC3_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;


###### RNA6 (vtd) vs control only #################
/*Importing the data*/
  FILENAME REFFILE 'path\\to\\file.csv';


PROC IMPORT DATAFILE=REFFILE
DBMS=CSV
REPLACE
OUT=WORK.norm_counts_forSAS;
GETNAMES=YES;
GUESSINGROWS=115880;
RUN;

/* sort by gene id*/
  PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY FlyBase_id;

/* format output tables*/
  PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_filter9657_vtd_01;
ods output LSMeans=LSMeans_filter9657_vtd_01;
ods output ClassLevels= ClassLevels_filter9657_vtd_01;
ods output NObs= NObs_filter9657_vtd_01;
ods output OverallANOVA= OverallANOVA_filter9657_vtd_01;
ods output ModelANOVA= ModelANOVA_filter9657_vtd_01;
ods output SlicedANOVA= SlicedANOVA_filter9657_vtd_01;
ods output Contrasts=Contrasts_filter9657_vtd_01;


/*Defining and running the model for each gene (geneid)*/
  proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by Flybase_id;
lsmeans line*sex / slice= sex;
run;




##### Post SAS FDR correction #####
# Created by Rebecca MacPherson 2022-01_26
# Last edited by Rebecca MacPherson on 2022-06-01

#grab the initial ANOVA output file from the glm SAS analysis
####### css only ########
setwd("path\\to\\directory")


#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_CSS_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_CSS_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657css.csv", quote=FALSE)



###### cdls only #####
setwd("path\\to\\directory")


#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_CDLS_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_CDLS_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657cdls.csv", quote=FALSE)


###### brm only #####
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_BRM_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_BRM_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)


#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("G:\\My Drive\\Research\\CS_BB\\RNAseq\\Analysis-merged\\SAS\\subset\\03_post_sas_analysis\\filter9657_brm")

write.csv(id_pvalues, file="01_pvalues_FDR_9657brm.csv", quote=FALSE)

###### osa only ##### 
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_OSA_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_OSA_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")

write.csv(id_pvalues, file="01_pvalues_FDR_9657osa.csv", quote=FALSE)



###### Snr1 only #####
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_SNR1_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_SNR1_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657snr1.csv", quote=FALSE)



###### SMC1 only #####
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_SMC1_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_SMC1_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657smc1.csv", quote=FALSE)




###### SMC3 only #####
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_SMC3_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_SMC3_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657smc3.csv", quote=FALSE)



###### vtd only #####
setwd("path\\to\\directory")

#output file from SAS contains TYPE I and Type III
model.anova.df<-read.csv("MODELANOVA_FILTER9657_VTD_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType=="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex'& model.anova.df$HypothesisType=="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex'& model.anova.df$HypothesisType=="3"),]

sliced.anova.df<-read.csv("SLICEDANOVA_FILTER9657_VTD_01.csv")
sliced.male<-sliced.anova.df[which(sliced.anova.df$sex=='M'),]
sliced.female<-sliced.anova.df[which(sliced.anova.df$sex=='F'),]
# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

FDR_male<-p.adjust(sliced.male$ProbF, method="fdr")
FDR_female<-p.adjust(sliced.female$ProbF, method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
FBid<-as.data.frame(model.line$ï..Flybase_id)

#make sure the order is the same
table(model.line$ï..Flybase_id == model.sex$ï..Flybase_id)
table(model.line$ï..Flybase_id == model.linexsex$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.female$ï..Flybase_id)
table(model.line$ï..Flybase_id == sliced.male$ï..Flybase_id)

#combine into one df and assign proper column names
id_pvalues<-cbind(FBid, model.line$ProbF, FDR_line,
                  model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex,
                  sliced.female$ProbF, FDR_female, sliced.male$ProbF, FDR_male)
colnames(id_pvalues)<-c("Flybase_id","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex",
                        "ProbF_linexsex", "FDR_linexsex", "ProbF_bysex_female", "FDR_bysex_female",
                        "ProbF_bysex_male", "FDR_bysex_male")     

#write table
setwd("path\\to\\directory")
write.csv(id_pvalues, file="01_pvalues_FDR_9657vtd.csv", quote=FALSE)





