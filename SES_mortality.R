####################################################################################################
####   block 0: set directory
####################################################################################################
rm(list=ls())
setwd("~/Desktop/SES_mortality_code_for_share/")
require("survey")
require("openxlsx")
require("mediation")
require("Hmisc")
set.seed(3543269)
extract_mediation_summary <- function (x) {
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) ||
                   (inherits(x$model.y, "glm") && x$model.y$family$family ==
                      "gaussian" && x$model.y$family$link == "identity") ||
                   (inherits(x$model.y, "survreg") && x$model.y$dist ==
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)",
                        "ADE (control)", "ADE (treated)", "Total Effect",
                        "Prop. Mediated (control)", "Prop. Mediated (treated)",
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""),
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}


clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge",
                "SkinBloodAge","ZhangAge",
                "PhenoAge","GrimAgeMort","GrimAge2Mort",
                "DunedinPoAm",
                "YangCell","HorvathTelo")
# clock_list <- c("GDF15Mort","B2MMort","CystatinCMort","TIMP1Mort","ADMMort","PAI1Mort","LeptinMort",
#                 "PACKYRSMort","CRPMort","logA1CMort")

# supplementary bio markers
if ("CystatinCMort" %in% clock_list){namelist <- "supplementary biomarkers"} else {namelist <- "main"}

tt_tic <- proc.time()
####################################################################################################


####################################################################################################
####    block 1: loading datasets
####################################################################################################
nhBioA <- read.csv("Raw_clean_SESmortality_01_21_2025.csv")

# Creat death age
nhBioA$deathage <- nhBioA$permth_int/12+ nhBioA$age
table(nhBioA$dead,useNA = "a")
nhBioA$undead <- NA
nhBioA$undead[nhBioA$dead==1] <- 0
nhBioA$undead[nhBioA$dead==0] <- 1
table(nhBioA$dead,nhBioA$undead,useNA = "a")

# Create binary exposures
table(paste0(nhBioA$otherHisp,nhBioA$mexican,nhBioA$latinovswhite))
nhBioA$latino <- 0
nhBioA$latino[nhBioA$latinovswhite==1] <- 1
table(paste0(nhBioA$otherHisp,nhBioA$mexican,nhBioA$latino))

nhBioA$otherHisp<-ifelse(nhBioA$ridreth1==2 ,1,0)
nhBioA$otherRace<-ifelse(nhBioA$ridreth1==5 ,1,0)
table(paste0(nhBioA$mexican,nhBioA$white,
             nhBioA$black,nhBioA$otherHisp,nhBioA$otherRace),useNA = "a")
nhBioA$otherRacevswhite <- NA
nhBioA$otherRacevswhite[nhBioA$otherRace==1] <- 1
nhBioA$otherRacevswhite[nhBioA$white==1] <- 0
table(nhBioA$otherRacevswhite,useNA = "a")

nhBioA$drinker[!is.na(nhBioA$abstainer)] <- 0
nhBioA$drinker[nhBioA$abstainer==0] <- 1
table(nhBioA$abstainer,nhBioA$drinker,useNA = "a")

nhBioA$Smoker[!is.na(nhBioA$nonSmoker)] <- 0
nhBioA$Smoker[nhBioA$nonSmoker==0] <- 1
table(nhBioA$Smoker,nhBioA$nonSmoker,useNA = "a")
table(nhBioA$nonSmoker,nhBioA$packyrs0,useNA = "a")
nhBioA$packyrslt30[nhBioA$packyrs0==1&nhBioA$nonSmoker==0&!is.na(nhBioA$packyrslt30)] <- NA 
nhBioA$packyrs30_59[nhBioA$packyrs0==1&nhBioA$nonSmoker==0&!is.na(nhBioA$packyrs0)] <- NA 
nhBioA$packyrs60[nhBioA$packyrs0==1&nhBioA$nonSmoker==0&!is.na(nhBioA$packyrs0)] <- NA 
nhBioA$packyrs0[nhBioA$packyrs0==1&nhBioA$nonSmoker==0&!is.na(nhBioA$packyrs0)] <- NA 


nhBioA$lthsvscoll <- NA
nhBioA$lthsvscoll[nhBioA$lths==1] <- 1
nhBioA$lthsvscoll[nhBioA$coll==1] <- 0
table(nhBioA$lthsvscoll,useNA = "a")
nhBioA$lths_hsvscoll <- NA
nhBioA$lths_hsvscoll[nhBioA$lths==1] <- 1
nhBioA$lths_hsvscoll[nhBioA$hs==1] <- 1
nhBioA$lths_hsvscoll[nhBioA$coll==1] <- 0
table(nhBioA$lths_hsvscoll,useNA = "a")
nhBioA$hsvscoll <- NA
nhBioA$hsvscoll[nhBioA$hs==1] <- 1
nhBioA$hsvscoll[nhBioA$coll==1] <- 0
table(nhBioA$hsvscoll,useNA = "a")
nhBioA$somecollvscoll <- NA
nhBioA$somecollvscoll[nhBioA$somecoll==1] <- 1
nhBioA$somecollvscoll[nhBioA$coll==1] <- 0
table(nhBioA$somecollvscoll,useNA = "a")
nhBioA$hsvssomecoll <- NA
nhBioA$hsvssomecoll[nhBioA$hs==1] <- 1
nhBioA$hsvssomecoll[nhBioA$somecoll==1] <- 0
table(nhBioA$hsvssomecoll,useNA = "a")

nhBioA$below1vsabove5 <- NA
nhBioA$below1vsabove5[nhBioA$pir_below1==1] <- 1
nhBioA$below1vsabove5[nhBioA$pir_above5==1] <- 0
table(nhBioA$below1vsabove5,useNA = "a")
nhBioA$below2vsabove5 <- NA
nhBioA$below2vsabove5[nhBioA$pir_below1==1] <- 1
nhBioA$below2vsabove5[nhBioA$pir_1_2==1] <- 1
nhBioA$below2vsabove5[nhBioA$pir_above5==1] <- 0
table(nhBioA$below2vsabove5,useNA = "a")
nhBioA$pir_1_2vsabove5 <- NA
nhBioA$pir_1_2vsabove5[nhBioA$pir_1_2==1] <- 1
nhBioA$pir_1_2vsabove5[nhBioA$pir_above5==1] <- 0
table(nhBioA$pir_1_2vsabove5,useNA = "a")
nhBioA$pir_2_5vsabove5 <- NA
nhBioA$pir_2_5vsabove5[nhBioA$pir_2_5==1] <- 1
nhBioA$pir_2_5vsabove5[nhBioA$pir_above5==1] <- 0
table(nhBioA$pir_2_5vsabove5,useNA = "a")

nhBioA$lowwhitevshiwhite <- NA
nhBioA$lowwhitevshiwhite[nhBioA$lowwhite==1] <- 1
nhBioA$lowwhitevshiwhite[nhBioA$hiwhite==1] <- 0
table(nhBioA$lowwhitevshiwhite,useNA = "a")
nhBioA$hibluevshiwhite <- NA
nhBioA$hibluevshiwhite[nhBioA$hiblue==1] <- 1
nhBioA$hibluevshiwhite[nhBioA$hiwhite==1] <- 0
table(nhBioA$hibluevshiwhite,useNA = "a")
nhBioA$lowbluevshiwhite <- NA
nhBioA$lowbluevshiwhite[nhBioA$lowblue==1] <- 1
nhBioA$lowbluevshiwhite[nhBioA$hiwhite==1] <- 0
table(nhBioA$lowbluevshiwhite,useNA = "a")
nhBioA$noworkvshiwhite <- NA
nhBioA$noworkvshiwhite[nhBioA$nowork==1] <- 1
nhBioA$noworkvshiwhite[nhBioA$hiwhite==1] <- 0
table(nhBioA$noworkvshiwhite,useNA = "a")
####################################################################################################
#### Load DNA methylation data
library(haven)
DNAm_Meta <- read_sas("dnmepi.sas7bdat")
names(DNAm_Meta)
head(DNAm_Meta)
analysis <- merge(nhBioA,DNAm_Meta,by.x="seqn",by.y="SEQN",all = T)

# Create indicator variable for clock data
analysis$overallSample<- 0
analysis$overallSample[!is.na(analysis$HorvathAge)] <- 1
table(analysis$overallSample,useNA = "a")
tapply(analysis$WTDN4YR,analysis$overallSample,median,na.rm=T)

# Create Z-SCORE for GrimAge2 variables
for (clock in c("GDF15Mort","B2MMort","CystatinCMort","TIMP1Mort","ADMMort","PAI1Mort","LeptinMort",
                "PACKYRSMort","CRPMort","logA1CMort")){
  #clock <- "GDF15Mort"
  print(clock)
  print(range(analysis[,clock],na.rm=T))
  tempmean <- mean(analysis[,clock],na.rm=T)
  tempsd <- sd(analysis[,clock],na.rm=T)
  analysis[,clock] <- (analysis[,clock]- tempmean)/tempsd
  print(range(analysis[,clock],na.rm=T))
}
####################################################################################################


####################################################################################################
####    block 2: multiple imputation
####################################################################################################
multiimpu_analysis <- read.csv("Impute_clean_SES_01_22_2025.csv")
multiimpu_analysis$forborn <- 0
multiimpu_analysis$forborn[multiimpu_analysis$nativity==0] <- 1
multiimpu_analysis$sedentary <- NA
multiimpu_analysis$sedentary <- ifelse(multiimpu_analysis$active==0,1,0)
multiimpu_analysis$male <- NA
multiimpu_analysis$male <- ifelse(multiimpu_analysis$female==0,1,0)
multiimpu_analysis$nonSmoker[multiimpu_analysis$packyrs0==0] <- 0
multiimpu_analysis$nonSmoker[multiimpu_analysis$packyrs0==1] <- 1
multiimpu_analysis$Smoker[multiimpu_analysis$nonSmoker==0] <- 1
multiimpu_analysis$Smoker[multiimpu_analysis$nonSmoker==1] <- 0
table(multiimpu_analysis$Smoker,multiimpu_analysis$nonSmoker,useNA = "a")
multiimpu_analysis <- merge(multiimpu_analysis,analysis[,c("seqn","WTDN4YR","sdmvpsu","sdmvstra","overallSample",clock_list,
                                                           "dead","ucod_113","deathage",
                                                           "hrtmort","canmort","accimort",
                                                           "resmort","cermort","diabmort",
                                                           "hrtmort_competing","canmort_competing","accimort_competing",
                                                           "resmort_competing","cermort_competing","diabmort_competing"
)], all=T,by="seqn")

multiimpu_analysis$undead <- NA
multiimpu_analysis$undead[multiimpu_analysis$dead==1] <- 0
multiimpu_analysis$undead[multiimpu_analysis$dead==0] <- 1
table(multiimpu_analysis$dead,multiimpu_analysis$undead,useNA = "a")

multiimpu_analysis$drinker[!is.na(multiimpu_analysis$abstainer)] <- 0
multiimpu_analysis$drinker[multiimpu_analysis$abstainer==0] <- 1
table(multiimpu_analysis$abstainer,multiimpu_analysis$drinker,useNA = "a")

multiimpu_analysis$blackvswhite <- NA
multiimpu_analysis$blackvswhite[multiimpu_analysis$black==1] <- 1
multiimpu_analysis$blackvswhite[multiimpu_analysis$white==1] <- 0
table(multiimpu_analysis$blackvswhite,useNA = "a")
multiimpu_analysis$latinovswhite <- NA
multiimpu_analysis$latinovswhite[multiimpu_analysis$mexican==1] <- 1
multiimpu_analysis$latinovswhite[multiimpu_analysis$otherHisp==1] <- 1
multiimpu_analysis$latinovswhite[multiimpu_analysis$white==1] <- 0
table(multiimpu_analysis$latinovswhite,useNA = "a")
multiimpu_analysis$latinovsblack <- NA
multiimpu_analysis$latinovsblack[multiimpu_analysis$mexican==1] <- 1
multiimpu_analysis$latinovsblack[multiimpu_analysis$otherHisp==1] <- 1
multiimpu_analysis$latinovsblack[multiimpu_analysis$black==1] <- 0
table(multiimpu_analysis$latinovsblack,useNA = "a")
multiimpu_analysis$otherRacevswhite <- NA
multiimpu_analysis$otherRacevswhite[multiimpu_analysis$otherRace==1] <- 1
multiimpu_analysis$otherRacevswhite[multiimpu_analysis$white==1] <- 0
table(multiimpu_analysis$otherRacevswhite,useNA = "a")

multiimpu_analysis$latino <- 0
multiimpu_analysis$latino[multiimpu_analysis$latinovswhite==1] <- 1
table(paste0(multiimpu_analysis$otherHisp,multiimpu_analysis$mexican,multiimpu_analysis$latino))


multiimpu_analysis$lthsvscoll <- NA
multiimpu_analysis$lthsvscoll[multiimpu_analysis$lths==1] <- 1
multiimpu_analysis$lthsvscoll[multiimpu_analysis$coll==1] <- 0
table(multiimpu_analysis$lthsvscoll,useNA = "a")
multiimpu_analysis$lths_hsvscoll <- NA
multiimpu_analysis$lths_hsvscoll[multiimpu_analysis$lths==1] <- 1
multiimpu_analysis$lths_hsvscoll[multiimpu_analysis$hs==1] <- 1
multiimpu_analysis$lths_hsvscoll[multiimpu_analysis$coll==1] <- 0
table(multiimpu_analysis$lths_hsvscoll,useNA = "a")
multiimpu_analysis$hsvscoll <- NA
multiimpu_analysis$hsvscoll[multiimpu_analysis$hs==1] <- 1
multiimpu_analysis$hsvscoll[multiimpu_analysis$coll==1] <- 0
table(multiimpu_analysis$hsvscoll,useNA = "a")
multiimpu_analysis$somecollvscoll <- NA
multiimpu_analysis$somecollvscoll[multiimpu_analysis$somecoll==1] <- 1
multiimpu_analysis$somecollvscoll[multiimpu_analysis$coll==1] <- 0
table(multiimpu_analysis$somecollvscoll,useNA = "a")
multiimpu_analysis$hsvssomecoll <- NA
multiimpu_analysis$hsvssomecoll[multiimpu_analysis$hs==1] <- 1
multiimpu_analysis$hsvssomecoll[multiimpu_analysis$somecoll==1] <- 0
table(multiimpu_analysis$hsvssomecoll,useNA = "a")

multiimpu_analysis$below1vsabove5 <- NA
multiimpu_analysis$below1vsabove5[multiimpu_analysis$pir_below1==1] <- 1
multiimpu_analysis$below1vsabove5[multiimpu_analysis$pir_above5==1] <- 0
table(multiimpu_analysis$below1vsabove5,useNA = "a")
multiimpu_analysis$below2vsabove5 <- NA
multiimpu_analysis$below2vsabove5[multiimpu_analysis$pir_below1==1] <- 1
multiimpu_analysis$below2vsabove5[multiimpu_analysis$pir_1_2==1] <- 1
multiimpu_analysis$below2vsabove5[multiimpu_analysis$pir_above5==1] <- 0
table(multiimpu_analysis$below2vsabove5,useNA = "a")
multiimpu_analysis$pir_1_2vsabove5 <- NA
multiimpu_analysis$pir_1_2vsabove5[multiimpu_analysis$pir_1_2==1] <- 1
multiimpu_analysis$pir_1_2vsabove5[multiimpu_analysis$pir_above5==1] <- 0
table(multiimpu_analysis$pir_1_2vsabove5,useNA = "a")
multiimpu_analysis$pir_2_5vsabove5 <- NA
multiimpu_analysis$pir_2_5vsabove5[multiimpu_analysis$pir_2_5==1] <- 1
multiimpu_analysis$pir_2_5vsabove5[multiimpu_analysis$pir_above5==1] <- 0
table(multiimpu_analysis$pir_2_5vsabove5,useNA = "a")

multiimpu_analysis$lowwhitevshiwhite <- NA
multiimpu_analysis$lowwhitevshiwhite[multiimpu_analysis$lowwhite==1] <- 1
multiimpu_analysis$lowwhitevshiwhite[multiimpu_analysis$hiwhite==1] <- 0
table(multiimpu_analysis$lowwhitevshiwhite,useNA = "a")
multiimpu_analysis$hibluevshiwhite <- NA
multiimpu_analysis$hibluevshiwhite[multiimpu_analysis$hiblue==1] <- 1
multiimpu_analysis$hibluevshiwhite[multiimpu_analysis$hiwhite==1] <- 0
table(multiimpu_analysis$hibluevshiwhite,useNA = "a")
multiimpu_analysis$lowbluevshiwhite <- NA
multiimpu_analysis$lowbluevshiwhite[multiimpu_analysis$lowblue==1] <- 1
multiimpu_analysis$lowbluevshiwhite[multiimpu_analysis$hiwhite==1] <- 0
table(multiimpu_analysis$lowbluevshiwhite,useNA = "a")
multiimpu_analysis$noworkvshiwhite <- NA
multiimpu_analysis$noworkvshiwhite[multiimpu_analysis$nowork==1] <- 1
multiimpu_analysis$noworkvshiwhite[multiimpu_analysis$hiwhite==1] <- 0
table(multiimpu_analysis$noworkvshiwhite,useNA = "a")
###################################################################################################


####################################################################################################
####    block 3: creating raw data
####################################################################################################
rawanalysis <- analysis[analysis$overallSample==1&!is.na(analysis$overallSample),]
rawanalysis <- rawanalysis[!rawanalysis$WTDN4YR==0,]
rawanalysis <- rawanalysis[rawanalysis$age<85,]

rawanalysis_multiimpu <- multiimpu_analysis[multiimpu_analysis$overallSample==1&!is.na(multiimpu_analysis$overallSample),]
rawanalysis_multiimpu <- rawanalysis_multiimpu[!rawanalysis_multiimpu$WTDN4YR==0,]
rawanalysis_multiimpu <- rawanalysis_multiimpu[rawanalysis_multiimpu$age<85,]
####################################################################################################


####################################################################################################
####    block 4: survey data for final estimates
####################################################################################################
names(DNAm_Meta)
analysis <- analysis[!is.na(analysis$WTDN4YR),]
svyNHE <- survey::svydesign(id = ~sdmvpsu , strata = ~sdmvstra , nest = TRUE ,
                            weights = ~WTDN4YR, data = analysis)
svyNHEanalysis <- subset(svyNHE,overallSample==1&!is.na(overallSample))
svyNHEanalysis <- subset(svyNHEanalysis,age<85)
dim(svyNHEanalysis$variables)

multiimpu_analysis <- multiimpu_analysis[!is.na(multiimpu_analysis$WTDN4YR),]
svyNHE_multiimpu <- survey::svydesign(id = ~sdmvpsu , strata = ~sdmvstra , nest = TRUE ,
                                      weights = ~WTDN4YR, data = multiimpu_analysis)
svyNHEanalysis_multiimpu <- subset(svyNHE_multiimpu,overallSample==1&!is.na(overallSample))
svyNHEanalysis_multiimpu <- subset(svyNHEanalysis_multiimpu,age<85)
dim(svyNHEanalysis_multiimpu$variables)
####################################################################################################


####################################################################################################
####    block 5: sample descriptive statistics
####################################################################################################

####################################################################################################
## block 5.1: sample descriptive for categorical variables without missing
####################################################################################################
variable_list_sample_categorical <- c("white","black","mexican","otherHisp","otherRace","",
                                      "female","male","", "nativity","forborn","",
                                      "lths","hs","somecoll","coll","",
                                      "pir_below1","pir_1_2","pir_2_5","pir_above5","",
                                      "hiwhite","lowwhite","hiblue","lowblue","nowork","",
                                      "active","sedentary","","abstainer","drinker","",
                                      "dead","hrtmort","canmort")

## descriptive tables
i=1; descriptive_categorical <- NULL
for (var in variable_list_sample_categorical){
  explore_dataset <- get("rawanalysis") ## get dataset
  explore_surveydata <- get("svyNHEanalysis")
  
  if (var==""){
    temp_2 <- as.data.frame(matrix(nrow=1,ncol = 8))
  } else{
    print(var)
    
    temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
    temp_1$Var1 <- as.character(temp_1$Var1)
    temp_1$Var1[temp_1$Var1==1] <- "Yes"
    temp_1$Var1[temp_1$Var1==0] <- "No"
    
    temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])
    
    temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))
    
    temp_1 <- temp_1[order(temp_1$Var1),]
    temp_1 <- temp_1[!is.na(temp_1$Var1),]
    
    svytt <- as.data.frame(survey::svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
    names(svytt)[2] <- "total_se"
    svyper <- as.data.frame(survey::svymean(design=explore_surveydata,make.formula(var),na.rm = T))
    names(svyper)[2] <- "mean_se"
    
    temp_1 <- cbind(temp_1,svytt,svyper)
    temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
  }
  
  if (!is.null(names(descriptive_categorical))){
    names(temp_2) <- names(descriptive_categorical)}
  
  descriptive_categorical <- rbind(descriptive_categorical,temp_2)
  
  i=i+1
}
descriptive_categorical$Var1 <- NULL
names(descriptive_categorical) <- c("variable",paste0("sample_",c("count","percentage")),
                                    paste0("survey_",c("total","total_se","percentage","percentage_se")))
row.names(descriptive_categorical) <- NULL
wb <- createWorkbook()
addWorksheet(wb, "Table 1_nonmissing_raw")
writeData(wb,x=descriptive_categorical,sheet = "Table 1_nonmissing_raw",
          rowNames = F)
####################################################################################################
############## multiple imputed data
####################################################################################################
i=1; descriptive_categorical <- NULL
for (var in variable_list_sample_categorical){
  explore_dataset <- get("rawanalysis_multiimpu") ## get dataset
  explore_surveydata <- get("svyNHEanalysis_multiimpu")
  
  if (var==""){
    temp_2 <- as.data.frame(matrix(nrow=1,ncol = 8))
  } else if (!var %in% names(explore_dataset)){
    next
  } else{
    print(var)
    
    temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
    temp_1$Var1 <- as.character(temp_1$Var1)
    temp_1$Var1[temp_1$Var1==1] <- "Yes"
    temp_1$Var1[temp_1$Var1==0] <- "No"
    
    temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])
    
    temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))
    
    temp_1 <- temp_1[order(temp_1$Var1),]
    temp_1 <- temp_1[!is.na(temp_1$Var1),]
    
    svytt <- as.data.frame(survey::svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
    names(svytt)[2] <- "total_se"
    svyper <- as.data.frame(survey::svymean(design=explore_surveydata,make.formula(var),na.rm = T))
    names(svyper)[2] <- "mean_se"
    
    temp_1 <- cbind(temp_1,svytt,svyper)
    temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
  }
  
  if (!is.null(names(descriptive_categorical))){
    names(temp_2) <- names(descriptive_categorical)}
  
  descriptive_categorical <- rbind(descriptive_categorical,temp_2)
  
  i=i+1
}
descriptive_categorical$Var1 <- NULL
names(descriptive_categorical) <- c("variable",paste0("sample_",c("count","percentage")),
                                    paste0("survey_",c("total","total_se","percentage","percentage_se")))
row.names(descriptive_categorical) <- NULL
addWorksheet(wb, "Table 1_nonmissing_multiimpu")
writeData(wb,x=descriptive_categorical,sheet = "Table 1_nonmissing_multiimpu",
          rowNames = F)
####################################################################################################

####################################################################################################
## block 5.2: sample descriptive for continuous variables
####################################################################################################
variable_list_sample_continuous <- c("age","agesq",
                                     "lbxlypct","lbxmopct","lbxnepct","lbxeopct","lbxbapct",
                                     "hei","packyrs","drinkvol","waist2thigh","BMI",
                                     "lbdtcsi","HDL","LDL","Glucose","CRP")
## descriptive tables
i=1; descriptive_continuous <- NULL
for (var in variable_list_sample_continuous){
  explore_dataset <- get("rawanalysis") ## get dataset
  explore_surveydata <- get("svyNHEanalysis")
  
  print(var)
  
  if (length(table(is.na(explore_dataset[,var])))>1) {
    temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
  } else {temp_0 <- 0}
  
  temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
  temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
  temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
                     c(.25,.5,.75))
  temp_4 <- cbind(var,as.data.frame(temp_0),
                  as.data.frame(temp_1),as.data.frame(temp_2),
                  t(as.data.frame(temp_3)),table(is.na(explore_dataset[,var]))[1])
  
  svytt <- as.data.frame(survey::svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
  names(svytt)[2] <- "total_se"
  svyper <- as.data.frame(survey::svymean(design=explore_surveydata,make.formula(var),na.rm = T))
  svyvariance <- as.data.frame(survey::svyvar(design=explore_surveydata,make.formula(var),na.rm = T))
  svyiqr <- as.data.frame(survey::svyquantile(design=explore_surveydata,make.formula(var),
                                              quantiles = c(0.25,0.5,0.75), na.rm = T)[var])
  names(svyper)[2] <- "mean_se"
  
  temp_4 <- cbind(temp_4,svytt,svyper,t(svyiqr[,1]))
  
  descriptive_continuous <- rbind(descriptive_continuous,temp_4)
  i=i+1
}
names(descriptive_continuous) <- c("variable",paste0("sample_",c("missing_rate","mean","se",
                                                                 "lower_quartile","median","upper_quartile","count")),
                                   paste0("survey_",c("total","total_se","mean","mean_se","lower_quartile",
                                                      "median","upper_quartile")))

addWorksheet(wb, "Table 1_continuous_raw")
writeData(wb,x=descriptive_continuous,sheet = "Table 1_continuous_raw",
          rowNames = F)
####################################################################################################
############## multiple imputed data
####################################################################################################
i=1; descriptive_continuous <- NULL
for (var in variable_list_sample_continuous){
  explore_dataset <- get("rawanalysis_multiimpu") ## get dataset
  explore_surveydata <- get("svyNHEanalysis_multiimpu")
  
  print(var)
  
  if (length(table(is.na(explore_dataset[,var])))>1) {
    temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
  } else {temp_0 <- 0}
  
  temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
  temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
  temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
                     c(.25,.5,.75))
  temp_4 <- cbind(var,as.data.frame(temp_0),
                  as.data.frame(temp_1),as.data.frame(temp_2),
                  t(as.data.frame(temp_3)),table(is.na(explore_dataset[,var]))[1])
  
  svytt <- as.data.frame(survey::svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
  names(svytt)[2] <- "total_se"
  svyper <- as.data.frame(survey::svymean(design=explore_surveydata,make.formula(var),na.rm = T))
  svyvariance <- as.data.frame(survey::svyvar(design=explore_surveydata,make.formula(var),na.rm = T))
  svyiqr <- as.data.frame(survey::svyquantile(design=explore_surveydata,make.formula(var),
                                              quantiles = c(0.25,0.5,0.75), na.rm = T)[var])
  names(svyper)[2] <- "mean_se"
  
  temp_4 <- cbind(temp_4,svytt,svyper,t(svyiqr[,1]))
  
  descriptive_continuous <- rbind(descriptive_continuous,temp_4)
  i=i+1
}
names(descriptive_continuous) <- c("variable",paste0("sample_",c("missing_rate","mean","se",
                                                                 "lower_quartile","median","upper_quartile","count")),
                                   paste0("survey_",c("total","total_se","mean","mean_se","lower_quartile",
                                                      "median","upper_quartile")))

addWorksheet(wb, "Table 1_continuous_multiimpu")
writeData(wb,x=descriptive_continuous,sheet = "Table 1_continuous_multiimpu",
          rowNames = F)
####################################################################################################

####################################################################################################
## block 5.3: correlation matrix for clocks
####################################################################################################
tcorrelation <- rcorr(as.matrix(rawanalysis[,c(clock_list,"age","sedentary","hei","packyrs","drinker",
                                               "waist2thigh","BMI",
                                               "lbdtcsi","HDL","LDL","Glucose","CRP")]))
tcorrelation_write <- as.matrix(tcorrelation$r)
library(ggcorrplot)
pdf(paste0("Results/correlation plot",format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
    width = 10,height = 5)
ggcorrplot(tcorrelation_write, type = "lower", lab_size = 1.5, lab=T, 
           colors = c("#6D9EC1", "white", "#E46726"))+ theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1,size=5),  # Rotate x-axis labels
    axis.text.y.right = element_text(size=5),  # Right y-axis labels
    axis.ticks.y.right = element_line(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) + scale_y_discrete(position = "right")+
  guides(fill = guide_colorbar(title = "Correlation coefficients"))+
  labs(fill = "Correlation", x = NULL, y = NULL)
dev.off()

test <- rawanalysis[rawanalysis$overallSample==1,c(clock_list,"age")]
test <- data.frame(scale(test, center=TRUE, scale=TRUE))
EFAresult1 =  factanal(~ ., data=test, factors = 7, rotation = "none", 
                       na.action = na.exclude)
EFAresult1#1, 2, 6, 7, 8
# library(lavaan)
# efa.model <- '
#     efa("efa")*f1 + 
#     efa("efa")*f2 + 
#     efa("efa")*f3 + 
#     efa("efa")*f4 + 
#     efa("efa")*f5 + 
#     efa("efa")*f6 + 
#     efa("efa")*f7 + 
#     efa("efa")*f8 =~ HannumAge + HorvathAge + WeidnerAge + LinAge + VidalBraloAge + 
#     SkinBloodAge + ZhangAge + PhenoAge + GrimAgeMort+GrimAge2Mort+DunedinPoAm+YangCell+
#     HorvathTelo+age+sedentary+hei+packyrs+drinker+waist2thigh+BMI+lbdtcsi+
# #   HDL+LDL+Glucose+CRP
# '
# fit <- cfa(efa.model, data = test)
# summary(fit, standardized = TRUE)

EFA_results <- as.matrix(EFAresult1$loadings)
addWorksheet(wb, "EFA_raw")

writeData(wb,x=EFA_results,sheet = "EFA_raw",rowNames = T)
####################################################################################################
saveWorkbook(wb, file=paste0("Results/Table 1",format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
####    block 6: Main analysis 
####################################################################################################
behav_biomarker <- c("sedentary","hei","packyrs","drinker",
                     "waist2thigh","BMI",
                     "lbdtcsi","HDL","LDL","Glucose","CRP")
clock_list <- c(clock_list,behav_biomarker)



SES_tablelist <- list(c("blackvswhite","latinovswhite","otherRacevswhite"),#"otherHispvswhite","otherRacevswhite",
                      c("lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll"),
                      c("below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5"),
                      c("lowwhitevshiwhite","hibluevshiwhite",
                        "lowbluevshiwhite","noworkvshiwhite"))#

# covariates <- "+female+age+agesq+forborn+lbxlypct+lbxmopct+lbxnepct+lbxeopct+lbxbapct"
####################################################################################################
## block 6.1: Survival mediation main models
####################################################################################################
loopgroup <- c("overallSample","female","male","black","white","latino")#"female","male","black","white","latino"


for (Group in loopgroup){
  #Group <- "black"
  print(Group)
  
  if (Group %in% c("overallSample","female","male")){
    var_list <- c("blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
                  "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
                  "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
                  "lowwhitevshiwhite","hibluevshiwhite",
                  "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  } else {
    var_list <- c(#"blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
      "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
      "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
      "lowwhitevshiwhite","hibluevshiwhite",
      "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  }
  
  if (Group %in% c("female","male")){
    covariates <- "+age+agesq+forborn"
  } else {
    covariates <- "+female+age+agesq+forborn"
  }

  for (var in var_list){
    # var <- "lthsvscoll"
    
    explore_surveydata <- get("svyNHEanalysis_multiimpu")  ## get imputed survey object
    # explore_surveydata <- get("svyNHEanalysis") ## get completed cases survey object
    
    explore_surveydata$variables$test <- explore_surveydata$variables[,var]
    explore_surveydata <- subset(explore_surveydata,!is.na(test)&WTDN4YR>0)
    
    explore_surveydata$variables$Group <- explore_surveydata$variables[,Group]
    explore_surveydata <- subset(explore_surveydata,Group==1)
    #print(var)
    
    ##############  MODELS
    model_printtable <-  NULL
    
    for (clock in clock_list){
      #clock <- "Horvath"; clock <- "DNAmGrimAge2"
      
      print(clock)
      #### mediator model
      med.fit <- survey::svyglm(design=explore_surveydata,paste0(clock,"~",var,covariates),family="gaussian")
      summary(med.fit)
      med.fit$df.null
      
      #### mediator_outcome model
      med_out.fit <- survey::svysurvreg(design=explore_surveydata,
                                        as.formula(paste0("Surv(deathage,dead)~",clock,
                                                          covariates)), dist="weibull")
      summary(med_out.fit)
      
      
      #### outcome model without mediators
      #table(rawanalysis$deathage)
      # model_1_nomed <- survey::svycoxph(design=explore_surveydata,
      #                                   as.formula(paste0("Surv(deathage,dead)~",var,
      #                                                     covariates)))
      # model_1 <- survey::svycoxph(design=explore_surveydata,
      #                             as.formula(paste0("Surv(deathage,dead)~",var,"+",clock,
      #                                               covariates)))
      # summary(model_1)
      
      
      
      model_3_nomed <- survreg(as.formula(paste0("Surv(deathage,dead)~",var,
                                                 covariates)), dist="weibull",data=rawanalysis)
      
      model_3 <- survreg(as.formula(paste0("Surv(deathage,dead)~",var,"+",clock,
                                           covariates)), dist="weibull",data=rawanalysis)
      summary(model_3)
      
      
      out.fit_nomed <- survey::svysurvreg(design=explore_surveydata,
                                          as.formula(paste0("Surv(deathage,dead)~",var,
                                                            covariates)), dist="weibull")
      tryCatch({      
        out.fit <- survey::svysurvreg(design=explore_surveydata,
                                      as.formula(paste0("Surv(deathage,dead)~",var,"+",clock,
                                                        covariates)), dist="weibull")
        summary(out.fit)
        
        #### causal mediation models
        med.out <- mediate(med.fit,out.fit,treat=var,mediator=clock)
        summary(med.out)
        
        ###################################################################################################
        #### Create table
        temp <- as.data.frame(extract_mediation_summary(summary(med.out)))[c(5,8:10),]
        temp$label <- rownames(temp)
        temp <- rbind(names(temp),temp)
        temp <- cbind(temp,matrix(NA,nrow = nrow(temp),ncol=1))
        temp$sample <- med.out$nobs
        names(temp) <- 1:7 },error=function(e){
          cat("Error",conditionMessage(e),"\n");temp <- as.data.frame(matrix(data = NA,nrow = 7,ncol = 4 ))
        })
      
      #################################
      weib_results_1 <- as.data.frame(cbind(summary(out.fit_nomed, df.resid=Inf)$table[,1],(summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit_nomed, df.resid=Inf)$table[,1]*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale)))
      weib_results_1$label <- rownames(weib_results_1)
      names(weib_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_1 <- rbind(names(weib_results_1),weib_results_1)
      names(weib_results_1) <- 1:7
      #################################
      weib_results_2 <- as.data.frame(cbind(summary(out.fit, df.resid=Inf)$table[,1],(summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit, df.resid=Inf)$table[,1]*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale)))
      weib_results_2$label <- rownames(weib_results_2)
      names(weib_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_2 <- rbind(names(weib_results_2),weib_results_2)
      names(weib_results_2) <- 1:7
      #################################
      weib_results_3 <- as.data.frame(cbind(summary(med_out.fit, df.resid=Inf)$table[,1],(summary(med_out.fit, df.resid=Inf)$table[,1]-summary(med_out.fit, df.resid=Inf)$table[,2]*1.96),
                                            (summary(med_out.fit, df.resid=Inf)$table[,1]+summary(med_out.fit, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(med_out.fit, df.resid=Inf)$table[,1]*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit, df.resid=Inf)$table[,1]+summary(med_out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit, df.resid=Inf)$table[,1]-summary(med_out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/med_out.fit$scale)))
      weib_results_3$label <- rownames(weib_results_3)
      names(weib_results_3) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_3 <- rbind(names(weib_results_3),weib_results_3)
      names(weib_results_3) <- 1:7
      #################################
      unweight_results_1 <- as.data.frame(cbind(summary(model_3_nomed)$table[,1],(summary(model_3_nomed)$table[,1]-summary(model_3_nomed)$table[,2]*1.96),
                                                (summary(model_3_nomed)$table[,1]+summary(model_3_nomed)$table[,2]*1.96),
                                                exp(summary(model_3_nomed)$table[,1]*(-1)/model_3_nomed$scale),
                                                exp((summary(model_3_nomed)$table[,1]+summary(model_3_nomed)$table[,2]*1.96)*(-1)/model_3_nomed$scale),
                                                exp((summary(model_3_nomed)$table[,1]-summary(model_3_nomed)$table[,2]*1.96)*(-1)/model_3_nomed$scale)))
      unweight_results_1$label <- rownames(unweight_results_1)
      names(unweight_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                     "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      unweight_results_1 <- rbind(names(unweight_results_1),unweight_results_1)
      names(unweight_results_1) <- 1:7
      #################################
      exp_m_results <-  as.data.frame(rbind(c(summary(med.fit, df.resid=Inf)$df.null+1,NA,NA,NA),
                                            cbind(coef(summary(med.fit, df.resid=Inf))[,1],
                                                  coef(summary(med.fit, df.resid=Inf))[,1]-1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                                  coef(summary(med.fit, df.resid=Inf))[,1]+1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                                  coef(summary(med.fit, df.resid=Inf))[,4])))
      row.names(exp_m_results)[1] <- "Raw_sample_size"
      names(exp_m_results) <- c("Beta","lowCI","highCI","pvalue")
      exp_m_results$label <- rownames(exp_m_results)
      exp_m_results <- cbind(exp_m_results,matrix(NA,nrow = nrow(exp_m_results),ncol=2))
      names(exp_m_results) <- 1:7
      #################################
      unweight_results_2 <- as.data.frame(cbind(summary(model_3)$table[,1],(summary(model_3)$table[,1]-summary(model_3)$table[,2]*1.96),
                                                (summary(model_3)$table[,1]+summary(model_3)$table[,2]*1.96),
                                                exp(summary(model_3)$table[,1]*(-1)/model_3$scale),
                                                exp((summary(model_3)$table[,1]+summary(model_3)$table[,2]*1.96)*(-1)/model_3$scale),
                                                exp((summary(model_3)$table[,1]-summary(model_3)$table[,2]*1.96)*(-1)/model_3$scale)))
      unweight_results_2$label <- rownames(unweight_results_2)
      names(unweight_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                     "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      unweight_results_2 <- rbind(names(unweight_results_2),unweight_results_2)
      names(unweight_results_2) <- 1:7
      #################################
      # cox_results_1 <- as.data.frame(cbind(coef(summary(model_1_nomed))[,1],(coef(summary(model_1_nomed))[,1]-coef(summary(model_1_nomed))[,4]*1.96),
      #                                      (coef(summary(model_1_nomed))[,1]+coef(summary(model_1_nomed))[,4]*1.96),coef(summary(model_1_nomed))[,2],
      #                                      exp(coef(summary(model_1_nomed))[,1]-coef(summary(model_1_nomed))[,4]*1.96),
      #                                      exp(coef(summary(model_1_nomed))[,1]+coef(summary(model_1_nomed))[,4]*1.96)))
      # cox_results_1$label <- rownames(cox_results_1)
      # names(cox_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
      #                           "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      # cox_results_1 <- rbind(names(cox_results_1),cox_results_1)
      # names(cox_results_1) <- 1:7
      # #################################
      # cox_results_2 <- as.data.frame(cbind(coef(summary(model_1))[,1],(coef(summary(model_1))[,1]-coef(summary(model_1))[,4]*1.96),
      #                                      (coef(summary(model_1))[,1]+coef(summary(model_1))[,4]*1.96),coef(summary(model_1))[,2],
      #                                      exp(coef(summary(model_1))[,1]-coef(summary(model_1))[,4]*1.96),
      #                                      exp(coef(summary(model_1))[,1]+coef(summary(model_1))[,4]*1.96)))
      # cox_results_2$label <- rownames(cox_results_2)
      # names(cox_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
      #                           "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      # cox_results_2 <- rbind(names(cox_results_2),cox_results_2)
      # names(cox_results_2) <- 1:7
      # #################################
      
      
      ####################################################################################################
      printout_table <-  as.data.frame(rbind(c("Causal mediation model:",rep(NA,6)),temp,rep(NA,7),rep(NA,7),rep(NA,7),
                                             c("Weibull model:",rep(NA,6)),weib_results_2,rep(NA,7),
                                             rep(NA,7),rep(NA,7),
                                             c("Weibull model without mediator:",rep(NA,6)),weib_results_1,rep(NA,7),
                                             c("Mediator Weibull model without exposure:",rep(NA,6)),weib_results_3,rep(NA,7),
                                             c("Exposure mediator model:",rep(NA,6)),exp_m_results,rep(NA,7),
                                             rep(NA,7),rep(NA,7),
                                             c("Unweighted Weibull model with mediator:",rep(NA,6)),unweight_results_2,rep(NA,7),
                                             c("Unweighted Weibull model without mediator:",rep(NA,6)),unweight_results_1,rep(NA,7)#,
                                             #c("Cox model:",rep(NA,6)),cox_results_2,rep(NA,7),
                                             #c("Cox model without mediator:",rep(NA,6)),cox_results_1 
                                             ))
      
      if (clock %in% c("HannumAge", "GDF15Mort")){
        wb <- createWorkbook()
      }
      addWorksheet(wb, paste0(clock))
      writeData(wb,x=printout_table,sheet = paste0(clock),colNames = F,
                rowNames = F)
    }
    saveWorkbook(wb, file=paste0("Results/SES main models for ",Group,"_",var,
                                 format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
  }
}
# dev.off()
####################################################################################################


####################################################################################################
## block 6.2: Bar plots for main models
####################################################################################################
library(ggplot2)

for (Group in loopgroup){
  pdf(paste0("Results/SES bar plots for the main models among ",Group,
             #"_07_13_2025.pdf"),
      format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
      
      width = 12,height = 5)
  
  if (Group %in% c("overallSample","female","male")){
    var_list <- c("blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
                  "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
                  "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
                  "lowwhitevshiwhite","hibluevshiwhite",
                  "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  } else {
    var_list <- c(#"blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
      "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
      "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
      "lowwhitevshiwhite","hibluevshiwhite",
      "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  }
  
  for (tt in var_list){
    #tt="blackvswhite"
    print(tt)
    
    dataforplots <-NULL
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/SES main models for ",Group,"_",tt,
                                         #"_07_13_2025.xlsx"),
                                  format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                  sheet = i)
      work <- as.data.frame(work[5,])
      names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value")
      work$set <- clock_list[i]
      
      dataforplots <- rbind(dataforplots,work)
    }
    
    if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                               " ","PhenoAge","GrimAgeMort","GrimAge2Mort",
                                               "   ","DunedinPoAm","      ",
                                               "YangCell","HorvathTelo","            ",
                                               "sedentary","hei","packyrs","drinker","       ",
                                               "waist2thigh","BMI","         ",
                                               "lbdtcsi","HDL","LDL","Glucose","CRP")
    
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA," "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"   "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"      "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"       "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"         "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"            "))
    } else {plot_clock_list <- clock_list}
    
    dataforplots <- as.data.frame(dataforplots)
    
    dataforplots$set <- factor(dataforplots$set,levels=plot_clock_list)
    dataforplots$sign <- ""
    dataforplots$sign[dataforplots$p_value<(0.05)&!is.na(dataforplots$p_value)] <- "*"
    dataforplots$sign2 <- ""
    dataforplots$sign2[dataforplots$p_value<(0.05/24)&!is.na(dataforplots$p_value)] <- "**"
    
    for (temp in c("Estimate","coef_lowCI","coef_highCI")){
      dataforplots[,temp] <- as.numeric(dataforplots[,temp])*100
    }
    print(range(c(dataforplots$coef_lowCI,dataforplots$coef_highCI),na.rm = T))
    dataforplots$starposition <- ifelse(dataforplots$coef_highCI<130,dataforplots$coef_highCI,140)
    
    
    print(ggplot(data=dataforplots,aes(x=set,y=Estimate))+
            coord_cartesian(ylim = c(-150,150))+
            theme_classic()+ xlab(paste0("Proportion of mediated effects of DNA methylation clocks and biomarkers on ",tt))+
            geom_errorbar(aes(x=set,y=Estimate,ymin=coef_lowCI,ymax=coef_highCI),
                          position = position_dodge(0.3), width= 0.4)+
            geom_point(aes(x=set,y=Estimate),position = position_dodge(0.3),size=2.5)+
            geom_hline(yintercept = 0)+geom_hline(yintercept = 50,linetype=3)+geom_hline(yintercept = 100,linetype=4)+
            geom_hline(yintercept = -50,linetype=3)+geom_hline(yintercept = -100,linetype=4)+
            geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 12,color="purple",linetype="dashed")+
            geom_vline(xintercept = 14,color="purple",linetype="dashed")+geom_vline(xintercept = 17,color="purple",linetype="dashed")+
            geom_vline(xintercept = 22,color="purple",linetype="dashed")+geom_vline(xintercept = 25,color="purple",linetype="dashed")+
            geom_text(aes(y=starposition,label=sign), vjust=-0.01, color="red",size = 10)+
            geom_text(aes(y=starposition,label=sign2), vjust=-0.01, color="red",size = 10)+
            labs(y="Proportion of mediated effects")+theme(axis.text.x = element_text(angle=30,hjust = 1))   )
  }
  dev.off()
}

####################################################################################################



####################################################################################################
## block 6.3: Bar plots putting together for main models
# ####################################################################################################
# library(ggplot2)
# pdf(paste0("Results/SES bar plots for the main models putting together",format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
#     width = 20,height = 5)
# 
# for (tt in SES_tablelist){
#   #tt="blackvswhite"
#   print(tt)
#   jj=1
#   if (jj==1){title <- "race"  } else
#     if (jj==2){title <- "education" } else
#       if (jj==3){title <- "income" } else
#         if (jj==4){title <- "occupation" }
#   
#   dataforplots <- NULL
#   for (tttt_tt in tt){
#     for (i in 1:length(clock_list)){
#       #i=1;i=13
#       print(i)
#       print(clock_list[i])
#       work <-  readxl::read_excel(paste0("Results/SES main models for ",tttt_tt, format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
#                                   sheet = i)
#       work <- as.data.frame(work[5,])
#       names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value")
#       work$set <- clock_list[i]
#       work$var <- tttt_tt
#       
#       dataforplots <- rbind(dataforplots,work)
#     }
#   }
#   
#   if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
#                                              " ","PhenoAge","GrimAgeMort","GrimAge2Mort",
#                                              "   ","DunedinPoAm","      ",
#                                              "YangCell","HorvathTelo","            ",
#                                              "sedentary","hei","packyrs","drinker","       ",
#                                              "waist2thigh","BMI","         ",
#                                              "lbdtcsi","HDL","LDL","Glucose","CRP")
#   
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA," ",tttt_tt))
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"   ",tttt_tt))
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"      ",tttt_tt))
#   
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"       ",tttt_tt))
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"         ",tttt_tt))
#   dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"            ",tttt_tt))
#   } else {plot_clock_list <- clock_list}
#   
#   dataforplots <- as.data.frame(dataforplots)
#   
#   dataforplots$set <- factor(dataforplots$set,levels=plot_clock_list)
#   dataforplots$var <- factor(dataforplots$var,levels=tt)
#   
#   dataforplots$sign <- ""
#   dataforplots$sign[dataforplots$p_value<(0.05)&!is.na(dataforplots$p_value)] <- "*"
#   
#   for (temp in c("Estimate","coef_lowCI","coef_highCI")){
#     dataforplots[,temp] <- as.numeric(dataforplots[,temp])*100
#   }
#   print(range(c(dataforplots$coef_lowCI,dataforplots$coef_highCI),na.rm = T))
#   
#   print(ggplot(data=dataforplots,aes(x=set,y=Estimate,colour = var))+
#           coord_cartesian(ylim = c(-150,150))+
#           theme_classic()+ xlab(paste0("Proportion of mediated effects of DNA methylation clocks and biomarkers on ",title))+
#           geom_errorbar(aes(x=set,y=Estimate,ymin=coef_lowCI,ymax=coef_highCI),
#                         position = position_dodge(0.3), width= 0.4)+
#           geom_point(aes(x=set,y=Estimate),position = position_dodge(0.3),size=2.5)+
#           geom_hline(yintercept = 0)+geom_hline(yintercept = 50,linetype=3)+geom_hline(yintercept = 100,linetype=4)+
#           geom_hline(yintercept = -50,linetype=3)+geom_hline(yintercept = -100,linetype=4)+
#           geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 12,color="purple",linetype="dashed")+
#           geom_vline(xintercept = 14,color="purple",linetype="dashed")+geom_vline(xintercept = 17,color="purple",linetype="dashed")+
#           geom_vline(xintercept = 22,color="purple",linetype="dashed")+geom_vline(xintercept = 25,color="purple",linetype="dashed")+
#           scale_color_manual(values = c("red","blue","orange","grey"))+
#           geom_text(aes(y=coef_highCI,label=sign), vjust=-0.01, color="red",size = 10)+
#           labs(y="Proportion of mediated effects")+theme(axis.text.x = element_text(angle=30,hjust = 1))   )
#   jj=jj+1
# }
# dev.off()
# ####################################################################################################


####################################################################################################
##  block 6.4: Printing out tables for main models
####################################################################################################
wb <- createWorkbook()

for (Group in loopgroup){
  dataforplots_final <- NULL
  
  if (Group %in% c("overallSample","female","male")){
    var_list <- c("blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
                  "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
                  "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
                  "lowwhitevshiwhite","hibluevshiwhite",
                  "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  } else {
    var_list <- c(#"blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
      "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
      "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
      "lowwhitevshiwhite","hibluevshiwhite",
      "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  }
  
  for (tt in var_list){
    #tt <- "blackvswhite"
    
    dataforplots <- NULL
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/SES main models for ",Group,"_",tt, 
                                         format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                  sheet = i)
      
      work <- as.data.frame(work[5,c(1:4,7)])
      names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value","sample_size")
      work$clock <- clock_list[i]
      
      dataforplots <- rbind(dataforplots,work)
    }
    dataforplots$set <- tt
    dataforplots$gap <- NA
    
    if (is.null(dataforplots_final)){
      dataforplots_final <- dataforplots
    } else {dataforplots_final<- cbind(dataforplots_final,dataforplots)}
  }
  addWorksheet(wb, paste0(Group,"mediation"))
  writeData(wb,x=dataforplots_final,sheet = paste0(Group,"mediation"),
            rowNames = F)
}

dataforsurvival <- NULL
for (SES_pp in var_list){
  #SES_pp <- "blackvswhite"
  
  work <-  readxl::read_excel(paste0("Results/SES main models for overallSample_",SES_pp, 
                                     format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                              sheet = 1)
  
  work <- as.data.frame(work[25,])
  work <- work[,c(4:6)]
  names(work) <- c("HR","HR_lowCI","HR_highCI")
  work$set <- SES_pp
  
  dataforsurvival <- rbind(dataforsurvival,work)
}
addWorksheet(wb, "Survival")
writeData(wb,x=dataforsurvival,sheet = "Survival",
          rowNames = F)


dataformediator_final <- NULL
for (tt in var_list){
  #tt <- "Smoker"
  
  dataformediator <- NULL
  for (i in 1:length(clock_list)){
    #i=1;i=13
    print(i)
    print(clock_list[i])
    work <-  readxl::read_excel(paste0("Results/SES main models for overallSample_",tt, 
                                       format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                sheet = i)
    work <- as.data.frame(work[35,c(4:6,7)])
    names(work) <- c("HR","HR_lowCI","HR_highCI","clock")
    
    dataformediator <- rbind(dataformediator,work)
  }
  dataformediator$set <- tt
  dataformediator$gap <- NA
  
  if (is.null(dataformediator_final)){
    dataformediator_final <- dataformediator
  } else {dataformediator_final<- cbind(dataformediator_final,dataformediator)}
}
addWorksheet(wb, paste0("mediator_survival"))
writeData(wb,x=dataformediator_final,sheet = paste0("mediator_survival"),
          rowNames = F)
saveWorkbook(wb, file=paste0("Results/main models print tables.xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
## block 6.5: Survival mediation specific-cause models
####################################################################################################
specific_reason_list <- c("hrtmort","canmort")#"accimort","diabmort","cermort","resmort"

for (specific_reason in specific_reason_list){
  #specific_reason <- "accimort"; specific_reason <- "diabmort"
  print(specific_reason)
  
  if (specific_reason %in% c("hrtmort","canmort")){
    var_list <-  c("blackvswhite","latinovswhite","latinovsblack","otherRacevswhite",#"otherHispvswhite","otherRacevswhite",
                   "lthsvscoll","hsvscoll","somecollvscoll","hsvssomecoll",
                   "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
                   "lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  } else if (specific_reason %in% c("accimort","diabmort","cermort")){
    var_list <-  c("blackvswhite","latinovswhite","latinovsblack",
                   #"lthsvscoll","lths_hsvscoll",
                   "below1vsabove5","below2vsabove5","lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite","noworkvshiwhite")#
  } else if (specific_reason %in% c("resmort")){
    var_list <-  c("blackvswhite","latinovswhite","latinovsblack",
                   "lthsvscoll","lths_hsvscoll",
                   "below1vsabove5","below2vsabove5","lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite")#,"noworkvshiwhite"
  }
  
  for (var in var_list){
    # var <- "latinovswhite"; 
    
    explore_surveydata <- get("svyNHEanalysis_multiimpu")  ## get imputed survey object
    # explore_surveydata <- get("svyNHEanalysis") ## get completed cases survey object
    
    
    explore_surveydata$variables$test <- explore_surveydata$variables[,var]
    explore_surveydata <- subset(explore_surveydata,!is.na(test)&WTDN4YR>0)
    #explore_surveydata <- subset(explore_surveydata,!is.na(forborn))
    
    print(var)
    # table(explore_surveydata$variables[,var],explore_surveydata$variables[,specific_reason])
    # table(paste0(explore_surveydata$variables[,var],explore_surveydata$variables[,specific_reason],explore_surveydata$variables$dead))
    #if  (min(unlist(table(explore_surveydata$variables[,var],explore_surveydata$variables[,specific_reason]))) >15 ){
    
    ##############   models
    model_printtable <-  NULL
    
    for (clock in clock_list){
      #clock <- "HorvathAge"; clock <- "DNAmGrimAge2"
      clock_accl <- paste0(clock,"_DNAmAge_Acc")
      
      print(clock)
      #### mediator model
      med.fit <- survey::svyglm(design=explore_surveydata,paste0(clock,"~",var,covariates),family="gaussian")
      summary(med.fit)#summary(med.fit)$df.null
      
      #### mediator_outcome model
      med_out.fit <- survey::svysurvreg(design=explore_surveydata,
                                        as.formula(paste0("Surv(deathage,",specific_reason,"==1)~",clock,
                                                          covariates)), dist="weibull")
      summary(med_out.fit)
      
      #### outcome model without mediators
      #table(rawanalysis$deathage)
      model_1_nomed <- survey::svycoxph(design=explore_surveydata,
                                        as.formula(paste0("Surv(deathage,",specific_reason,"==1)~",var,
                                                          covariates)))
      model_1 <- survey::svycoxph(design=explore_surveydata,
                                  as.formula(paste0("Surv(deathage,",specific_reason,"==1)~",var,"+",clock,
                                                    covariates)))
      #summary(model_1)
      
      
      competing_model_nomed <- cmprsk::crr(cencode = 0,failcode = 1,ftime = rawanalysis$deathage,
                                           fstatus = rawanalysis[,paste0(specific_reason,"_competing")],
                                           cov1 = rawanalysis[,c(var,"female","age","agesq","forborn")])
      competing_model <- cmprsk::crr(cencode = 0,failcode = 1,ftime = rawanalysis$deathage,
                                     fstatus = rawanalysis[,paste0(specific_reason,"_competing")],
                                     cov1 = rawanalysis[,c(var,clock,"female","age","agesq","forborn")])
      summary(competing_model)
      
      
      out.fit_nomed <- survey::svysurvreg(design=explore_surveydata,
                                          as.formula(paste0("Surv(deathage,",specific_reason,"==1)~",var,
                                                            covariates)), dist="weibull")
      tryCatch({      
        out.fit <- survey::svysurvreg(design=explore_surveydata,
                                      as.formula(paste0("Surv(deathage,",specific_reason,"==1)~",var,"+",clock,
                                                        covariates)), dist="weibull")
        summary(out.fit)
        
        
        #### causal mediation models
        
        med.out <- mediate(med.fit,out.fit,treat=var,mediator=clock)
        summary(med.out)
        
        ###################################################################################################
        #### Create table
        temp <- as.data.frame(extract_mediation_summary(summary(med.out)))[c(5,8:10),]
        temp$label <- rownames(temp)
        temp <- rbind(names(temp),temp)
        temp <- cbind(temp,matrix(NA,nrow = nrow(temp),ncol=1))
        temp$sample <- med.out$nobs
        names(temp) <- 1:7
      },error=function(e){
        cat("Error",conditionMessage(e),"\n");temp <- as.data.frame(matrix(data = NA,nrow = 7,ncol = 4 ))
      })
      
      ###################################################################################################
      #### Create table
      #################################
      weib_results_1 <- as.data.frame(cbind(summary(out.fit_nomed, df.resid=Inf)$table[,1],(summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit_nomed, df.resid=Inf)$table[,1]*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale)))
      weib_results_1$label <- rownames(weib_results_1)
      names(weib_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_1 <- rbind(names(weib_results_1),weib_results_1)
      names(weib_results_1) <- 1:7
      #################################
      weib_results_2 <- as.data.frame(cbind(summary(out.fit, df.resid=Inf)$table[,1],(summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit, df.resid=Inf)$table[,1]*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale)))
      weib_results_2$label <- rownames(weib_results_2)
      names(weib_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_2 <- rbind(names(weib_results_2),weib_results_2)
      names(weib_results_2) <- 1:7
      #################################
      weib_results_3 <- as.data.frame(cbind(summary(med_out.fit, df.resid=Inf)$table[,1],(summary(med_out.fit, df.resid=Inf)$table[,1]-summary(med_out.fit, df.resid=Inf)$table[,2]*1.96),
                                            (summary(med_out.fit, df.resid=Inf)$table[,1]+summary(med_out.fit, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(med_out.fit, df.resid=Inf)$table[,1]*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit, df.resid=Inf)$table[,1]+summary(med_out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit, df.resid=Inf)$table[,1]-summary(med_out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/med_out.fit$scale)))
      weib_results_3$label <- rownames(weib_results_3)
      names(weib_results_3) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_3 <- rbind(names(weib_results_3),weib_results_3)
      names(weib_results_3) <- 1:7
      #################################
      cox_results_1 <- as.data.frame(cbind(coef(summary(model_1_nomed))[,1],(coef(summary(model_1_nomed))[,1]-coef(summary(model_1_nomed))[,4]*1.96),
                                           (coef(summary(model_1_nomed))[,1]+coef(summary(model_1_nomed))[,4]*1.96),coef(summary(model_1_nomed))[,2],
                                           exp(coef(summary(model_1_nomed))[,1]-coef(summary(model_1_nomed))[,4]*1.96),
                                           exp(coef(summary(model_1_nomed))[,1]+coef(summary(model_1_nomed))[,4]*1.96)))
      cox_results_1$label <- rownames(cox_results_1)
      names(cox_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      cox_results_1 <- rbind(names(cox_results_1),cox_results_1)
      names(cox_results_1) <- 1:7
      #################################
      cox_results_2 <- as.data.frame(cbind(coef(summary(model_1))[,1],(coef(summary(model_1))[,1]-coef(summary(model_1))[,4]*1.96),
                                           (coef(summary(model_1))[,1]+coef(summary(model_1))[,4]*1.96),coef(summary(model_1))[,2],
                                           exp(coef(summary(model_1))[,1]-coef(summary(model_1))[,4]*1.96),
                                           exp(coef(summary(model_1))[,1]+coef(summary(model_1))[,4]*1.96)))
      cox_results_2$label <- rownames(cox_results_2)
      names(cox_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      cox_results_2 <- rbind(names(cox_results_2),cox_results_2)
      names(cox_results_2) <- 1:7
      #################################
      competing_results_1 <- as.data.frame(cbind(summary(competing_model_nomed)$coef[,1],(summary(competing_model_nomed)$coef[,1]-summary(competing_model_nomed)$coef[,3]*1.96),
                                                 (summary(competing_model_nomed)$coef[,1]+summary(competing_model_nomed)$coef[,3]*1.96),exp(summary(competing_model_nomed)$coef[,1]),
                                                 exp(summary(competing_model_nomed)$coef[,1]-summary(competing_model_nomed)$coef[,3]*1.96),
                                                 exp(summary(competing_model_nomed)$coef[,1]+summary(competing_model_nomed)$coef[,3]*1.96)))
      competing_results_1$label <- rownames(competing_results_1)
      names(competing_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                      "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      competing_results_1 <- rbind(names(competing_results_1),competing_results_1)
      names(competing_results_1) <- 1:7
      #################################
      competing_results_2 <- as.data.frame(cbind(summary(competing_model)$coef[,1],(summary(competing_model)$coef[,1]-summary(competing_model)$coef[,3]*1.96),
                                                 (summary(competing_model)$coef[,1]+summary(competing_model)$coef[,3]*1.96),exp(summary(competing_model)$coef[,1]),
                                                 exp(summary(competing_model)$coef[,1]-summary(competing_model)$coef[,3]*1.96),
                                                 exp(summary(competing_model)$coef[,1]+summary(competing_model)$coef[,3]*1.96)))
      competing_results_2$label <- rownames(competing_results_2)
      names(competing_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                      "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      competing_results_2 <- rbind(names(competing_results_2),competing_results_2)
      names(competing_results_2) <- 1:7
      #################################
      exp_m_results <-  as.data.frame(rbind(c(summary(med.fit, df.resid=Inf)$df.null+1,NA,NA,NA),
                                            cbind(coef(summary(med.fit, df.resid=Inf))[,1],
                                                  coef(summary(med.fit, df.resid=Inf))[,1]-1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                                  coef(summary(med.fit, df.resid=Inf))[,1]+1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                                  coef(summary(med.fit, df.resid=Inf))[,4])))
      row.names(exp_m_results)[1] <- "Raw_sample_size"
      names(exp_m_results) <- c("Beta","lowCI","highCI","pvalue")
      exp_m_results$label <- rownames(exp_m_results)
      exp_m_results <- cbind(exp_m_results,matrix(NA,nrow = nrow(exp_m_results),ncol=2))
      names(exp_m_results) <- 1:7
      #################################
      
      ###################################################################################################
      printout_table <-  as.data.frame(rbind(c("Causal mediation model:",rep(NA,6)),temp,rep(NA,7),rep(NA,7),rep(NA,7),
                                             c("Weibull model:",rep(NA,6)),weib_results_2,rep(NA,7),
                                             c("Weibull model without mediator:",rep(NA,6)),weib_results_1,rep(NA,7),
                                             c("Weibull model without exposure:",rep(NA,6)),weib_results_3,rep(NA,7),
                                             c("Exposure mediator model:",rep(NA,6)),exp_m_results,rep(NA,7),
                                             rep(NA,7),rep(NA,7),
                                             #c("Unweighted Weibull model with mediator:",rep(NA,6)),unweight_results_2,rep(NA,7),
                                             #c("Unweighted Weibull model without mediator:",rep(NA,6)),unweight_results_1,rep(NA,7),
                                             c("Cox model:",rep(NA,6)),cox_results_2,rep(NA,7),
                                             c("Cox model without mediator:",rep(NA,6)),cox_results_1,rep(NA,7),
                                             c("Competing model:",rep(NA,6)),competing_results_2,rep(NA,7),
                                             c("Competing model without mediator:",rep(NA,6)),competing_results_1 ))
      
      if (clock %in% c("HannumAge", "GDF15Mort")){
        wb <- createWorkbook()
      }
      addWorksheet(wb, paste0(clock))
      writeData(wb,x=printout_table,sheet = paste0(clock),colNames = F,
                rowNames = F)
    }
    saveWorkbook(wb, file=paste0("Results/SES specific reason models for ",var," with specific death of ",specific_reason,
                                 format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
  }
}
# dev.off()
####################################################################################################

####################################################################################################
## block 6.6: Bar plots for survival mediation of specific-cause models
####################################################################################################
library(ggplot2)
for (specific_reason in specific_reason_list){
  #specific_reason <- "cermort"
  print(specific_reason)
  pdf(paste0("Results/SES bar plots for specific cause models with specific death of ",
             specific_reason,format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
      width = 12,height = 5)
  
  
  if (specific_reason %in% c("hrtmort","canmort","accimort","diabmort")){
    var_list <-  c("blackvswhite","latinovswhite","otherRacevswhite",#"otherHispvswhite","latinovsblack",
                   "lthsvscoll","hsvscoll","somecollvscoll",
                   "below1vsabove5","pir_1_2vsabove5","pir_2_5vsabove5",
                   "lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite","noworkvshiwhite","Smoker")#
  } else if (specific_reason %in% c("cermort")){
    var_list <-  c("blackvswhite","latinovswhite","latinovsblack",
                   #"lthsvscoll","lths_hsvscoll",
                   "below1vsabove5","below2vsabove5","lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite","noworkvshiwhite")#
  } else if (specific_reason %in% c("resmort")){
    var_list <-  c("blackvswhite","latinovswhite","latinovsblack",
                   "lthsvscoll","lths_hsvscoll",
                   "below1vsabove5","below2vsabove5","lowwhitevshiwhite","hibluevshiwhite",
                   "lowbluevshiwhite")#,"noworkvshiwhite"
  }
  
  for (tt in var_list){
    #tt="blackvswhite"
    print(tt)
    dataforplots <-NULL
    
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/SES specific reason models for ",tt," with specific death of ",
                                         #specific_reason,"_11_07_2024.xlsx"),
                                         specific_reason,format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                  sheet = i)
      work <- as.data.frame(work[5,])
      names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value")
      work$set <- clock_list[i]
      
      dataforplots <- rbind(dataforplots,work)
    }
    
    
    if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                               " ","PhenoAge","GrimAgeMort","GrimAge2Mort",
                                               "   ","DunedinPoAm","      ",
                                               "YangCell","HorvathTelo","            ",
                                               "sedentary","hei","packyrs","drinker","       ",
                                               "waist2thigh","BMI","         ",
                                               "lbdtcsi","HDL","LDL","Glucose","CRP")
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA," "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"   "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"      "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"       "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"         "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"            "))
    } else {plot_clock_list <- clock_list}
    
    
    dataforplots <- as.data.frame(dataforplots)
    
    dataforplots$set <- factor(dataforplots$set,levels=plot_clock_list)
    dataforplots$sign <- ""
    dataforplots$sign[dataforplots$p_value<(0.05)&!is.na(dataforplots$p_value)] <- "*"
    
    
    for (temp in c("Estimate","coef_lowCI","coef_highCI")){
      dataforplots[,temp] <- as.numeric(dataforplots[,temp])*100
    }
    print(range(c(dataforplots$coef_lowCI,dataforplots$coef_highCI),na.rm = T))
    dataforplots$starposition <- ifelse(dataforplots$coef_highCI<130,dataforplots$coef_highCI,140)
    
    print(ggplot(data=dataforplots,aes(x=set,y=Estimate))+
            theme_classic()+ xlab(paste0("Proportion of mediated effects of DNA methylation clocks and biomarkers on ",tt,
                                         " with specific death of ",specific_reason))+
            coord_cartesian(ylim = c(-150,150))+
            geom_errorbar(aes(x=set,y=Estimate,ymin=coef_lowCI,ymax=coef_highCI),
                          position = position_dodge(0.3), width= 0.4)+
            geom_point(aes(x=set,y=Estimate),position = position_dodge(0.3),size=2.5)+
            geom_hline(yintercept = 0)+geom_hline(yintercept = 50,linetype=3)+geom_hline(yintercept = 100,linetype=4)+
            geom_hline(yintercept = -50,linetype=3)+geom_hline(yintercept = -100,linetype=4)+
            geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 12,color="purple",linetype="dashed")+
            geom_vline(xintercept = 14,color="purple",linetype="dashed")+geom_vline(xintercept = 17,color="purple",linetype="dashed")+
            geom_vline(xintercept = 22,color="purple",linetype="dashed")+geom_vline(xintercept = 25,color="purple",linetype="dashed")+
            geom_text(aes(y=starposition,label=sign), vjust=-0.01, color="red",size = 10)+
            labs(y="Proportion of mediated effects")+theme(axis.text.x = element_text(angle=30,hjust = 1))   )
  }
  dev.off()
}
###################################################################################################


####################################################################################################
## block 6.7: Bar plots putting together for specific-cause models
####################################################################################################
library(ggplot2)


specific_reason_list <- c("hrtmort","canmort")#,"accimort","diabmort","cermort","resmort"

for (specific_reason in specific_reason_list){
  pdf(paste0("Results/SES bar plots for the specific cause_",specific_reason," models putting together",format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
      width = 20,height = 5)
  
  for (tt in SES_tablelist){
    #tt="blackvswhite"
    print(tt)
    jj=1
    if (jj==1){title <- "race"  } else
      if (jj==2){title <- "education" } else
        if (jj==3){title <- "income" } else
          if (jj==4){title <- "occupation" }
    
    dataforplots <- NULL
    for (tttt_tt in tt){
      for (i in 1:length(clock_list)){
        #i=1;i=13
        print(i)
        print(clock_list[i])
        work <-  readxl::read_excel(paste0("Results/SES specific reason models for ",tttt_tt," with specific death of ",specific_reason,
                                           format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                    sheet = i)
        work <- as.data.frame(work[5,])
        names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value")
        work$set <- clock_list[i]
        work$var <- tttt_tt
        
        dataforplots <- rbind(dataforplots,work)
      }
    }
    
    if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                               " ","PhenoAge","GrimAgeMort","GrimAge2Mort",
                                               "   ","DunedinPoAm","      ",
                                               "YangCell","HorvathTelo","            ",
                                               "sedentary","hei","packyrs","drinker","       ",
                                               "waist2thigh","BMI","         ",
                                               "lbdtcsi","HDL","LDL","Glucose","CRP")
    
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA," ",tttt_tt))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"   ",tttt_tt))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"      ",tttt_tt))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"       ",tttt_tt))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"         ",tttt_tt))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"            ",tttt_tt))
    } else {plot_clock_list <- clock_list}
    
    dataforplots <- as.data.frame(dataforplots)
    
    dataforplots$set <- factor(dataforplots$set,levels=plot_clock_list)
    dataforplots$var <- factor(dataforplots$var,levels=tt)
    
    dataforplots$sign <- ""
    dataforplots$sign[dataforplots$p_value<(0.05)&!is.na(dataforplots$p_value)] <- "*"
    dataforplots$sign2 <- ""
    dataforplots$sign2[dataforplots$p_value<(0.05/24)&!is.na(dataforplots$p_value)] <- "**"
    
    for (temp in c("Estimate","coef_lowCI","coef_highCI")){
      dataforplots[,temp] <- as.numeric(dataforplots[,temp])*100
    }
    print(range(c(dataforplots$coef_lowCI,dataforplots$coef_highCI),na.rm = T))
    
    print(ggplot(data=dataforplots,aes(x=set,y=Estimate,colour = var))+
            coord_cartesian(ylim = c(-150,380))+
            theme_classic()+ xlab(paste0("Proportion of mediated effects of DNA methylation clocks and biomarkers on ",title))+
            geom_errorbar(aes(x=set,y=Estimate,ymin=coef_lowCI,ymax=coef_highCI),
                          position = position_dodge(0.3), width= 0.4)+
            geom_point(aes(x=set,y=Estimate),position = position_dodge(0.3),size=2.5)+
            geom_hline(yintercept = 0)+geom_hline(yintercept = 50,linetype=3)+geom_hline(yintercept = 100,linetype=4)+
            geom_hline(yintercept = -50,linetype=3)+geom_hline(yintercept = -100,linetype=4)+
            geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 12,color="purple",linetype="dashed")+
            geom_vline(xintercept = 14,color="purple",linetype="dashed")+geom_vline(xintercept = 17,color="purple",linetype="dashed")+
            geom_vline(xintercept = 22,color="purple",linetype="dashed")+geom_vline(xintercept = 25,color="purple",linetype="dashed")+
            scale_color_manual(values = c("red","blue","orange","grey"))+
            geom_text(aes(y=coef_highCI,label=sign), vjust=-0.01, color="red",size = 10)+
            geom_text(aes(y=coef_highCI,label=sign2), vjust=-0.01, color="red",size = 10)+
            labs(y="Proportion of mediated effects")+theme(axis.text.x = element_text(angle=30,hjust = 1))   )
    jj=jj+1
  }
  dev.off()
}
#####################################################################################################



####################################################################################################
##  block 6.8: Printing out tables for specific-cause models
####################################################################################################
specific_reason_list <- c("hrtmort","canmort")#,"accimort","diabmort","cermort","resmort"

wb <- createWorkbook()
for (specific_reason in specific_reason_list){
  dataforplots_final <- NULL
  
  for (tt in var_list){
    #tt <- "blackvswhite"
    
    dataforplots <- NULL
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/SES specific reason models for ",tt," with specific death of ",specific_reason,
                                         format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                  sheet = i)
      work <- as.data.frame(work[5,c(1:4,7)])
      names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value","sample_size")
      work$clock <- clock_list[i]
      
      dataforplots <- rbind(dataforplots,work)
    }
    dataforplots$set <- tt
    dataforplots$gap <- NA
    
    if (is.null(dataforplots_final)){
      dataforplots_final <- dataforplots
    } else {dataforplots_final<- cbind(dataforplots_final,dataforplots)}
  }
  addWorksheet(wb, paste0(specific_reason,"_mediation"))
  writeData(wb,x=dataforplots_final,sheet = paste0(specific_reason,"_mediation"),
            rowNames = F)
}

dataforsurvival_final <- NULL
for (specific_reason in specific_reason_list){
  
  dataforsurvival <- NULL
  for (SES_pp in var_list){
    #SES_pp <- "Smoker"
    
    work <-  readxl::read_excel(paste0("Results/SES specific reason models for ",SES_pp," with specific death of ",specific_reason,
                                       format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                sheet = 1)
    
    work <- as.data.frame(work[23,])
    work <- work[,c(4:7)]
    names(work) <- c("HR","HR_lowCI","HR_highCI","set")
    
    dataforsurvival <- rbind(dataforsurvival,work)
  }
  dataforsurvival$new <- specific_reason
  dataforsurvival$gap <- NA
  
  if (is.null(dataforsurvival_final)){
    dataforsurvival_final <- dataforsurvival
  } else {dataforsurvival_final<- cbind(dataforsurvival_final,dataforsurvival)}
}
addWorksheet(wb, "Survival")
writeData(wb,x=dataforsurvival_final,sheet = "Survival",
          rowNames = F)


for (specific_reason in specific_reason_list){
  dataformediator_final <- NULL
  for (tt in var_list){
    #tt <- "blackvswhite"
    
    dataformediator <- NULL
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/SES specific reason models for ",tt," with specific death of ",specific_reason,
                                         format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                  sheet = i)
      work <- as.data.frame(work[33,c(4:6,7)])
      names(work) <- c("HR","HR_lowCI","HR_highCI","clock")
      
      dataformediator <- rbind(dataformediator,work)
    }
    dataformediator$set <- tt
    dataformediator$gap <- NA
    
    if (is.null(dataformediator_final)){
      dataformediator_final <- dataformediator
    } else {dataformediator_final<- cbind(dataformediator_final,dataformediator)}
  }
  addWorksheet(wb, paste0("mediator_survival_",specific_reason))
  writeData(wb,x=dataformediator_final,sheet = paste0("mediator_survival_",specific_reason),
            rowNames = F)
}


saveWorkbook(wb, file=paste0("Results/specific reason models print tables.xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
##  block 6.9: sample descriptive for categorical variables without missing by group
####################################################################################################
variable_list_sample_categorical <- c("white","black","mexican","otherHisp","otherRace","",
                                      "lths","hs","somecoll","coll","",
                                      "pir_below1","pir_1_2","pir_2_5","pir_above5","",
                                      "hiwhite","lowwhite","hiblue","lowblue","nowork","",
                                      "female","male","", "nativity","forborn","",
                                      "active","sedentary","","abstainer","drinker")#,"","dead","hrtmort","canmort"


## descriptive tables
descriptive_categorical_bygroup <- NULL

for (tt in c("overallSample","dead","undead")){
  #tt="undead"
  print(tt)
  i=1;descriptive_categorical <- NULL
  
  for (var in variable_list_sample_categorical){
    # var="white"
    
    explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
    # explore_dataset <- get("rawanalysis") ## get completed cases dataset
    
    explore_dataset <- explore_dataset[explore_dataset[,tt]==1,]
    
    if (var==""){
      temp_2 <- as.data.frame(matrix(nrow=1,ncol = 4))
    } else{
      #print(var)
      
      temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
      temp_1$Var1 <- as.character(temp_1$Var1)
      temp_1$Var1[temp_1$Var1==1] <- "Yes"
      temp_1$Var1[temp_1$Var1==0] <- "No"
      if (!"Yes" %in% temp_1$Var1){temp_1 <- rbind(temp_1,c("Yes",0))}
      
      temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])
      
      temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))
      
      temp_1 <- temp_1[order(temp_1$Var1),]
      temp_1 <- temp_1[!is.na(temp_1$Var1),]
      
      temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
    }
    
    if (!is.null(names(descriptive_categorical))){
      names(temp_2) <- names(descriptive_categorical)}
    
    descriptive_categorical <- rbind(descriptive_categorical,temp_2)
    
    i=i+1
  }
  
  descriptive_categorical$variables <- descriptive_categorical$`c(var, rep("", nrow(temp_1) - 1))`
  descriptive_categorical$`c(var, rep("", nrow(temp_1) - 1))` <- NULL
  descriptive_categorical$Var1 <- NULL
  names(descriptive_categorical) <- paste0( tt,"_", c("count","percentage"))
  
  if (is.null(descriptive_categorical_bygroup)) {
    descriptive_categorical_bygroup <- descriptive_categorical } else {
      descriptive_categorical_bygroup <- cbind(descriptive_categorical_bygroup,descriptive_categorical)
    }
}

wb <- createWorkbook()
addWorksheet(wb, "Table 1_nonmissing_raw")
writeData(wb,x=descriptive_categorical_bygroup,sheet = "Table 1_nonmissing_raw",
          rowNames = F)

####################################################################################################
##  block 6.10: sample descriptive for continuous variables by group
####################################################################################################
variable_list_sample_continuous <- c("age","hei","packyrs","waist2thigh","BMI",
                                     "lbdtcsi","HDL","LDL","Glucose","CRP",
                                     "lbxlypct","lbxmopct","lbxnepct","lbxeopct","lbxbapct")
## descriptive tables
descriptive_continuous_bygroup <- as.data.frame(variable_list_sample_continuous)

for (tt in c("overallSample","dead","undead")){
  i=1; descriptive_continuous <- NULL
  for (var in variable_list_sample_continuous){
    #var <- "age"

    explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
    # explore_dataset <- get("rawanalysis") ## get completed cases dataset

    explore_dataset <- explore_dataset[explore_dataset[,tt]==1,]

    print(var)

    if (length(table(is.na(explore_dataset[,var])))>1) {
      temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
    } else {temp_0 <- 0}

    temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
    temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
    # temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
    #                    c(.25,.5,.75))
    temp_4 <- cbind(var,table(is.na(explore_dataset[,var]))[1],#as.data.frame(temp_0),
                    as.data.frame(temp_1),as.data.frame(temp_2)#,
                    #t(as.data.frame(temp_3)),
    )

    descriptive_continuous <- rbind(descriptive_continuous,temp_4)
    i=i+1
  }
  descriptive_continuous$var <- NULL
  names(descriptive_continuous) <- paste0( tt,"_", c(#"missing_rate",
    "count","mean","se"#,
    #"lower_quartile","median","upper_quartile",
  ))


  descriptive_continuous_bygroup <- cbind(descriptive_continuous_bygroup,descriptive_continuous)
}

addWorksheet(wb, "Table 1_continuous")
writeData(wb,x=descriptive_continuous_bygroup,sheet = "Table 1_continuous",
          rowNames = F)
saveWorkbook(wb, file=paste0("Results/Table 1 by group",
                             format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
tt_toc <- proc.time()
tt_tic_toc_time <- tt_toc-tt_tic
tt_tic_toc_time/3600
####################################################################################################







