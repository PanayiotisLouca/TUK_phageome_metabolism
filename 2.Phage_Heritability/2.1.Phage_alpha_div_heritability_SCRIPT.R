## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  Estimate heritability of alpha diversity of phageome 

## Clear environment 
rm(list = ls()) 

## Set seed 
set.seed(1234)

## load up packages: 

### core 
library(tidyverse)
library(skimr)
library(readxl)

# mixed effect modelling 
library(lme4)
library(lmerTest)
library(broom.mixed)

# mets for heritability 
library(mets)
library(lava)

    # -------------------------------------------------------------------------- # 
  
      # ************************* # 
      #   IMPORT & PREP DATA   ---- 
      # ************************* # 
   
   # -------------------------------------------------------------------------- #  

      ##   Dataset ---- 
path <-
  file.path(
    "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Data/phages_microbes_metabs_updt_DATASET.rds"
  )

df <- read_rds(path) %>% 
  as.data.frame()

covar_columns = c('gm_age','gm_sex','gm_BMI','gm_batch', 'fid')

# -------------------------------------------------------------------------- #  

df <- df %>%
  mutate(
    gm_zygosity = factor(gm_zygosity),
    gm_sex = factor(gm_sex),
    gm_batch = factor(gm_batch))

    # -------------------------------------------------------------------------- #    

# create new dataset with just needed variables 
df_analyse <- df %>%
  select(
    iid,
    all_of(covar_columns),
    fid,
    gm_zygosity,
    # phage alpha div 
    phage_alpha_div_shannon,
    phage_alpha_div_taxonomic_richness
  )

    # -------------------------------------------------------------------------- #  
  
            ####   RUN FOR SHANNON DIVERISTY  ---- 

# split data into phenotype data and outcome data 
phen <- data.frame(
  fid = df_analyse$fid,
  iid = df_analyse$iid,
  zygosity = df_analyse$gm_zygosity,
  age = df_analyse$gm_age,
  bmi = df_analyse$gm_BMI
)

dat <- t(data.frame( # transposed outcome data 
  outcome = df_analyse$phage_alpha_div_shannon
))

# assign colnames of transposed outcome data as family id's 
colnames(dat) <- df_analyse$fid

dat <- as.data.frame(dat)

    # -------------------------------------------------------------------------- #    

id <- phen$fid

age <- phen$age

zyg <- factor(phen$zygosity)

bmi <- phen$bmi

################################################ data input ###################################################

res<-{}
cor_sat <- {}
model_tab <-{}

for (i in 1:nrow(dat)) {
  
  methy<-	as.numeric(dat[i, ])
  
  data<-	data.frame(methy, zyg, id, age, bmi)
  
  formula<- as.formula("methy ~ age + bmi")
  #############################################################################
  
  ace<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace"))
  ae<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae"))
  ce<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce"))
  
  e<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e"))
  
  ade <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ade"))
  sat <- summary(twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat"))
  
  ############################### P value for models ##############################################
  
  ace2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace")
  ae2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae")
  ce2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce")
  
  e2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e")
  
  ade2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ade")
  sat2<- twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat")
  
  #--------------------------------------------------------------------------------------------------------------------------#
  
   comparison_model<- sat2
  
  pval_ACE<- compare(ace2,sat2)
  Rpval_ACE<- pval_ACE$p.value
  
  pval_AE<- compare(ae2,comparison_model)
  Rpval_AE<- pval_AE$p.value
  
  pval_CE<- compare(ce2,comparison_model)
  Rpval_CE<- pval_CE$p.value
  
  pval_E<- compare(e2,comparison_model)
  Rpval_E<- pval_E$p.value
  
  pval_ADE<- compare(ade2,comparison_model)
  Rpval_ADE<- pval_ADE$p.value

  ##############################################################################################################
  #################################   generate comparison metrics    ###########################################
  ##############################################################################################################
  
  ### Saturated model ###
  
  sat_2LL<- -2*logLik(sat2)
  sat_AIC<- sat$AIC 
  sat_df<- attr(logLik(sat2), "df")
  
  ### ACE ###
  
  ACE_2LL<- -2*logLik(ace2)
  ACE_AIC<- ace$AIC 
  ACE_df<- attr(logLik(ace2), "df")
  ACE_DEL_2LL<- ACE_2LL - sat_2LL
  ACE_DEL_df<- pval_ACE$parameter[1]
  
  ACE_CHI<- pval_ACE$statistic[1]
  
  ### ADE ###
  
  ADE_2LL<- -2*logLik(ade2)
  ADE_AIC<- ade$AIC 
  ADE_df<- attr(logLik(ade2), "df")
  ADE_DEL_2LL<- ADE_2LL - sat_2LL
  ADE_DEL_df<- pval_ADE$parameter[1]
  ADE_CHI<- pval_ADE$statistic[1]
  
  ### AE ###
  
  AE_2LL<- -2*logLik(ae2)
  AE_AIC<- ae$AIC 
  AE_df<- attr(logLik(ae2), "df")
  AE_DEL_2LL<- AE_2LL - sat_2LL
  AE_DEL_df<- pval_AE$parameter[1]
  AE_CHI<- pval_AE$statistic[1]
  
  # ---- CE ---- 
  CE_2LL <- -2*logLik(ce2)
  CE_AIC <- ce$AIC
  CE_df <- attr(logLik(ce2),"df")
  CE_DEL_2LL <- CE_2LL - sat_2LL
  CE_DEL_df <- pval_CE$parameter[1]
  CE_CHI <- pval_CE$statistic[1]
  
  ### E ###
  
  E_2LL<- -2*logLik(e2)
  E_AIC<- e$AIC 
  E_df<- attr(logLik(e2), "df")
  E_DEL_2LL<- E_2LL - sat_2LL
  E_DEL_df<- pval_E$parameter[1]
  E_CHI<- pval_E$statistic[1]
  
  
  ################################# generate and format result table ###########################################
  
  
  res1<- rbind(
    c(row.names(dat)[i],"ACE",ace$acde[1,],ace$acde[2,],"-","-","-",ace$acde[3,],ace$corMZ,ace$corDZ,ace$zyg[1],ace$zyg[2],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
    c(row.names(dat)[i],"ADE",ade$acde[1,],"-","-","-",ade$acde[2,],ade$acde[3,],ade$corMZ,ade$corDZ,ade$zyg[1],ade$zyg[2],"ADE",ADE_2LL,ADE_df,ADE_AIC,ADE_DEL_2LL,ADE_DEL_df,ADE_CHI,Rpval_ADE),
    c(row.names(dat)[i],"AE", ae$acde[1,],"-","-","-","-","-","-",ae$acde[2,],ae$corMZ,ae$corDZ,ae$zyg[1],ae$zyg[2],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
    c(row.names(dat)[i],"E", "-","-","-","-","-","-","-","-","-",e$acde[1,],e$corMZ,e$corDZ,e$zyg[1],e$zyg[2],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E),
    c(row.names(dat)[i],"CE", "-","-","-",ce$acde[1,],"-","-","-",ce$acde[2,],ce$corMZ,ce$corDZ,ce$zyg[1],ce$zyg[2],"CE",CE_2LL,CE_df,CE_AIC,CE_DEL_2LL,CE_DEL_df,CE_CHI,Rpval_CE)
  )
  
  res<- rbind(res, res1)
  
  cor_sat1<- rbind(c(sat$corMZ,sat$corDZ))
  
  cor_sat<- rbind(cor_sat,c(row.names(dat)[i],cor_sat1))
  
  mod_tab1 <- rbind(
    c(row.names(dat)[i],"saturated",sat_2LL,sat_df,sat_AIC,"-","-","-","-"),
    c (row.names(dat)[i],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
    c(row.names(dat)[i],"ADE",ADE_2LL,ADE_df,ADE_AIC,ADE_DEL_2LL,ADE_DEL_df,ADE_CHI,Rpval_ADE),
    c(row.names(dat)[i],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
    c(row.names(dat)[i],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E)
  )
  
  model_tab<-rbind(model_tab,mod_tab1)
  
}

########################################## output result table ###############################################

head <- c("Trait",
         "Model","A","A_2.50%","A_97.50%","C","C_2.50","C_97.50%","D","D_2.50","D_97.50%","E","E_2.50%","E_97.50%","Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%","Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%","MZ:pairs/singletons","DZ:pairs/singletons","Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
resF <- rbind(head,res)


head_sat <- c("Trait",
              "Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%",
              "Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%")
res_sat <- rbind(head_sat,cor_sat)


head_mod <- c("Trait",
              "Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
res_mod <- rbind(head_mod,model_tab)

########################################## save result table ###############################################

write.table(resF,
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Shannon/Output_heritibility_phage_shannon_diversity.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_sat, 
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Shannon/Output_heritibility_phage_shannon_diversity_CORR_SAT.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_mod, 
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Shannon/Output_heritibility_phage_shannon_diversity_MOD.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

        # -------------------------------------------------------------------------- #  
  
            ####   RUN FOR TAXONOMIC RICHNESS  ---- 

# split data into phenotype data and outcome data 

phen <- data.frame(
  fid = df_analyse$fid,
  iid = df_analyse$iid,
  zygosity = df_analyse$gm_zygosity,
  age = df_analyse$gm_age,
  bmi = df_analyse$gm_BMI
)

dat <- t(data.frame( # transposed outcome data 
  outcome = df_analyse$phage_alpha_div_taxonomic_richness
))

# assign colnames of transposed outcome data as family id's 
colnames(dat) <- df_analyse$fid

dat <- as.data.frame(dat)

# -------------------------------------------------------------------------- #    

id <- phen$fid

age <- phen$age

zyg <- factor(phen$zygosity)

bmi <- phen$bmi


################################################ data input ###################################################

res<-{}
cor_sat <- {}
model_tab <-{}

for (i in 1:nrow(dat)) {
  
  methy<-	as.numeric(dat[i, ])
  
  data<-	data.frame(methy, zyg, id, age, bmi)
  
  formula<- as.formula("methy ~ age + bmi")
  #############################################################################
  
  ace<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace"))
  ae<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae"))
  ce<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce"))
  
  e<- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e"))
  
  ade <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ade"))
  sat <- summary(twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat"))
  
  ############################### P value for models ##############################################
  
  ace2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace")
  ae2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae")
  ce2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce")
  
  e2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e")
  
  ade2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ade")
  sat2<- twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat")
  
  #--------------------------------------------------------------------------------------------------------------------------#
  
  comparison_model<- sat2
  
  pval_ACE<- compare(ace2,sat2)
  Rpval_ACE<- pval_ACE$p.value
  
  pval_AE<- compare(ae2,comparison_model)
  Rpval_AE<- pval_AE$p.value
  
  pval_CE<- compare(ce2,comparison_model)
  Rpval_CE<- pval_CE$p.value
  
  pval_E<- compare(e2,comparison_model)
  Rpval_E<- pval_E$p.value
  
  pval_ADE<- compare(ade2,comparison_model)
  Rpval_ADE<- pval_ADE$p.value
  
  
  ##############################################################################################################
  #################################   generate comparison metrics    ###########################################
  ##############################################################################################################
  
  ### Saturated model ###
  
  sat_2LL<- -2*logLik(sat2)
  sat_AIC<- sat$AIC 
  sat_df<- attr(logLik(sat2), "df")
  
  ### ACE ###
  
  ACE_2LL<- -2*logLik(ace2)
  ACE_AIC<- ace$AIC 
  ACE_df<- attr(logLik(ace2), "df")
  ACE_DEL_2LL<- ACE_2LL - sat_2LL
  ACE_DEL_df<- pval_ACE$parameter[1]
  
  ACE_CHI<- pval_ACE$statistic[1]
  
  ### ADE ###
  
  ADE_2LL<- -2*logLik(ade2)
  ADE_AIC<- ade$AIC 
  ADE_df<- attr(logLik(ade2), "df")
  ADE_DEL_2LL<- ADE_2LL - sat_2LL
  ADE_DEL_df<- pval_ADE$parameter[1]
  ADE_CHI<- pval_ADE$statistic[1]
  
  ### AE ###
  
  AE_2LL<- -2*logLik(ae2)
  AE_AIC<- ae$AIC 
  AE_df<- attr(logLik(ae2), "df")
  AE_DEL_2LL<- AE_2LL - sat_2LL
  AE_DEL_df<- pval_AE$parameter[1]
  AE_CHI<- pval_AE$statistic[1]
  
  # ---- CE ---- 
  CE_2LL <- -2*logLik(ce2)
  CE_AIC <- ce$AIC
  CE_df <- attr(logLik(ce2),"df")
  CE_DEL_2LL <- CE_2LL - sat_2LL
  CE_DEL_df <- pval_CE$parameter[1]
  CE_CHI <- pval_CE$statistic[1]
  
  ### E ###
  
  E_2LL<- -2*logLik(e2)
  E_AIC<- e$AIC 
  E_df<- attr(logLik(e2), "df")
  E_DEL_2LL<- E_2LL - sat_2LL
  E_DEL_df<- pval_E$parameter[1]
  E_CHI<- pval_E$statistic[1]
  
  ################################# generate and format result table ###########################################
  
  
  res1<- rbind(
    c(row.names(dat)[i],"ACE",ace$acde[1,],ace$acde[2,],"-","-","-",ace$acde[3,],ace$corMZ,ace$corDZ,ace$zyg[1],ace$zyg[2],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
    c(row.names(dat)[i],"ADE",ade$acde[1,],"-","-","-",ade$acde[2,],ade$acde[3,],ade$corMZ,ade$corDZ,ade$zyg[1],ade$zyg[2],"ADE",ADE_2LL,ADE_df,ADE_AIC,ADE_DEL_2LL,ADE_DEL_df,ADE_CHI,Rpval_ADE),
    c(row.names(dat)[i],"AE", ae$acde[1,],"-","-","-","-","-","-",ae$acde[2,],ae$corMZ,ae$corDZ,ae$zyg[1],ae$zyg[2],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
    c(row.names(dat)[i],"E", "-","-","-","-","-","-","-","-","-",e$acde[1,],e$corMZ,e$corDZ,e$zyg[1],e$zyg[2],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E),
    c(row.names(dat)[i],"CE", "-","-","-",ce$acde[1,],"-","-","-",ce$acde[2,],ce$corMZ,ce$corDZ,ce$zyg[1],ce$zyg[2],"CE",CE_2LL,CE_df,CE_AIC,CE_DEL_2LL,CE_DEL_df,CE_CHI,Rpval_CE)
  )
  
  
  
  res<- rbind(res, res1)
  
  cor_sat1<- rbind(c(sat$corMZ,sat$corDZ))
  
  cor_sat<- rbind(cor_sat,c(row.names(dat)[i],cor_sat1))
  
  mod_tab1 <- rbind(
    c(row.names(dat)[i],"saturated",sat_2LL,sat_df,sat_AIC,"-","-","-","-"),
    c (row.names(dat)[i],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
    c(row.names(dat)[i],"ADE",ADE_2LL,ADE_df,ADE_AIC,ADE_DEL_2LL,ADE_DEL_df,ADE_CHI,Rpval_ADE),
    c(row.names(dat)[i],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
    c(row.names(dat)[i],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E)
  )
  
  model_tab<-rbind(model_tab,mod_tab1)
  
  
}

########################################## output result table ###############################################

head <- c("Trait",
          "Model","A","A_2.50%","A_97.50%","C","C_2.50","C_97.50%","D","D_2.50","D_97.50%","E","E_2.50%","E_97.50%","Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%","Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%","MZ:pairs/singletons","DZ:pairs/singletons","Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
resF <- rbind(head,res)


head_sat <- c("Trait",
              "Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%",
              "Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%")
res_sat <- rbind(head_sat,cor_sat)


head_mod <- c("Trait",
              "Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
res_mod <- rbind(head_mod,model_tab)

########################################## save result table ###############################################

write.table(resF,
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Taxonomic_richness/Output_heritibility_phage_taxonomic_richness.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_sat, 
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Taxonomic_richness/Output_heritibility_phage_taxonomic_richness_CORR_SAT.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_mod, 
            file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Database_independent_analysis/Analysis/Phages_heritability/Alpha_diversity/Taxonomic_richness/Output_heritibility_phage_taxonomic_richness_MOD.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


# -------------------------------------------------------------------------- #  
