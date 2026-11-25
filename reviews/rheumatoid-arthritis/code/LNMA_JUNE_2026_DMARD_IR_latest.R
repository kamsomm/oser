######################################################################################
######################################################################################
##BY KAMSO MOHAMMED MUJAAB
##THIS IS AN R-SCRIPT FOR RUNNING THE ANALYSIS OF THE LNMA FOR DMARD-BIOLOGIC
##
##
##
######################################################################################
######################################################################################
#################################


#setwd("/Users/mujaab/Documents/OneDrive - University of Calgary/DMARD-Biologic")
#setwd("/Users/mujaab/Documents/OneDrive - University of Calgary/BIO-IR:TNF-IR")
oldwd<-getwd()
setwd('/Users/mujaab/Documents/OneDrive - University of Calgary/DMARD-IR')
source("All functions for LNMA_2026.R")


#dmard_ir_new<-read.csv("review_60381_extracted_data_csv_20260618112552.csv", header = T)
#dmard_ir_new<-read.csv("review_60381_extracted_data_csv_20260625121704.csv", header = T)
#dmard_ir_new<-read.csv("review_60381_extracted_data_csv_20260627112741.csv", header = T)

dmard_ir_new<-read.csv("review_60381_extracted_data_csv_20260901033527.csv", header = T)

View(dmard_ir_new)
unique(dmard_ir_new$Intervention)
dim(dmard_ir_new)



row.names(dmard_ir_new)<-NULL

# #Solving wrongly written intervention
# dmard_biologic_new[dmard_biologic_new$Study.Identifier=="Gottenberg 2016",2][1:2]<-c("Non-TNF Biologic","Second Anti-TNF")

length(unique(dmard_ir_new$Study.Identifier))
# 
# dims2<-matrix(1:nrow(dmard_ir_new),10,nrow(dmard_ir_new)%/%10)
# ncol(dims2)
# 
# 
# dmard_list_new_ir<-list()
# S_ID_new2<-c()
# 
# for(i in 1:ncol(dims2)){
#   if(any(S_ID_new2 == unique(dmard_ir_new[dims2[,i],1]))){
#     dmard_ir_new[dims2[,i],1]<-paste0(dmard_ir_new[dims2[,i],1],"_",i)
#     dmard_list_new_ir[[i]]<-cbind(Study_no= rep(i,nrow(dmard_ir_new[dims2[which(substr(dmard_ir_new[dims2[,i],2],1,3)!="Arm"),i],]))
#                                         ,dmard_ir_new[dims2[which(substr(dmard_ir_new[dims2[,i],2],1,3)!="Arm"),i],])
#   }else{
#     S_ID_new2<-c(S_ID_new2,unique(dmard_ir_new[dims2[,i],1]))
#     dmard_list_new_ir[[i]]<-cbind(Study_no= rep(i,nrow(dmard_ir_new[dims2[which(substr(dmard_ir_new[dims2[,i],2],1,3)!="Arm"),i],]))
#                                         ,dmard_ir_new[dims2[which(substr(dmard_ir_new[dims2[,i],2],1,3)!="Arm"),i],])
#   }
#   
# }

sub_dmard_new_ir <- dmard_ir_new %>%
  filter(!grepl("^Arm\\s*\\d+$", Intervention))


# Rename duplicate study identifiers with .a, .b, .c suffixes
# Strategy: detect when the same study name appears in separate blocks
# (non-consecutive rows), then label subsequent blocks with .a, .b, etc.

library(dplyr)

sub_dmard_new_ir <- sub_dmard_new_ir %>%
  mutate(
    .block = cumsum(Study.Identifier != lag(Study.Identifier, default = first(Study.Identifier)))
  ) %>%
  group_by(Study.Identifier) %>%
  mutate(
    .study_block_rank = dense_rank(.block),
    .n_blocks = n_distinct(.block)
  ) %>%
  ungroup() %>%
  mutate(
    Study.Identifier = if_else(
      .n_blocks > 1,
      paste0(Study.Identifier, ".", letters[.study_block_rank]),
      Study.Identifier
    )
  ) %>%
  select(-.block, -.study_block_rank, -.n_blocks)%>%
  mutate(Study_no = match(Study.Identifier, unique(Study.Identifier))) %>%
  relocate(Study_no, .before = 1)




length(unique(sub_dmard_new_ir$Study.Identifier))
unique(sub_dmard_new_ir$Study_no)
dmard_ir_read_new<-sub_dmard_new_ir



#ORIGINAL FIRST DMARD NAIVE CLEAN
write.csv(dmard_ir_read_new,"First dmard IR clean JUNE 2026_V1.csv",row.names = F)


intervention_dmardir_list<-read_xlsx("first trial master intervention dmard-IR list new JUNE 2026_V4.xlsx")
#intervention_dmardir_list<-intervention_dmardir_list[,-9]
dmard_ir_read_new_update3_included2<-cbind(dmard_ir_read_new,intervention_dmardir_list%>% select(DoseLevel01,Unique,Intervention2,Intervention3,Intervention4))
#dmard_ir_read_new_update3_included2_main<-dmard_ir_read_new_update3_included2[dmard_ir_read_new_update3_included2$Main.record..Y.N.=="Y",]
dmard_ir_read_new_update3_included2_main<-dmard_ir_read_new_update3_included2
creating_master_intervention_dmard_ir<-dmard_ir_read_new_update3_included2_main[,c("Study.Identifier","Unique","Intervention4","DoseLevel01")]

#Standardize intervention names
for(i in 1:length(unique(creating_master_intervention_dmard_ir$Unique))){
  eeee<-creating_master_intervention_dmard_ir[creating_master_intervention_dmard_ir$Unique==i,3][1]
  creating_master_intervention_dmard_ir[which(creating_master_intervention_dmard_ir$Unique==i),3]<-eeee 
}

##adding variable for dropping similar interventions and low dose levels
dmard_ir<-dmard_ir_read_new_update3_included2_main
dmard_ir1<-dmard_ir[dmard_ir$DoseLevel01==1,]

#Renumber studies
dmard_ir1$Study_no<-renumber(dmard_ir1,"Study_no")
row.names(dmard_ir1)<-NULL
View(dmard_ir1)

####Creating dataset for Indirectness assessment
Dmard_ir_dataset_for_indirectness_assessment<-dmard_ir1[,c("Study_no","Study.Identifier","Intervention","Year.of.publication","NMA.population",
                                                           "X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.",
                                                           "Baseline.Mean.age..years.","Baseline.Mean.disease.duration..years.","Baseline.Mean.MTX.dose..mg.",
                                                           "Baseline.Proportion.female..0.1.","DAS.28.Baseline.Mean","HAQ.DI.Baseline.Mean" ,
                                                           "Baseline.Daily.prednisone.equivalent.dose.actual..mg.","Baseline.Prednisone.dose.max..mg.",
                                                           "Baseline.Proportion.maintenance.of.oral.steroid..0.1.",
                                                           "Baseline.Proportion.prior.non.TNF.biologic.use..0.1.","Baseline.Proportion.prior.TNF.use..0.1.")]

Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-gsub(" weeks","",Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.)
Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-ifelse(Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.=="6 months",Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.,round(as.numeric(Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.)/4.345,0)) 
Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-as.numeric(gsub(" months","",Dmard_ir_dataset_for_indirectness_assessment$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.))


 
colnames(Dmard_ir_dataset_for_indirectness_assessment)<-c("Study_no","Study Identifier","Intervention","Year of publication", "NMA Population",
                                                          "Trial duration in months","Baseline mean age (years)","Baseline mean disease duration (years)","MTX Dose (mg/wk)",
                                                          "Baseline proportion female (0,1)","DAS28 baseline (mean)","HAQ-DI baseline (mean)",
                                                          "Baseline daily prednisone equivalent dose actual (mg)",
                                                          "Baseline prednisone dose max (mg)","Baseline proportion maintenance of oral steroid (0,1)",
                                                          "Baseline proportion prior non-TNF biologic use (0,1)","Baseline proportion prior TNF use (0,1)")
write.csv(Dmard_ir_dataset_for_indirectness_assessment,"Indirectness data for Dmard IR dataset criteria_JUNE2026.csv",row.names = F)



#Changing some intervention names to have no dosages 
dmard_ir1$Intervention<-creating_master_intervention_dmard_ir[creating_master_intervention_dmard_ir$DoseLevel01==1,3]

#####Creating dmard_ir ACR50 dataset
#Formatting disease duration to months
final_dmard_irDataset<-dmard_ir1
row.names(final_dmard_irDataset)<-NULL
final_dmard_irDataset[final_dmard_irDataset$Study.Identifier=="Keystone 2009","X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy."]<-24
final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-gsub(" weeks","",final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.)
final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-ifelse(final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.=="6 months",final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.,round(as.numeric(final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.)/4.345,0)) 
final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.<-as.numeric(gsub(" months","",final_dmard_irDataset$X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.))
#colnames(final_biologicDataset)[23]<-"Trial.duration..in.months"
final_dmard_irDataset<- final_dmard_irDataset%>%
  dplyr::rename("Trial.duration..in.months" = X4..Trial.duration..in.weeks..end.of.trial...end.of.randomized.treatment.period..e.g..switch.to.open.label..or..for.strategy.trials..time.of.the.first.switch.in.therapy.)
View(final_dmard_irDataset)

write.csv(final_dmard_irDataset,"Final DMARD-IR dataset V1_26JUNE2026.csv",row.names = F)


#Creating medtable
master.medtable.temp_dmard_ir <-cbind(as.numeric(unique(creating_master_intervention_dmard_ir$Unique)),
                                      unique(creating_master_intervention_dmard_ir$Intervention4))

master.medtable.temp_dmard_ir2<-as.data.frame(master.medtable.temp_dmard_ir)
master.medtable.temp_dmard_ir2$V1<-as.numeric(master.medtable.temp_dmard_ir2$V1)
master.medtable.temp_dmard_ir2<-master.medtable.temp_dmard_ir2[order(master.medtable.temp_dmard_ir2$V1),]
colnames(master.medtable.temp_dmard_ir2) <- c('mednum','medname')
master.medtable_dmard_ir<- master.medtable.temp_dmard_ir2
row.names(master.medtable_dmard_ir)<-NULL

#Inspect the table
master.medtable_dmard_ir
View(master.medtable_dmard_ir)
write.csv(master.medtable_dmard_ir,"Master medtable for dmardIR_JUNE2026.csv",row.names = F)
View(final_dmard_irDataset)

####
table(final_dmard_irDataset$NMA.population)
# ##DMARD-IR review
# final_dmard_irDataset_actual<- final_dmard_irDataset %>%
#   filter(NMA.population == "DMARD-IR")
# 
# final_dmard_irDataset_actual1 <- final_dmard_irDataset_actual %>%
#   mutate(
#     `ACR.50.End.of.trial.n` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
#         NA,
#         `ACR.50.End.of.trial.n`
#       ),
#     
#     `ACR.50.End.of.trial.N` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
#         NA,
#         `ACR.50.End.of.trial.N`
#       ),
#     
#     `ACR.50.12.months.n` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,6,23,24),
#         NA,
#         `ACR.50.12.months.n`
#       ),
#     
#     `ACR.50.12.months.N` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,6,23,24),
#         NA,
#         `ACR.50.12.months.N`
#       ),
#     
#     `ACR.50.6.months.n` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24),
#         NA,
#         `ACR.50.6.months.n`
#       ),
#     
#     `ACR.50.6.months.N` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24),
#         NA,
#         `ACR.50.6.months.N`
#       ),
#     
#     `ACR.50.3.months.n` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24),
#         NA,
#         `ACR.50.3.months.n`
#       ),
#     `ACR.50.3.months.N` =
#       ifelse(
#         `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24),
#         NA,
#         `ACR.50.3.months.N`
#       )
#   )
# 
# df_acr50 <- final_dmard_irDataset_actual1%>% 
#   select(matches("^ACR\\.50\\..*\\.(n|N)$"))%>%
#   select(-c("ACR.50.4.months.n", "ACR.50.4.months.N" ,"ACR.50.End.of.trial.22.Weeks..n" ,"ACR.50.End.of.trial.22.Weeks..N"))
#   
# 
# ACR_dataset_dmardir<-cbind(final_dmard_irDataset_actual1[,c(1:27,30:37,40:43,47:57)],df_acr50, final_dmard_irDataset_actual1[,c("DoseLevel01")])
# 
# ACR_dataset_dmardir<- ACR_dataset_dmardir %>%
#   dplyr:: rename("DoseLevel01" =`final_dmard_irDataset_actual1[, c("DoseLevel01")]`)
# View(ACR_dataset_dmardir)
# 
# #Dropping empty columns of ACR
# ACR_dataset_dmardir2<-ACR_dataset_dmardir
# 
# ##ACR50
# ACR50_dataset_dmardir<-ACR_dataset_dmardir2[,c(1:3,51:58)]
# ACR50_dataset_dmardir<-as.data.frame(ACR50_dataset_dmardir)
# View(ACR50_dataset_dmardir)
# 
# ACR50_datalist_dmardir <-list()
# 
# n_study1<-unique(ACR50_dataset_dmardir$Study_no)
# 
# for(i in 1:length(n_study1) ){
#   eee<-ACR50_dataset_dmardir[ACR50_dataset_dmardir$Study_no==n_study1[i],]
#   eee[eee=="NA"]<-NA
#   ggg<-eee[colSums(!is.na(eee)) > 0]
#   if(ncol(ggg)==11){
#     ggg<-ggg[,-c(4:9)]
#     ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
#     colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
#     ACR50_datalist_dmardir[[i]]<-ggg
#   } else if(ncol(ggg)==9){
#     ggg<-ggg[,-c(4:7)]
#     ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
#     colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
#     ACR50_datalist_dmardir[[i]]<-ggg
#   }else if(ncol(ggg)==7){
#     ggg<-ggg[,-c(4:5)]
#     ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
#     colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
#     ACR50_datalist_dmardir[[i]]<-ggg
#   }else if(ncol(ggg)==5){
#     ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
#     colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
#     ACR50_datalist_dmardir[[i]]<-ggg
#   }else if((ncol(ggg)<=4) || (ncol(ggg)==6) || (ncol(ggg)==8) || (ncol(ggg)==10)){
#     ggg<-ggg[,c(1,2,3)]
#     aaaa<-matrix(NA,nrow(ggg),3)
#     ggg<-cbind(ggg,aaaa)
#     colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
#     ACR50_datalist_dmardir[[i]]<-ggg
#   }
#   
#   
# }
# 
# 
# 
# 
# ACR50_datalist_dmardir
# ACR50_big_data_dmardir= do.call(rbind, ACR50_datalist_dmardir)
# ACR50_big_data_dmardir
# ACR50_big_data_dmardir<-cbind(ACR50_big_data_dmardir,DoseLevel01=ACR_dataset_dmardir2$DoseLevel01,TrialDuration=ACR_dataset_dmardir2$Trial.duration..in.months)
# View(ACR50_big_data_dmardir)
# dummy <- NULL
# dummy[is.na(ACR50_big_data_dmardir$month)] = 0
# dummy[ACR50_big_data_dmardir$month=="ACR.50.3.months.n"] = 1
# #dummy[ACR50_big_data_biologic$month=="ACR.50.4.months.n"] = 2
# dummy[ACR50_big_data_dmardir$month=="ACR.50.6.months.n"] = 2
# dummy[ACR50_big_data_dmardir$month=="ACR.50.12.months.n"] = 3
# dummy[ACR50_big_data_dmardir$month=="ACR.50.End.of.trial.n"] = 4
# ACR50_big_data_dmardir<-cbind(ACR50_big_data_dmardir,dummy)
# 
# ACR50_big_data_dmardir<- ACR50_big_data_dmardir%>%
#   dplyr::rename("Time_point" = dummy)
# ACR50_big_data_dmardir
# 
# ##getting the actual contributed studies for acr50
# ACR50_big_data_dmardir2<-ACR50_big_data_dmardir[ACR50_big_data_dmardir$Time_point!=0,]
# ACR50_big_data_dmardir2$month <- sub("^ACR\\.50\\.(.*)\\.(n|N)$", "\\1", ACR50_big_data_dmardir2$month)
# 
# ##Solving error in Blanco 2017 placebo arm
# #ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Blanco 2017",][1,c(4,5)]<-c(15,138)#responders not reported for placebo arm acr50
# 
# # ##Solving error in Burnmester 2013 placebo arm
# # ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",][1,c(4,5)]<-c(11,131)
# # ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",][2,c(4,5)]<-c(35,132)
# # ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",]$month<-"3.months"
# 
# ##Drop furst study due to high ROB
# ACR50_big_data_dmardir2<- ACR50_big_data_dmardir2 %>%
#   filter(!study_name %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
#                        "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
#                        "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023"))
# 
# #ACR50_big_data_dmardir2<-na.omit(ACR50_big_data_dmardir2)
# 
# 
# 
# ACR50_big_data_dmardir2$study<-renumber(ACR50_big_data_dmardir2,"study")
# 
# ##Moreland 2022 graphic
# ACR50_big_data_dmardir2[ACR50_big_data_dmardir2$study_name=="Moreland 2002","responders"]<-paste0("8888",round(c(0.063482*32,0.158089*32),0))
# ##Fleischmann 2023.b, placebo arm missing, created a fictitious data
# ACR50_big_data_dmardir2[ACR50_big_data_dmardir2$study_name=="Fleischmann 2023.b",c("responders","sampleSize")][1,]<- c(round(0.12*270,0),270)
# 
# #Drop NA rows
# ACR50_big_data_dmardir2<- ACR50_big_data_dmardir2%>%
#   filter(!is.na(responders))
# 
# row.names(ACR50_big_data_dmardir2)<-NULL
# View(ACR50_big_data_dmardir2)
# 
# library(dplyr)
# for(i in 1:nrow(ACR50_big_data_dmardir2)){
#   ppp<- master.medtable_dmard_ir %>% filter(ACR50_big_data_dmardir2$treatment[i] ==master.medtable_dmard_ir$medname )
#   ACR50_big_data_dmardir2$treatment[i]<-ppp[1,1]
# }
# 
# 
# ACR50_big_data_dmardir2$responders<-as.numeric(gsub(8888,'',ACR50_big_data_dmardir2$responders))
# 
# 
# #####Formatting the month column
# ACR50_big_data_dmardir2$month<-sub(".months","",ACR50_big_data_dmardir2$month)
# ACR50_big_data_dmardir2$month<-as.numeric(ifelse(ACR50_big_data_dmardir2$month=="End.of.trial",ACR50_big_data_dmardir2$TrialDuration,ACR50_big_data_dmardir2$month))
# View(ACR50_big_data_dmardir2)
# 
# ##Creating initial legend
# master.medtable<-master.medtable_dmard_ir
# dmardIR.medlegend <- new.medtable(ACR50_big_data_dmardir2,'treatment')
# dmardIR.medlegend
# 
# ##NMA plot reason for narrowing of network
# nma.plot(ACR50_big_data_dmardir2[,1],ACR50_big_data_dmardir2[,3],data = ACR50_big_data_dmardir2, trtname = as.vector(dmardIR.medlegend$medname))
# 
# ################################################################################################################
# #trying the narrowing a network function
# #############################################################################################
# ACR50_big_data_dmardir2
# 
# #Narrow 1
# dmardIR_narrow<-narrow.all(ACR50_big_data_dmardir2)
# dmardIR_narrow
# 
# # #Narrow 2
# dmardIR_narrow2<-narrow.all(dmardIR_narrow[[2]])
# dmardIR_narrow2
# 
# 
# #sub network plot 1
# dmardIR.d<-dmardIR_narrow[[1]]
# length(unique(dmardIR.d$treatment))
# dmardIR_subl<-new.medtable(dmardIR.d,'treatment')
# # ##Removing the Yamamoto 2014.b study
# # dmardIR.d<- dmardIR.d %>%
# #   filter(!study_name=="Yamamoto 2014.b")
# 
# nma.plot(dmardIR.d[,1],dmardIR.d[,3],data = dmardIR.d, trtname = as.vector(dmardIR_subl$medname), title = "ACR50 Network plot")
# 
# 
# #sub network plot 2
# dmardIR.d2<-dmardIR_narrow2[[1]]
# length(unique(dmardIR.d2$treatment))
# dmardIR_subl2<-new.medtable(dmardIR.d2,'treatment')
# nma.plot(dmardIR.d2[,1],dmardIR.d2[,3],data = dmardIR.d2, trtname = as.vector(dmardIR_subl2$medname), title = "ACR50 Network plot")
# 
# ##############################################################################
# #Performing sub network plot 1
# library(rjags)
# library(coda)
# library(gemtc)
# View(dmardIR.d)
# dmardIR_subl
# 
# # network_data_ACR50_dmardIR_1<-dmardIR.d%>%
# #   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
# #            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
# network_data_ACR50_dmardIR_1<-dmardIR.d%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month"))
# 
# network_data_ACR50_dmardIR_1$treatment<-renumber(network_data_ACR50_dmardIR_1,"treatment")
# ##To download acr50 data
# network_data_ACR50_dmardIR_data_download<-network_data_ACR50_dmardIR_1
# network_data_ACR50_dmardIR_data_download$treatment<-selected_treatment_name(network_data_ACR50_dmardIR_data_download$treatment,dmardIR_subl)
# write.csv(network_data_ACR50_dmardIR_data_download,"ACR50 DMARD-IR dataset_JUNE2026.csv",row.names = F)
# # 1. Open jpeg file
# jpeg("ACR50 -DMARD-IR Network plot 1 with named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_dmardIR_1$study,network_data_ACR50_dmardIR_1$treatment,data = network_data_ACR50_dmardIR_1, trtname = as.vector(dmardIR_subl$medname),title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# # 1. Open jpeg file
# jpeg("ACR50 -DMARD-IR Network plot 1 without named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_dmardIR_1$study,network_data_ACR50_dmardIR_1$treatment,data = network_data_ACR50_dmardIR_1,title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# 
# network_data_ACR50_dmardIR_1<-mtc.network(network_data_ACR50_dmardIR_1)
# 
# 
# ##RUN the model
# #mtc.hy.empirical.lor("semi-objective", "non-pharma")
# #mtc.hy.prior("std.dev", "dunif", 0, 2)
# #mtc.hy.prior("prec", "dgamma", 0.01, 0.01)
# set.seed(501)
# network_data_ACR50_dmardIR_1.model <- mtc.model(network_data_ACR50_dmardIR_1, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
# 
# # By default, the model will have 4 chains - generate a seed for each
# seeds <- sample.int(4, n = .Machine$integer.max)
# # Apply JAGS RNG settings to each chain
# network_data_ACR50_dmardIR_1.model$inits<- mapply(c, network_data_ACR50_dmardIR_1.model$inits, list(
#   list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
#   list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
#   list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
#   list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
# results_dmardIR<-mtc.run(network_data_ACR50_dmardIR_1.model,n.iter = 50000, n.adapt = 30000)
# summary(results_dmardIR)
# 
# View(relative.effect.table.OR(results_dmardIR))
# refOR_dmardIR<-as.data.frame(relative.effect.table.OR2(results_dmardIR, comparator = "9"))
# 
# #library for mapvalues
# library(plyr)
# row.names(refOR_dmardIR)<-mapvalues(row.names(refOR_dmardIR),
#                                 from = c(1:(nrow(dmardIR_subl))),
#                                 to = as.vector(dmardIR_subl$medname))
# colnames(refOR_dmardIR)<-mapvalues(colnames(refOR_dmardIR),
#                                from = c(1:(nrow(dmardIR_subl))),
#                                to = as.vector(dmardIR_subl$medname))
# View(refOR_dmardIR)
# # which(row(refOR)==col(refOR))
# # refOR[row(refOR)==col(refOR)]
# index_diagonal2<-data.frame(index1=1:nrow(refOR_dmardIR),index2=1:nrow(refOR_dmardIR))
# options(warn=-1)
# for(i in 1:nrow(refOR_dmardIR)){
#   refOR_dmardIR[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_dmardIR[row(refOR_dmardIR)==col(refOR_dmardIR)])[i],dmardIR_subl)
# }
# options(warn=0)
# View(refOR_dmardIR)
# 
# write.csv(refOR_dmardIR,"ACR50 -DMARD-IR Network 1 pairwise comparisons_JUNE2026.csv", row.names = T)
# # library(openxlsx)
# # write.xlsx(refOR_dmardIR,"Results for DMARD-IR_JUNE2026.xlsx",sheetName = "ACR50",rowNames = T)
# 
# 
# ##############################################################################
# #Performing sub network plot 2
# library(rjags)
# library(coda)
# library(gemtc)
# View(dmardIR.d2)
# dmardIR_subl2
# 
# # network_data_ACR50_dmardIR_1<-dmardIR.d%>%
# #   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
# #            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
# network_data_ACR50_dmardIR_2<-dmardIR.d2%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month"))
# 
# network_data_ACR50_dmardIR_2$treatment<-renumber(network_data_ACR50_dmardIR_2,"treatment")
# ##To download acr50 data
# network_data_ACR50_dmardIR_data_download2<-network_data_ACR50_dmardIR_2
# network_data_ACR50_dmardIR_data_download2$treatment<-selected_treatment_name(network_data_ACR50_dmardIR_data_download2$treatment,dmardIR_subl2)
# write.csv(network_data_ACR50_dmardIR_data_download2,"ACR50 DMARD-IR_subnetwork dataset_JUNE2026.csv",row.names = F)
# # 1. Open jpeg file
# jpeg("ACR50 Network plot (subnetwork) with named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_dmardIR_2$study,network_data_ACR50_dmardIR_2$treatment,data = network_data_ACR50_dmardIR_2, trtname = as.vector(dmardIR_subl2$medname),title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# # 1. Open jpeg file
# jpeg("ACR50 Network plot (subnetwork) without named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_dmardIR_2$study,network_data_ACR50_dmardIR_2$treatment,data = network_data_ACR50_dmardIR_1,title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# 
# network_data_ACR50_dmardIR_2<-mtc.network(network_data_ACR50_dmardIR_2)
# 
# 
# ##RUN the model
# #mtc.hy.empirical.lor("semi-objective", "non-pharma")
# #mtc.hy.prior("std.dev", "dunif", 0, 2)
# #mtc.hy.prior("prec", "dgamma", 0.01, 0.01)
# set.seed(502)
# network_data_ACR50_dmardIR_2.model <- mtc.model(network_data_ACR50_dmardIR_2, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
# 
# # By default, the model will have 4 chains - generate a seed for each
# seeds <- sample.int(4, n = .Machine$integer.max)
# # Apply JAGS RNG settings to each chain
# network_data_ACR50_dmardIR_2.model$inits<- mapply(c, network_data_ACR50_dmardIR_2.model$inits, list(
#   list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
#   list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
#   list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
#   list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
# results_dmardIR2<-mtc.run(network_data_ACR50_dmardIR_2.model,n.iter = 50000, n.adapt = 30000)
# summary(results_dmardIR2)
# 
# View(relative.effect.table.OR(results_dmardIR2))
# refOR_dmardIR2<-as.data.frame(relative.effect.table.OR(results_dmardIR2))
# 
# #library for mapvalues
# library(plyr)
# row.names(refOR_dmardIR2)<-mapvalues(row.names(refOR_dmardIR2),
#                                     from = c(1:(nrow(dmardIR_subl2))),
#                                     to = as.vector(dmardIR_subl2$medname))
# colnames(refOR_dmardIR2)<-mapvalues(colnames(refOR_dmardIR2),
#                                    from = c(1:(nrow(dmardIR_subl2))),
#                                    to = as.vector(dmardIR_subl2$medname))
# View(refOR_dmardIR2)
# # which(row(refOR)==col(refOR))
# # refOR[row(refOR)==col(refOR)]
# index_diagonal2<-data.frame(index1=1:nrow(refOR_dmardIR2),index2=1:nrow(refOR_dmardIR2))
# options(warn=-1)
# for(i in 1:nrow(refOR_dmardIR2)){
#   refOR_dmardIR2[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_dmardIR2[row(refOR_dmardIR2)==col(refOR_dmardIR2)])[i],dmardIR_subl2)
# }
# options(warn=0)
# View(refOR_dmardIR2)
# 
# write.csv(refOR_dmardIR2,"ACR50 (subnetwork) pairwise comparisons_JUNE2026.csv", row.names = T)


#MTX-IR review
final_mtx_irDataset_actual<- final_dmard_irDataset %>%
  filter(NMA.population == "MTX-IR")

final_mtx_irDataset_actual1 <- final_mtx_irDataset_actual %>%
  mutate(
    `ACR.50.End.of.trial.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
        NA,
        `ACR.50.End.of.trial.n`
      ),
    
    `ACR.50.End.of.trial.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
        NA,
        `ACR.50.End.of.trial.N`
      ),
    
    `ACR.50.12.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,23,24),
        NA,
        `ACR.50.12.months.n`
      ),
    
    `ACR.50.12.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,23,24),
        NA,
        `ACR.50.12.months.N`
      ),
    
    `ACR.50.6.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24),
        NA,
        `ACR.50.6.months.n`
      ),
    
    `ACR.50.6.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24),
        NA,
        `ACR.50.6.months.N`
      ),
    
    `ACR.50.3.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24),
        NA,
        `ACR.50.3.months.n`
      ),
    `ACR.50.3.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24),
        NA,
        `ACR.50.3.months.N`
      )
  )

df_acr502 <- final_mtx_irDataset_actual1%>% 
  select(matches("^ACR\\.50\\..*\\.(n|N)$"))%>%
  select(-c("ACR.50.4.months.n", "ACR.50.4.months.N" ,"ACR.50.End.of.trial.22.Weeks..n" ,"ACR.50.End.of.trial.22.Weeks..N"))


ACR_dataset_mtxir<-cbind(final_mtx_irDataset_actual1[,c(1:27,30:37,40:43,47:57)],df_acr502, final_mtx_irDataset_actual1[,c("DoseLevel01")])

ACR_dataset_mtxir<- ACR_dataset_mtxir %>%
  dplyr:: rename("DoseLevel01" =`final_mtx_irDataset_actual1[, c("DoseLevel01")]`)
View(ACR_dataset_mtxir)

#Dropping empty columns of ACR
ACR_dataset_mtxir2<-ACR_dataset_mtxir

##ACR50
ACR50_dataset_mtxir<-ACR_dataset_mtxir2[,c(1:3,51:58)]
ACR50_dataset_mtxir<-as.data.frame(ACR50_dataset_mtxir)
View(ACR50_dataset_mtxir)

ACR50_datalist_mtxir <-list()
n_study2<-unique(ACR50_dataset_mtxir$Study_no)

for(i in 1:length(n_study2) ){
  eee<-ACR50_dataset_mtxir[ACR50_dataset_mtxir$Study_no==n_study2[i],]
  eee[eee=="NA"]<-NA
  ggg<-eee[colSums(!is.na(eee)) > 0]
  if(ncol(ggg)==11){
    ggg<-ggg[,-c(4:9)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    ACR50_datalist_mtxir[[i]]<-ggg
  } else if(ncol(ggg)==9){
    ggg<-ggg[,-c(4:7)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    ACR50_datalist_mtxir[[i]]<-ggg
  }else if(ncol(ggg)==7){
    ggg<-ggg[,-c(4:5)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    ACR50_datalist_mtxir[[i]]<-ggg
  }else if(ncol(ggg)==5){
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    ACR50_datalist_mtxir[[i]]<-ggg
  }else if((ncol(ggg)<=4) || (ncol(ggg)==6) || (ncol(ggg)==8) || (ncol(ggg)==10)){
    ggg<-ggg[,c(1,2,3)]
    aaaa<-matrix(NA,nrow(ggg),3)
    ggg<-cbind(ggg,aaaa)
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    ACR50_datalist_mtxir[[i]]<-ggg
  }
  
  
}




ACR50_datalist_mtxir
ACR50_big_data_mtxir= do.call(rbind, ACR50_datalist_mtxir)
ACR50_big_data_mtxir
ACR50_big_data_mtxir<-cbind(ACR50_big_data_mtxir,DoseLevel01=ACR_dataset_mtxir2$DoseLevel01,TrialDuration=ACR_dataset_mtxir2$Trial.duration..in.months)
View(ACR50_big_data_mtxir)
dummy <- NULL
dummy[is.na(ACR50_big_data_mtxir$month)] = 0
dummy[ACR50_big_data_mtxir$month=="ACR.50.3.months.n"] = 1
dummy[ACR50_big_data_mtxir$month=="ACR.50.6.months.n"] = 2
dummy[ACR50_big_data_mtxir$month=="ACR.50.12.months.n"] = 3
dummy[ACR50_big_data_mtxir$month=="ACR.50.End.of.trial.n"] = 4
ACR50_big_data_mtxir<-cbind(ACR50_big_data_mtxir,dummy)

ACR50_big_data_mtxir<- ACR50_big_data_mtxir%>%
  dplyr::rename("Time_point" = dummy)
ACR50_big_data_mtxir

##getting the actual contributed studies for acr50
ACR50_big_data_mtxir2<-ACR50_big_data_mtxir[ACR50_big_data_mtxir$Time_point!=0,]
ACR50_big_data_mtxir2$month <- sub("^ACR\\.50\\.(.*)\\.(n|N)$", "\\1", ACR50_big_data_mtxir2$month)

##Solving error in Blanco 2017 placebo arm
#ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Blanco 2017",][1,c(4,5)]<-c(15,138)#responders not reported for placebo arm acr50

# ##Solving error in Burnmester 2013 placebo arm
# ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",][1,c(4,5)]<-c(11,131)
# ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",][2,c(4,5)]<-c(35,132)
# ACR50_big_data_biologic2[ACR50_big_data_biologic2$study_name=="Burmester 2013",]$month<-"3.months"

##Drop furst study due to high ROB
ACR50_big_data_mtxir2<- ACR50_big_data_mtxir2 %>%
filter(!study_name %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                         "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                         "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                         "Kameda 2011","Li 2016","Tam 2012"))
# ACR50_big_data_mtxir2<-na.omit(ACR50_big_data_mtxir2)

#removing studies with NA in arms
ACR50_big_data_mtxir2 <- ACR50_big_data_mtxir2 %>%
  group_by(study_name) %>%
  filter(!any(is.na(responders))) %>%
  ungroup()
ACR50_big_data_mtxir2<- as.data.frame(ACR50_big_data_mtxir2)
ACR50_big_data_mtxir2$study<-renumber(ACR50_big_data_mtxir2,"study")

##Keystone 2008 graphic
ACR50_big_data_mtxir2[ACR50_big_data_mtxir2$study_name=="Keystone 2008","responders"]<-paste0("8888",round(c(0.076*199,0.371*393),0))

##Tanaka 2016 graphic
ACR50_big_data_mtxir2[ACR50_big_data_mtxir2$study_name=="Tanaka 2016","responders"]<-paste0("8888",round(c(0.070173*49,0.459204*24,0.577237*24),0))

row.names(ACR50_big_data_mtxir2)<-NULL
View(ACR50_big_data_mtxir2)
##Master table
master.medtable<-master.medtable_dmard_ir
master.medtable[1,"medname"]<- "MTX"
master.medtable[2,"medname"]<-   "PBO (csDMARD stopped)"

library(dplyr)
for(i in 1:nrow(ACR50_big_data_mtxir2)){
  ppp<- master.medtable%>% filter(ACR50_big_data_mtxir2$treatment[i] ==master.medtable$medname )
  ACR50_big_data_mtxir2$treatment[i]<-ppp[1,1]
}


ACR50_big_data_mtxir2$responders<-as.numeric(gsub(8888,'',ACR50_big_data_mtxir2$responders))


#####Formatting the month column
ACR50_big_data_mtxir2$month<-sub(".months","",ACR50_big_data_mtxir2$month)
ACR50_big_data_mtxir2$month<-as.numeric(ifelse(ACR50_big_data_mtxir2$month=="End.of.trial",ACR50_big_data_mtxir2$TrialDuration,ACR50_big_data_mtxir2$month))
View(ACR50_big_data_mtxir2)

##Creating initial legend
mtxIR.medlegend <- new.medtable(ACR50_big_data_mtxir2,'treatment')
mtxIR.medlegend

##NMA plot reason for narrowing of network
nma.plot(ACR50_big_data_mtxir2[,1],ACR50_big_data_mtxir2[,3],data = ACR50_big_data_mtxir2, trtname = as.vector(mtxIR.medlegend$medname))

################################################################################################################
#trying the narrowing a network function
#############################################################################################
ACR50_big_data_mtxir2

#Narrow 1
mtxIR_narrow<-narrow.all(ACR50_big_data_mtxir2)
mtxIR_narrow

# # #Narrow 2
# dmardIR_narrow2<-narrow.all(dmardIR_narrow[[2]])
# dmardIR_narrow2


#sub network plot 1
mtxIR.d<-mtxIR_narrow[[1]]
length(unique(mtxIR.d$treatment))
mtxIR_subl<-new.medtable(mtxIR.d,'treatment')
nma.plot(mtxIR.d[,1],mtxIR.d[,3],data = mtxIR.d, trtname = as.vector(mtxIR_subl$medname), title = "ACR50 Network plot")



##############################################################################
#Performing sub network plot 1
library(rjags)
library(coda)
library(gemtc)
View(mtxIR.d)
mtxIR_subl

# network_data_ACR50_dmardIR_1<-dmardIR.d%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
#            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
network_data_ACR50_mtxIR_1<-mtxIR.d%>%
  select(c("study", "study_name","treatment","responders", "sampleSize", "month"))

network_data_ACR50_mtxIR_1$treatment<-renumber(network_data_ACR50_mtxIR_1,"treatment")
##To download acr50 data
network_data_ACR50_mtxIR_data_download<-network_data_ACR50_mtxIR_1
network_data_ACR50_mtxIR_data_download$treatment<-selected_treatment_name(network_data_ACR50_mtxIR_data_download$treatment,mtxIR_subl)
network_data_ACR50_mtxIR_data_download$study<-renumber(network_data_ACR50_mtxIR_data_download,"study")
write.csv(network_data_ACR50_mtxIR_data_download,"ACR50 MTX-IR dataset_OCT2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("ACR50 Network plot-MTX-IR with named interventions_OCT2026.jpg", width = 3480, height = 2050, res=300)
nma.plot(network_data_ACR50_mtxIR_1$study,network_data_ACR50_mtxIR_1$treatment,data = network_data_ACR50_mtxIR_1, trtname = as.vector(mtxIR_subl$medname),title = "ACR50 Network plot")
# Close the pdf file
dev.off() 
# 1. Open jpeg file
jpeg("ACR50 Network plot-MTX-IR without named interventions_OCT2026.jpg", width = 3480, height = 2050, res=300)
nma.plot(network_data_ACR50_mtxIR_1$study,network_data_ACR50_mtxIR_1$treatment,data = network_data_ACR50_mtxIR_1,title = "ACR50 Network plot")
# Close the pdf file
dev.off() 

network_data_ACR50_mtxIR_1<-mtc.network(network_data_ACR50_mtxIR_1)


##RUN the model
#mtc.hy.empirical.lor("semi-objective", "non-pharma")
#mtc.hy.prior("std.dev", "dunif", 0, 2)
#mtc.hy.prior("prec", "dgamma", 0.01, 0.01)
set.seed(505)
network_data_ACR50_mtxIR_1.model <- mtc.model(network_data_ACR50_mtxIR_1, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_ACR50_mtxIR_1.model$inits<- mapply(c, network_data_ACR50_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR<-mtc.run(network_data_ACR50_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR)

View(relative.effect.table.OR(results_mtxIR))
refOR_mtxIR<-as.data.frame(relative.effect.table.OR2(results_mtxIR, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refOR_mtxIR)<-mapvalues(row.names(refOR_mtxIR),
                                    from = c(1:(nrow(mtxIR_subl))),
                                    to = as.vector(mtxIR_subl$medname))
colnames(refOR_mtxIR)<-mapvalues(colnames(refOR_mtxIR),
                                   from = c(1:(nrow(mtxIR_subl))),
                                   to = as.vector(mtxIR_subl$medname))
View(refOR_mtxIR)
# which(row(refOR)==col(refOR))
# refOR[row(refOR)==col(refOR)]
index_diagonal2<-data.frame(index1=1:nrow(refOR_mtxIR),index2=1:nrow(refOR_mtxIR))
options(warn=-1)
for(i in 1:nrow(refOR_mtxIR)){
  refOR_mtxIR[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_mtxIR[row(refOR_mtxIR)==col(refOR_mtxIR)])[i],mtxIR_subl)
}
options(warn=0)
View(refOR_mtxIR)

write.csv(refOR_mtxIR,"ACR50 -MTX-IR pairwise comparisons_OCTOBER2026.csv", row.names = T)

library(openxlsx)
write.xlsx(refOR_mtxIR,"Results for MTX-IR_OCT2026.xlsx",sheetName = "ACR50",rowNames = T)


########################
########################

nma_estimate_mtxir_dataframe<-nma.estimate.table.OR(results_mtxIR)
View(nma_estimate_mtxir_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_ACR50_mtxIR_1.ns<-mtc.nodesplit(network_data_ACR50_mtxIR_1, thin=50,likelihood="binom", link="logit",hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 
names(network_data_ACR50_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_ACR50_mtxIR_1)

# Overall summary and plot
summary.network_data_ACR50_mtxIR_1.ns <- summary(network_data_ACR50_mtxIR_1.ns)

#Direct effects
direct.data<-summary.network_data_ACR50_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_ACR50_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_ACR50_mtxIR_1.ns$cons.effect

#P-values
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}
P_values_nodesplit_acr50mtxIR<-summary.network_data_ACR50_mtxIR_1.ns$p.value
P_values_nodesplit_acr50mtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_acr50mtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_acr50mtxIR2<-P_values_nodesplit_acr50mtxIR
P_values_nodesplit_acr50mtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_acr50mtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_acr50mtxIR2$comb<-paste0(P_values_nodesplit_acr50mtxIR2$t1,":",P_values_nodesplit_acr50mtxIR2$t2)


direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(exp(direct.data[,c(3:5)]),2),round(exp(indirect.data[,c(3:5)]),2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","OR1","low1","high1","OR2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                round(1/direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","OR1","low1","high1","OR2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)
grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxir_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)


#Getting the quality assessment for each study
ROB_quality_assess_data_mtxIR<-read.csv("review-60381-consensus-only-quality-assessment-study-per-row-20260627113107.csv")
quality_assess_data_consensus_mtxIR<-ROB_quality_assess_data_mtxIR[ROB_quality_assess_data_mtxIR$Reviewer=="Consensus",]
#ROB
ROB_quality_assess_data_consensus_mtxIR<-quality_assess_data_consensus_mtxIR[,c(1,seq(3,42,by=3))]
row.names(ROB_quality_assess_data_consensus_mtxIR)<- NULL

ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id=="vanRiel 2006", c("OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
                                                                                                      "OVERALL..ALL.OTHER.OUTCOMES")] <- "unclear"

###Adding ROB data
ffff_mtxIR<-unique(mtxIR.d$study_name)
data_qq_mtxIR<-mtxIR.d
View(ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff_mtxIR,])
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff_mtxIR,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<- current_assessment_for_studies_mtxIR %>%
  select(c(col_names))
list_appraisal_mtxIR<-list()
for(i in 1:length(ffff_mtxIR)){
  temp11<-data_qq_mtxIR[data_qq_mtxIR$study_name==ffff_mtxIR[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==ffff_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal_mtxIR[[i]]<-temp13
}
ACR50_mtxIR_withROB<-do.call(rbind,list_appraisal_mtxIR)
View(ACR50_mtxIR_withROB)

ACR50_mtxIR_withROB2<-ACR50_mtxIR_withROB %>%
  select(c("study","study_name", "treatment", "responders", "sampleSize",
           "month", "DoseLevel01", "TrialDuration","Time_point","OVERALL..ALL.OTHER.OUTCOMES"))
View(ACR50_mtxIR_withROB2)

###Adding Indirectness data
#ACR50_bio_withROB2_indrectness<-cbind(ACR50_bio_withROB2,Indirectness="low")
study_unique<-unique(ACR50_mtxIR_withROB2$study_name)
data_nn<-ACR50_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp[!is.na(resolved_indirectness_consensus_temp$`Final Rating`),]

View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]
row.names(current_indirect_assessment_for_studies)<-NULL
##because of bari 4mg
# current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:17),],current_indirect_assessment_for_studies[c(17),],current_indirect_assessment_for_studies[c(18:28),])
# current_indirect_assessment_for_studies$intervention[18]<-"Baricitinib (4 mg)"
# current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:13),],current_indirect_assessment_for_studies[c(13),],current_indirect_assessment_for_studies[c(14:25),])
# current_indirect_assessment_for_studies$intervention[14]<-"Baricitinib (4 mg)"
row.names(current_indirect_assessment_for_studies)<- NULL

list_appraisal2<-list()
for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
ACR50_try_data_new_withROB2_indrectness<-do.call(rbind,list_appraisal2)
ACR50_try_data_new_withROB2_indrectness<-  ACR50_try_data_new_withROB2_indrectness%>%
  dplyr::rename("Indirectness" = `Final Rating`)

ACR50_try_data_new_withROB2_indrectness$Indirectness<-tolower(ACR50_try_data_new_withROB2_indrectness$Indirectness)
View(ACR50_try_data_new_withROB2_indrectness)

ACR50_mtxIR_withROB2_indrectness<-ACR50_try_data_new_withROB2_indrectness
#####Formatting the month column
ACR50_mtxIR_withROB2_indrectness$month<-sub(".months","",ACR50_mtxIR_withROB2_indrectness$month)
ACR50_mtxIR_withROB2_indrectness$month<-as.numeric(ifelse(ACR50_mtxIR_withROB2_indrectness$month=="End.of.trial",ACR50_mtxIR_withROB2_indrectness$TrialDuration,ACR50_mtxIR_withROB2_indrectness$month))
View(ACR50_mtxIR_withROB2_indrectness)

ACR50_mtxIR_data_new_withROB2_indrectness2<-ACR50_mtxIR_withROB2_indrectness%>%
  select(c("study", "study_name","treatment","responders", "sampleSize", "month",
           "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))

ACR50_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(ACR50_mtxIR_data_new_withROB2_indrectness2,"treatment")
ACR50_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES<-ifelse(ACR50_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES=="unclear","some concerns",ACR50_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES)

##Check this,the ROB for vanRiel 2006 it was empty. Fixing temporary rating
ACR50_mtxIR_data_new_withROB2_indrectness2[ACR50_mtxIR_data_new_withROB2_indrectness2$study_name=="vanRiel 2006", "OVERALL..ALL.OTHER.OUTCOMES"]<- "some concerns"


##Loading the contribution matrix package, netmeta
library(netmeta)
pairwise_ACR50_OR<-pairwise(treat=treatment,event = responders,data=network_data_ACR50_mtxIR_1$data.ab,
                            studlab = study,n=sampleSize,sm="OR")

netPairwiseAcr50<-netmeta(pairwise_ACR50_OR)
netPairwiseAcr50_split<-netsplit(netPairwiseAcr50)
View(netPairwiseAcr50_split$direct.common)
View(netPairwiseAcr50_split$indirect.common)
netPairwiseAcr50_contribmat<-netcontrib(netPairwiseAcr50, method = "r")
View(netPairwiseAcr50_contribmat$random)
direct_evidence_new<-cbind(netPairwiseAcr50_split$direct.common[,c(1:5)],netPairwiseAcr50_split$indirect.common[,c(2:5)])

#getting only direct estimates
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]
library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))

direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")

#####Renaming the number treatment to named treatment
selected_treatment_name<-function(dataset,sub.medtable){
  suppressWarnings(library(dplyr))
  for(i in 1:length(dataset)){
    selected_treatment_name<- sub.medtable %>% filter(dataset[i] ==sub.medtable$new_mednum )
    dataset[i]<-selected_treatment_name[1,2]
  }
  return(dataset)
}

#Setting up initial grade table
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]
##Getting direct evidence for intervention in only multi-arm trials that fails to get direct estimate
all_actual_direct_evidence<-t(apply(apply(str_split_fixed(colnames(netPairwiseAcr50_contribmat$random), ":", 2),2,as.numeric),1,switch_func))
colnames(all_actual_direct_evidence)<-c("treat1", "treat2")
inthreeArmtrial<-anti_join(as.data.frame(all_actual_direct_evidence), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]



grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))

#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
acr50_output_df2_list<-list()
I_square_list<-list()


for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-ACR50_mtxIR_data_new_withROB2_indrectness2[ACR50_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:5)
  
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select("study", "responders.1", "sampleSize.1","responders.2", "sampleSize.2" )
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"OVERALL..ALL.OTHER.OUTCOMES"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","n.t","x.c","n.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    model.file<-"Dichot_direct_onestudy.txt"
    n.chains=3
    set.seed(1453)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(#p.c=c(0.5),p.t=c(0.5),
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "n.t"=temp12$n.t,"x.c"=temp12$x.c,"n.c"=temp12$n.c),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("OR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(OR=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
  }
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    model.file<-"Dichot_direct_multiplestudy_RE.txt"
    n.chains=3
    #model.inits <- list(p.c=c(rep(0.5,length(temp12$NS))),p.t=c(rep(0.5,length(temp12$NS))))
    set.seed(8663)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "n.t"=temp12$n.t,"x.c"=temp12$x.c,"n.c"=temp12$n.c, "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("OR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select("study", "responders.1", "sampleSize.1","responders.2", "sampleSize.2" )
  #nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(Studies_contributing_to_direct_evidence_wide$study),c(1,4,5,13,14)]
  nma_data50_4<-nma_data50_3
  colnames(nma_data50_4)<-c("NS","x.t","n.t","x.c","n.c")
  # nma_data50_4 <- nma_data50_4 %>%
  #               mutate(across(c("x.t","n.t","x.c","n.c"), as.numeric))
  
  row.names(nma_data50_4)<-NULL
  
  m.bin<-metabin(event.e = x.t, 
                 n.e = n.t,
                 event.c = x.c,
                 n.c = n.c,
                 # studlab = study_name,
                 data = nma_data50_4,
                 sm = "OR",
                 method = "MH",
                 MH.exact = TRUE,
                 common = FALSE,
                 random = TRUE,
                 method.tau = "PM",
                 method.random.ci = TRUE,
                 title = "ACR50"
                 
  )
  
  summary_m.bin<-summary(m.bin)
  contribution_weight<-as.character(round((summary_m.bin$w.random/sum(summary_m.bin$w.random))*100,2))
  pooled_estimate2<-data.frame(OR=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.bin$tau2,4)
                          ,
                          Lower=round(summary_m.bin$lower.tau2,4),
                          Upper=round(summary_m.bin$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.bin$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.bin$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.bin$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.bin$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  
  
  #}
  
  OR_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1)
  base_data <- tibble(mean  = OR_individual_studies$OR, 
                      lower = OR_individual_studies$Lower,
                      upper = OR_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      # treatment1 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      # treatment2 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      treatment1 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      treatment2 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      OR = paste(format(round(OR_individual_studies$OR,2),nsmall=2),paste("(",format(round(OR_individual_studies$Lower,2),nsmall=2),", ",format(round(OR_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$OR, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   OR = paste(format(round(pooled_estimate2$OR,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("ACR50", "Study"),
                   # treatment1 = c("",selected_treatment_name(unique(temp11$treatment.1),subl)),
                   # treatment2 = c("",selected_treatment_name(unique(temp11$treatment.2),subl)),
                   treatment1 = c("",selected_treatment_name(unique(temp11$treatment.2),mtxIR_subl)),
                   treatment2 = c("",selected_treatment_name(unique(temp11$treatment.1),mtxIR_subl)),
                   #ROB=c(paste("Favours",selected_treatment_name(unique(temp11$treatment.1),subl)),"ROB"),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   # OR = c(paste("Favours",selected_treatment_name(unique(temp11$treatment.2),subl)), "OR"),
                   OR = c("", "OR"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  acr50_output_df <- bind_rows(header,
                               base_data,
                               empty_row,
                               pooled)
  
  
  
  
  acr50_output_df2<-rbind(acr50_output_df,acr50_output_df[c(nrow(acr50_output_df)),])
  acr50_output_df2[nrow(acr50_output_df),1]<-I2
  acr50_output_df2[nrow(acr50_output_df),2]<-NA
  acr50_output_df2[nrow(acr50_output_df),c(6:10)]<-c(NA,NA,NA,NA,NA)
  # acr50_output_df2[nrow(acr50_output_df),2]<-c("Prediction Interval")
  # acr50_output_df2[nrow(acr50_output_df),6]<-paste(format(round(prediction_interval$Predictive,2),nsmall=2),paste("(",format(round(prediction_interval$Lower,2),nsmall=2),", ",format(round(prediction_interval$Upper,2),nsmall=2),")",sep=""))
  # acr50_output_df2[nrow(acr50_output_df),8]<-c(round(prediction_interval$Predictive,2))
  # acr50_output_df2[nrow(acr50_output_df),9]<-c(round(prediction_interval$Lower,2))
  # acr50_output_df2[nrow(acr50_output_df),10]<-c(round(prediction_interval$Upper,2))
  
  acr50_output_df2_list[[i]]<-acr50_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-acr50_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, OR,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(0.1, 2.4), 
               xlog = TRUE, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               # grid = structure(c(1.2), 
               #                  gp = gpar(lty = 2, col = "#F73831")),
               #txt_gp = fpTxtGp(cex=0.85),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               # xlab = paste("Outcome: ACR50-Direct Comparison",selected_treatment_name(unique(temp11$treatment.1),subl),"and",selected_treatment_name(unique(temp11$treatment.2),subl))
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_subl))
    )
  
}

length(graph_list)
graph_list[[21]]
acr50_output_df2_list[[7]]
n_studies<-c()
n_participants<-c()
for(i in 1:length(graph_list)){
  df<-acr50_output_df2_list[[i]][grep("\\d+/\\d+", acr50_output_df2_list[[i]]$treatment1),c(2,3) ]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Select denominator values
  df$col1_after_slash <- sapply(strsplit(df$treatment1, "/"), "[[", 2)
  df$col2_after_slash <- sapply(strsplit(df$treatment2, "/"), "[[", 2)
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(unlist(df[,c(3,4)]))))
  
}


Grade_assessment_dataset_direct<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_subl),
                                            selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_subl),
                                            Outcome="ACR50",
                                            `Effect measure`="OR",
                                            `Number of studies`=n_studies,
                                            `Number of participant`=n_participants,
                                            ROB="",Indirectness="",Imprecision="",
                                            Inconsistency="",
                                            Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                   "Intevention",
                                                                   "Effect measure",
                                                                   "Number of studies",
                                                                   "Number of participant",
                                                                   "Risk of Bias",
                                                                   "Publication Bias",
                                                                   "Certainty of evidence (ignores imprecision)")
View(Grade_assessment_dataset_direct)
#Saving grade assessment details for ACR50 BIOLOGICS as an object of lists
grade_details_acr50_mtxIR<-list(graph_list,acr50_output_df2_list,Grade_assessment_dataset_direct,I_square_list)

saveRDS(grade_details_acr50_mtxIR,file="List of ACR50 MTX-IR data components for GRADE_AUG2026.rds")

acr50_output_one_study_list <- list()
acr50_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

# Count studies for each comparison
n_studies <- sapply(acr50_output_df2_list, function(x) {
  
  df <- x[
    grep("\\d+/\\d+", x$treatment1),
    c(2,3)
  ]
  
  nrow(df)
})

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
acr50_output_one_study_list  <- acr50_output_df2_list[one_idx]
acr50_output_more_than_2_list <- acr50_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct[
    Grade_assessment_dataset_direct$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct[
    Grade_assessment_dataset_direct$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_acr50_mtxIR_one_study <- list(
  graph_list_one_study_list,
  acr50_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_acr50_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  acr50_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_acr50_mtxIR_one_study,
  file = "List of ACR50 MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_acr50_mtxIR_more_than_2,
  file = "List of ACR50 MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)


###############################################################################
###############################################################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct$Treatment, 
                                                Grade_assessment_dataset_direct$Intevention), 
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
#Getting the actual direct estimate for missing three arm trial
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(acr50_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","OR"][2]
}

grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


#############
##Creating indirect assessment 
####Get direct evidence report 
# setwd("/Users/mujaab/Documents/OneDrive - University of Calgary/Manucripts for PhD papers/Creating various datasets for GRADE evaluation")
# kamso_direct_evidence_rating_acr50<-read.csv("BIO-IR_ACR50_directAssessment_April2023.csv")
kamso_direct_evidence_rating_acr50<-read.csv("Glen_2026-09-03 22_16_ACR50-MTXIR_complete.csv")
kamso_direct_evidence_rating_acr502<-kamso_direct_evidence_rating_acr50[,-c(13:15)]
colnames(kamso_direct_evidence_rating_acr502)[1:2]<-colnames(grade_table4)[1:2]

kamso_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],kamso_direct_evidence_rating_acr502,
                                     by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)


kamso_direct_evidence_rating4<-cbind(kamso_direct_evidence_rating3[,c(1,2,4:7)],
                                     kamso_direct_evidence_rating3[,c(3,8:13)])
#View(kamso_direct_evidence_rating4[order(kamso_direct_evidence_rating4$Outcome),])##Perhaps might use this in the long run
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Outcome),"Outcome"]<-"ACR50"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"OR"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"

kamso_direct_evidence_rating5<-cbind(kamso_direct_evidence_rating4,
                                     num1=as.numeric(changing_treatment_name2num(kamso_direct_evidence_rating4$treat1,mtxIR_subl)),
                                     num2=as.numeric(changing_treatment_name2num(kamso_direct_evidence_rating4$treat2,mtxIR_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}

##ratings with comb
kamso_direct_evidence_rating5[,14:15]<-t(apply(kamso_direct_evidence_rating5[,14:15],1,switch_func2))
kamso_direct_evidence_rating6<-cbind(kamso_direct_evidence_rating5,
                                     comb=paste(kamso_direct_evidence_rating5[,14],
                                                kamso_direct_evidence_rating5[,15],sep = ":"))
kamso_direct_evidence_rating6_onlyDirect<-kamso_direct_evidence_rating6[!is.na(kamso_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
kamso_direct_evidence_rating6$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
##Getting indirect and nma estimates
Temp_GRADE<-merge(kamso_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(kamso_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("3:13","3:5"),"Indirect"]<-"NE"
#write.csv(Temp_GRADE4,"Initial GRADE table for ACR50_APR24.csv",row.names = F)
write.csv(Temp_GRADE4,"Initial GRADE table for ACR50_MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for ACR50_MTX-IR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR ACR50 OUTCOME
initial_gradeTable_ACR50<-Temp_GRADE5

estimate_contribution_matrix<-netPairwiseAcr50_contribmat$random
all_nma_estimates<-row.names(estimate_contribution_matrix)
flip_to_smaller_to_larger <- function(vector) {
  # Extract values before and after ":" and split them into two vectors
  before_colon <- as.numeric(sub(":.*", "", vector))
  after_colon <- as.numeric(sub(".*:", "", vector))
  
  # Determine which values are smaller
  smaller_values <- before_colon > after_colon
  
  # Create a new vector with flipped values
  flipped_vector <- ifelse(smaller_values, paste(after_colon, ":", before_colon, sep = ""), vector)
  
  return(flipped_vector)
}

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_ACR502<-initial_gradeTable_ACR50[match(all_nma_estimates,initial_gradeTable_ACR50$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))


#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()
#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_ACR502),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_ACR502.1<-cbind(initial_gradeTable_ACR502[,c(1:14)],first_order_loop_contri,
                                   initial_gradeTable_ACR502[,-c(1:14)])

initial_gradeTable_ACR502.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_ACR502.1,
                                                              all_nma_estimates,
                                                              estimate_contribution_matrix,
                                                              kamso_direct_evidence_rating6_onlyDirect,
                                                              top_n = 2))
# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_ACR502.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_ACR502.1<-result_df
initial_gradeTable_ACR502.1$combI<-str_replace(initial_gradeTable_ACR502.1$comb, ":", " vs ")

initial_gradeTable_ACR502.1<-process_dataset(initial_gradeTable_ACR502.1, 2)
initial_gradeTable_ACR503<-initial_gradeTable_ACR502.1[match(initial_gradeTable_ACR50$comb,initial_gradeTable_ACR502.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for ACR50_MTX-IR outcome_AUG2026.csv",row.names = T)

###########Assessing Intransitivity
##Here, since all the indirectness assessments are low, intransitivity will be not serious for all.
##Else manually evaluate intransitivity of each direct evidence contributing to the indirect evidence.
check_vector <- function(vector) {
  intransitivity_rating <- ifelse("Serious" %in% vector, "Serious",
                                  ifelse("Some Concerns" %in% vector, "Some Concerns", "Not Serious"))
  degree <- ifelse((intransitivity_rating == "Not Serious")|(intransitivity_rating == "Some Concerns"), 0, -1)
  
  return(c(degree, intransitivity_rating))
}


sub_data_intransitivity <- initial_gradeTable_ACR503[initial_gradeTable_ACR503$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_ACR503[initial_gradeTable_ACR503$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_ACR503.1<-cbind(initial_gradeTable_ACR503[,c(1:23)],
                                   Intransitivity="",
                                   `Certainty of evidence (ignores imprecision)-Indirect`="",
                                   initial_gradeTable_ACR503[,-c(1:23)])
initial_gradeTable_ACR503.1[initial_gradeTable_ACR503.1$Indirect!="-",c("Intransitivity",
                                                                        "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                              "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_ACR504<-initial_gradeTable_ACR503.1


#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)
grade_temp<-cbind(initial_gradeTable_ACR504[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_ACR504[,-c(1:25)])
###Applying the P-values
P_values_nodesplit_acr50mtxIR2
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_acr50mtxIR3 <- P_values_nodesplit_acr50mtxIR2[match(drt$comb, P_values_nodesplit_acr50mtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_acr50mtxIR3$p,2)

subACR50SOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatACR50<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subACR50SOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(nrow(subcontriMatACR50))
IndirPer <- numeric(nrow(subcontriMatACR50))

# Loop to calculate dirPer and IndirPer
for (i in 1:nrow(subcontriMatACR50)) {
  dirPer[i] <- subcontriMatACR50[i, colnames(subcontriMatACR50) == row.names(subcontriMatACR50)[i]]
  IndirPer[i] <- 1- dirPer[i]
}

# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")

#This code is wrong, change for all outcomes
# grade_temp2[grade_temp2$comb%in%row.names(subcontriMatACR50),"DirectPerContri"]<-round(dirPer,2)
# grade_temp2[grade_temp2$comb%in%row.names(subcontriMatACR50),"IndirectPerContri"]<-round(IndirPer,2)

grade_temp2[match(row.names(subcontriMatACR50), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(row.names(subcontriMatACR50), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)

##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
#When both direct and indirect are present
fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg

####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho


initial_gradeTable_ACR505<-grade_temp2.2
initial_gradeTable_ACR505<-initial_gradeTable_ACR505[,-c(30,33)]
colnames(initial_gradeTable_ACR505)[c(29,30)]<-c("Initial NMA evidence rating",
                                                 "P-Value")

###After incoherence evaluation
initial_gradeTable_ACR505.1<-cbind(initial_gradeTable_ACR505[,c(1:31)],
                                   qoeAfterIncohenrence="",
                                   initial_gradeTable_ACR505[,-c(1:31)])
initial_gradeTable_ACR505.1[(initial_gradeTable_ACR505.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_ACR505.1[(initial_gradeTable_ACR505.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_ACR505.1[initial_gradeTable_ACR505.1$Incoherence!="",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_ACR505.1[(initial_gradeTable_ACR505.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence

#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_ACR505.1[(initial_gradeTable_ACR505.1$Direct!="-")&(initial_gradeTable_ACR505.1$Indirect=="-"),"NMA"]<-initial_gradeTable_ACR505.1[(initial_gradeTable_ACR505.1$Direct!="-")&(initial_gradeTable_ACR505.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_ACR506<-cbind(initial_gradeTable_ACR505.1,
                                 p.estimate=apply(as.data.frame(initial_gradeTable_ACR505.1$NMA),1,function(x){
                                   as.numeric(sub(" .*", "", x))
                                 }),
                                 t(apply(as.data.frame(initial_gradeTable_ACR505.1$NMA),1,function(x){
                                   as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                   
                                 })))

colnames(initial_gradeTable_ACR506)[c(36,37)]<-c("l.bound","u.bound")

##Getting baseline risk using the placebo intervention
ACR50_data<-network_data_ACR50_mtxIR_1$data.ab%>%
  #select(-c(month, OVERALL..ALL.OTHER.OUTCOMES, Indirectness))
  select(-c(month))
ACR50_data_wide2<-reshape(number_of_arms(ACR50_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(ACR50_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(ACR50_data_wide2[,c(3,6,9)],1,xcv)

ACR50_data_wide3<-cbind(ACR50_data_wide2,na=numberOfArms)
View(ACR50_data_wide3)
colnames(ACR50_data_wide3)<-c("study", "study_name", "t.1","x.1","n.1",
                              "t.2","x.2","n.2","t.3","x.3","n.3",
                             # "t.4","x.4","n.4",
                              "na" )
#######
MTCData_dichot <- as.data.frame(ACR50_data_wide3)
x <- data.matrix(MTCData_dichot[,c(4,7,10)])
n <- data.matrix(MTCData_dichot[,c(5,8,11)])
t <- data.matrix(MTCData_dichot[,c(3,6,9)])
na <- data.matrix(MTCData_dichot[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <-length(t[which(t==1)])
xPBO<-x[which(t==1)]
nPBO <- n[which(t==1)]

# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Wichmann-Hill", .RNG.seed=1)

#loading rate model txt file
model.file4<-"Dichot_multiarm_RE.txt"
n.chains=4
model.inits <- initv1 

dichot.model3<-jags.model(file=model.file4,
                          data= list("x"=x,"t"=t,"na" = na[,1], "NT" = NT, "NS" = NS, "n"=n,
                                     "NTpbo"= NTpbo,"xPBO"=xPBO,"nPBO"=nPBO ),
                          n.chains=n.chains,
                          inits = model.inits,
                          n.adapt=10000
)

update(object=dichot.model3,n.iter=20000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("OR","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( dichot.model3, params4, n.iter=50000 )
summary_samps<-summary(samps4)
View(summary_samps$statistics)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]
typeof(estimates_relative_to_baseline)

#Converting estimates relative to baseline to absolutes
#nma_estimate_bio_dataframe3<-as.data.frame(apply(t(do.call(rbind,nma_estimate_bio_dataframe2)),2,as.numeric))
estimates_relative_to_baseline<-initial_gradeTable_ACR506[1:(NT-1),c(1,2,35:37)]

absolutes_from_baseline<-cbind(baseline_risk,(estimates_relative_to_baseline[,3:5]*baseline_risk)/((1-baseline_risk)+(estimates_relative_to_baseline[,3:5]*baseline_risk)))
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(baseline_risk,estimates_relative_to_baseline_combine[,7]))

colnames(all_baseline)<-"base"
com.list<-list()
unique_treat2<-unique(initial_gradeTable_ACR506$treat2)
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_ACR506[initial_gradeTable_ACR506$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,(d.sub[,3:5]*base)/((1-base)+(d.sub[,3:5]*base)))
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9]*1000)
more_fewer_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"fewer"
  }else{
    text<-"more"
  }
  return(paste(round(abs(e1-e2),0),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,13]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,13])))

prefix.assumedrisk <- ''
prefix.corrrisk<-''
suffix.assumedrisk <- ' per 1000'
suffix.corrrisk <- ' per 1000'
result_com.dataframe_withtext<-cbind(com.dataframe_withtext[,1:2],
                                     paste(prefix.assumedrisk,round(com.dataframe_withtext[,6]*1000,0), suffix.assumedrisk, sep=''),
                                     paste(prefix.corrrisk, round(com.dataframe_withtext[,7]*1000,0), suffix.corrrisk, ' (', round(com.dataframe_withtext[,8]*1000,0), ' to ',round(com.dataframe_withtext[,9]*1000,0),')', sep=''),
                                     paste(prefix.corrrisk, com.dataframe_withtext[,13], ' (', com.dataframe_withtext[,14], ' to ',com.dataframe_withtext[,15],')', sep=''),
                                     # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                     com.dataframe_withtext[,10:15]
)

colnames(result_com.dataframe_withtext)[c(3:11)]<-c("Assumed risk (comparator)",	"Corresponding risk (intervention)",
                                                    "Risk difference","p.estimate","l.bound","u.bound","","","")
View(cbind(initial_gradeTable_ACR506[,-c(35:37)],
           result_com.dataframe_withtext[,-c(1,2)]))

initial_gradeTable_ACR507<-cbind(initial_gradeTable_ACR506[,-c(35:37)],
                                 result_com.dataframe_withtext[,-c(1,2)],pair="")
#colnames(initial_gradeTable_ACR507)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_ACR508<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_ACR507$treat1,mtxIR_subl)),
                                 `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_ACR507$treat2,mtxIR_subl)),
                                 initial_gradeTable_ACR507)
colnames(initial_gradeTable_ACR508)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_ACR508<-initial_gradeTable_ACR508[,-c(46)]
#write.csv(initial_gradeTable_ACR508,"ACR50 SoF table without imprecision assessment May2023.csv",row.names = F)
write.csv(initial_gradeTable_ACR508,"ACR50 MTX-IR SoF table without imprecision assessment AUG2026.csv",row.names = F)

effect.size.rating_list<-lapply(split(initial_gradeTable_ACR507, seq_len(nrow(initial_gradeTable_ACR507))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 50, large.effect = 200, framework = T)
})

as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_ACR509<-initial_gradeTable_ACR508
initial_gradeTable_ACR509<-cbind(initial_gradeTable_ACR509[,c(1:34,37:45,35:36)])
initial_gradeTable_ACR509.1<-cbind(initial_gradeTable_ACR509[,c(1:37)],
                                   Small.effect.threshold= 50,
                                   Large.effect.threshold= 200,
                                   Framework = "NULL",
                                   Imprecision.1 = "",
                                   Certainty.of.evidence = "",
                                   initial_gradeTable_ACR509[,c(44,45)])


initial_gradeTable_ACR509.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_ACR50_mtxIR<-initial_gradeTable_ACR509.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_ACR50_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_ACR50_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_ACR50_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_ACR50_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_ACR50_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_ACR50_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_ACR50_mtxIR$Certainty.of.evidence[i]<-rating_impre
}



View(detailed_SOF_ACR50_mtxIR[is.na(detailed_SOF_ACR50_mtxIR$Certainty.of.evidence),])
detailed_SOF_ACR50_mtxIR[is.na(detailed_SOF_ACR50_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

#detailed_SOF_ACR50_mtxIR$comb<-initial_gradeTable_ACR507$comb

write.csv(detailed_SOF_ACR50_mtxIR,"ACR50_MTX-IR_Detail_SoF_OCT2026.csv", row.names = F)



######################################################################################
######################################################################################
#############HAQ Analysis##################
############################
#######Graphical data for HAQ and errors
#######Kim 2000
##3 months
final_mtx_irDataset_actual[final_mtx_irDataset_actual$Study.Identifier=="Kim 2000","HAQ.DI.3.months...actual.value.Mean"]<-c(NA,NA)
final_mtx_irDataset_actual[final_mtx_irDataset_actual$Study.Identifier=="Kim 2000","HAQ.DI.3.months...change.from.baseline.Mean"]<-c(-0.30,0.10)
final_mtx_irDataset_actual[final_mtx_irDataset_actual$Study.Identifier=="Kim 2000","HAQ.DI.3.months...change.from.baseline.SD"] <-exactp.sd(0.001, 17, 17, 0.10, -0.30)[[2]]

##Creating HAQ-DI dataset
dataset_HAQ <- final_mtx_irDataset_actual%>% 
  select(matches("^HAQ\\.DI\\..*")) %>%
  select(-contains("Reported"))

colnames(dataset_HAQ)
dataset_HAQ[dataset_HAQ=="NA"]<-NA
View(dataset_HAQ)
actuals<-dataset_HAQ %>%
  select(matches("^HAQ.DI\\.12\\.months..actual*"))

View(actuals)
colnames(actuals)

library(tibble)# library to add column
actuals<-add_column(actuals, "HAQ.DI.12.months..actual.value.Mean_edit" = actuals$HAQ.DI.12.months..actual.value.Mean, .after = "HAQ.DI.12.months..actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
index<-which(is.na(actuals$HAQ.DI.12.months..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){#evaluate when entire dataset is NA
  actuals<-add_column(actuals, "HAQ.DI.12.months..actual.value.SD_edit" = actuals$HAQ.DI.12.months..actual.value.SD, .after = "HAQ.DI.12.months..actual.value.SD")
  changeFbaseline<-dataset_HAQ%>%
    select(matches("^HAQ.DI\\.12.months\\..change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.12.months..change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.12.months..change.from.baseline.Mean, .after = "HAQ.DI.12.months..change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.12.months..change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.12.months..change.from.baseline.SD, .after = "HAQ.DI.12.months..change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"HAQ.DI.12.months..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"HAQ.DI.12.months..actual.value.Mean_edit"]<-sub_actuals[index2,"HAQ.DI.12.months..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  actuals<-add_column(actuals, "HAQ.DI.12.months..actual.value.SD_edit" = actuals$HAQ.DI.12.months..actual.value.SD, .after = "HAQ.DI.12.months..actual.value.SD")
  
  #for missing variability
  index3<-which(is.na(actuals$HAQ.DI.12.months..actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"HAQ.DI.12.months..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"HAQ.DI.12.months..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"HAQ.DI.12.months..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.12.months..actual.value.Upper.CI))
      temp_sd[index7,"HAQ.DI.12.months..actual.value.SD_edit"]<-((temp_sd[index7,"HAQ.DI.12.months..actual.value.Upper.CI"]-temp_sd[index7,"HAQ.DI.12.months..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.12.months..actual.value.N"])
      sub_actuals_sd[index6,] <-temp_sd
      actuals[index3,]<-sub_actuals_sd
    } 
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"HAQ.DI.12.months..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.12.months..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.12.months..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.12.months..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.12.months..actual.value.N"][which.max(temp_sd2[index9,"HAQ.DI.12.months..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.12.months..actual.value.N"][which.min(temp_sd2[index9,"HAQ.DI.12.months..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.12.months..actual.value.p.value"]))
      values<-as.numeric(c(p.val, Nrx, Npbo, meanhigh, meanlow))
      temp_sd2[index9,"HAQ.DI.12.months..actual.value.SD_edit"]<- exactp.sd(p.val=values[1], Nrx=values[2], Npbo=values[3], meanhigh=values[4], meanlow=values[5])[[2]]
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####12 months Change from baseline
  changeFbaseline<-dataset_HAQ%>%
    select(matches("^HAQ.DI\\.12.months\\..change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.12.months..change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.12.months..change.from.baseline.Mean , .after = "HAQ.DI.12.months..change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$HAQ.DI.12.months..change.from.baseline.Mean))
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"HAQ.DI.12.months..change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"HAQ.DI.12.months..change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"HAQ.DI.12.months..change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.12.months..change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.12.months..change.from.baseline.SD, .after = "HAQ.DI.12.months..change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$HAQ.DI.12.months..change.from.baseline.SD))
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"HAQ.DI.12.months..change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"HAQ.DI.12.months..change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"HAQ.DI.12.months..change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.12.months..change.from.baseline.Upper.CI))
      temp_sd[index7,"HAQ.DI.12.months..change.from.baseline.SD_edit"]<-((as.numeric(temp_sd[index7,"HAQ.DI.12.months..change.from.baseline.Upper.CI"])-as.numeric(temp_sd[index7,"HAQ.DI.12.months..change.from.baseline.Lower.CI"]))/3.92)*sqrt(as.numeric(temp_sd[index7,"HAQ.DI.12.months..change.from.baseline.N"]))
      sub_changeFbaseline_sd[index6,] <-temp_sd
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.12.months..change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.12.months..change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.N"][which.max(temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.N"][which.min(temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.p.value"]))
      values<-as.numeric(c(p.val, Nrx, Npbo, meanhigh, meanlow))
      temp_sd2[index9,"HAQ.DI.12.months..change.from.baseline.SD_edit"]<- exactp.sd(p.val=values[1], Nrx=values[2], Npbo=values[3], meanhigh=values[4], meanlow=values[5])[[2]]
      # temp_sd2[index9,"DAS.12.months..change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.12.months..change.from.baseline.p.value"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.Mean"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}


##Replacing processed HAQ outcome-12months
HAQ_12months_continuous_outcome<-cbind(actuals,changeFbaseline)


####3 months Actuals
dataset_HAQ_3months<-dataset_HAQ%>%
  select(matches("^HAQ.DI\\.3\\.months*"))
dataset_HAQ_3months[dataset_HAQ_3months=="NA"]<-NA
View(dataset_HAQ_3months)
colnames(dataset_HAQ_3months)
actuals<-dataset_HAQ_3months%>%
  select(matches("^HAQ.DI\\.3\\.months...actual*"))
View(actuals)
colnames(actuals)
#for missing means
actuals<-add_column(actuals, "HAQ.DI.3.months...actual.value.Mean_edit" = actuals$HAQ.DI.3.months...actual.value.Mean, .after = "HAQ.DI.3.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$HAQ.DI.3.months...actual.value.Mean))

if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "HAQ.DI.3.months...actual.value.SD_edit" = actuals$HAQ.DI.3.months...actual.value.SD, .after = "HAQ.DI.3.months...actual.value.SD")
  changeFbaseline<-dataset_HAQ_3months%>%
    select(matches("^HAQ.DI\\.3.months\\...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.3.months...change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.3.months...change.from.baseline.Mean, .after = "HAQ.DI.3.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.3.months...change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.3.months...change.from.baseline.SD, .after = "HAQ.DI.3.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"HAQ.DI.3.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"HAQ.DI.3.months...actual.value.Mean_edit"]<-sub_actuals[index2,"HAQ.DI.3.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "HAQ.DI.3.months...actual.value.SD_edit" = actuals$HAQ.DI.3.months...actual.value.SD, .after = "HAQ.DI.3.months...actual.value.SD")
  index3<-which(is.na(actuals$HAQ.DI.3.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"HAQ.DI.3.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"HAQ.DI.3.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"HAQ.DI.3.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.3.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.3.months...actual.value.SD_edit"]<-((temp_sd[index7,"HAQ.DI.3.months...actual.value.Upper.CI"]-temp_sd[index7,"HAQ.DI.3.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.3.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"HAQ.DI.3.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.3.months...actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.3.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.3.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.3.months...actual.value.N"][which.max(temp_sd2[index9,"HAQ.DI.3.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.3.months...actual.value.N"][which.min(temp_sd2[index9,"HAQ.DI.3.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.3.months...actual.value.p.value"]))
      temp_sd2[index9,"HAQ.DI.3.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####3 months Change from baseline
  changeFbaseline<-dataset_HAQ_3months%>%
    select(matches("^HAQ.DI\\.3.months\\...change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  # View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.3.months...change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.3.months...change.from.baseline.Mean, .after = "HAQ.DI.3.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$HAQ.DI.3.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"HAQ.DI.3.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"HAQ.DI.3.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[as.character(index2),"HAQ.DI.3.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.3.months...change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.3.months...change.from.baseline.SD, .after = "HAQ.DI.3.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$HAQ.DI.3.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"HAQ.DI.3.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"HAQ.DI.3.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"HAQ.DI.3.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.3.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.3.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"HAQ.DI.3.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"HAQ.DI.3.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.3.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.3.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.3.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.N"][which.max(temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.N"][which.min(temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"HAQ.DI.3.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}




##Replacing processed HAQ outcome-3months
HAQ_3months_continuous_outcome<-cbind(actuals,changeFbaseline)


#################HAQ 6 months
####6 months Actuals
dataset_HAQ_6months<-dataset_HAQ%>%
  select(matches("^HAQ.DI\\.6\\.months*"))
dataset_HAQ_6months[dataset_HAQ_6months=="NA"]<-NA
View(dataset_HAQ_6months)
colnames(dataset_HAQ_6months)
actuals<-dataset_HAQ_6months%>%
  select(matches("^HAQ.DI\\.6\\.months...actual*"))
View(actuals)
colnames(actuals)

#for missing means
actuals<-add_column(actuals, "HAQ.DI.6.months...actual.value.Mean_edit" = actuals$HAQ.DI.6.months...actual.value.Mean, .after = "HAQ.DI.6.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$HAQ.DI.6.months...actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "HAQ.DI.6.months...actual.value.SD_edit" = actuals$HAQ.DI.6.months...actual.value.SD, .after = "HAQ.DI.6.months...actual.value.SD")
  changeFbaseline<-dataset_HAQ_6months%>%
    select(matches("^HAQ.DI\\.6.months\\...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.6.months...change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.6.months...change.from.baseline.Mean, .after = "HAQ.DI.6.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.6.months...change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.6.months...change.from.baseline.SD, .after = "HAQ.DI.6.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"HAQ.DI.6.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"HAQ.DI.6.months...actual.value.Mean_edit"]<-sub_actuals[index2,"HAQ.DI.6.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "HAQ.DI.6.months...actual.value.SD_edit" = actuals$HAQ.DI.6.months...actual.value.SD, .after = "HAQ.DI.6.months...actual.value.SD")
  index3<-which(is.na(actuals$HAQ.DI.6.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"HAQ.DI.6.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"HAQ.DI.6.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"HAQ.DI.6.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.6.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.6.months...actual.value.SD_edit"]<-((temp_sd[index7,"HAQ.DI.6.months...actual.value.Upper.CI"]-temp_sd[index7,"HAQ.DI.6.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.6.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"HAQ.DI.6.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.6.months...actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.6.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.6.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.6.months...actual.value.N"][which.max(temp_sd2[index9,"HAQ.DI.6.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.6.months...actual.value.N"][which.min(temp_sd2[index9,"HAQ.DI.6.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.6.months...actual.value.p.value"]))
      temp_sd2[index9,"HAQ.DI.6.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####6 months Change from baseline
  changeFbaseline<-dataset_HAQ_6months%>%
    select(matches("^HAQ.DI\\.6.months\\...change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.6.months...change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.6.months...change.from.baseline.Mean, .after = "HAQ.DI.6.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$HAQ.DI.6.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"HAQ.DI.6.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"HAQ.DI.6.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"HAQ.DI.6.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.6.months...change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.6.months...change.from.baseline.SD, .after = "HAQ.DI.6.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$HAQ.DI.6.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"HAQ.DI.6.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"HAQ.DI.6.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"HAQ.DI.6.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.6.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.6.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"HAQ.DI.6.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"HAQ.DI.6.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.6.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.6.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.6.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.N"][which.max(temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.N"][which.min(temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"HAQ.DI.6.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
  
  
}

##Replacing processed DAS28 outcome-6months
HAQ_6months_continuous_outcome<-cbind(actuals,changeFbaseline)


options(warn = defaultW)


################HAQ baseline
####HAQ baseline Actuals
dataset_HAQBaseline<-dataset_HAQ%>%
  select(matches("^HAQ.DI\\.Baseline*"))
dataset_HAQBaseline[dataset_HAQBaseline=="NA"]<-NA
View(dataset_HAQBaseline)
actuals<-dataset_HAQBaseline
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "HAQ.DI.Baseline.Mean_edit" = actuals$HAQ.DI.Baseline.Mean, .after = "HAQ.DI.Baseline.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$HAQ.DI.Baseline.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "HAQ.DI.Baseline.SD_edit" = actuals$HAQ.DI.Baseline.SD, .after = "HAQ.DI.Baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"HAQ.DI.Baseline.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"HAQ.DI.Baseline.Mean_edit"]<-sub_actuals[index2,"HAQ.DI.Baseline.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "HAQ.DI.Baseline.SD_edit" = actuals$HAQ.DI.Baseline.SD, .after = "HAQ.DI.Baseline.SD")
  index3<-which(is.na(actuals$HAQ.DI.Baseline.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"HAQ.DI.Baseline.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"HAQ.DI.Baseline.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"HAQ.DI.Baseline.SD_edit"]<-(sub_actuals_sd[index4,"HAQ.DI.Baseline.SE"])*sqrt(sub_actuals_sd[index4,"HAQ.DI.Baseline.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.Baseline.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"HAQ.DI.Baseline.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.Baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.Baseline.SD_edit"]<-((temp_sd[index7,"HAQ.DI.Baseline.Upper.CI"]-temp_sd[index7,"HAQ.DI.Baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.Baseline.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.Baseline.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"HAQ.DI.Baseline.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.Baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.Baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.Baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.Baseline.N"][which.max(temp_sd2[index9,"HAQ.DI.Baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.Baseline.N"][which.min(temp_sd2[index9,"HAQ.DI.Baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.Baseline.p.value"]))
      temp_sd2[index9,"HAQ.DI.Baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
}


options(warn = defaultW)
##Replacing processed DAS28 outcome-Baseline
HAQBaseline_continuous_outcome<-actuals


#################HAQ End of trial
####End of trial Actuals
dataset_HAQ_Endoftrial<-dataset_HAQ%>%
  select(matches("^HAQ.DI\\.End.of.trial*"))
dataset_HAQ_Endoftrial[dataset_HAQ_Endoftrial=="NA"]<-NA
View(dataset_HAQ_Endoftrial)
actuals<-dataset_HAQ_Endoftrial%>%
  select(matches("^HAQ.DI\\.End.of.trial\\..actual*"))
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "HAQ.DI.End.of.trial..actual.value.Mean_edit" = actuals$HAQ.DI.End.of.trial..actual.value.Mean, .after = "HAQ.DI.End.of.trial..actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
index<-which(is.na(actuals$HAQ.DI.End.of.trial..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "HAQ.DI.End.of.trial..actual.value.SD_edit" = actuals$HAQ.DI.End.of.trial..actual.value.SD, .after = "HAQ.DI.End.of.trial..actual.value.SD")
  changeFbaseline<-dataset_HAQ_Endoftrial%>%
    select(matches("^HAQ.DI\\.End.of.trial\\.change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.Mean, .after = "HAQ.DI.End.of.trial.change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.SD, .after = "HAQ.DI.End.of.trial.change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"HAQ.DI.End.of.trial..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"HAQ.DI.End.of.trial..actual.value.Mean_edit"]<-sub_actuals[index2,"HAQ.DI.End.of.trial..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "HAQ.DI.End.of.trial..actual.value.SD_edit" = actuals$HAQ.DI.End.of.trial..actual.value.SD, .after = "HAQ.DI.End.of.trial..actual.value.SD")
  index3<-which(is.na(actuals$HAQ.DI.End.of.trial..actual.value.SD))
  
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"HAQ.DI.End.of.trial..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"HAQ.DI.End.of.trial..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"HAQ.DI.End.of.trial..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.End.of.trial..actual.value.Upper.CI))
      
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.End.of.trial..actual.value.SD_edit"]<-((temp_sd[index7,"HAQ.DI.End.of.trial..actual.value.Upper.CI"]-temp_sd[index7,"HAQ.DI.End.of.trial..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.End.of.trial..actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"HAQ.DI.End.of.trial..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.End.of.trial..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.N"][which.max(temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.N"][which.min(temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.p.value"]))
      temp_sd2[index9,"HAQ.DI.End.of.trial..actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  options(warn = defaultW)
  
  
  ####End of trial Change from baseline
  changeFbaseline<-dataset_HAQ_Endoftrial%>%
    select(matches("^HAQ.DI\\.End.of.trial\\.change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #colnames(changeFbaseline)
  defaultW <- getOption("warn") 
  options(warn = -1) 
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.Mean, .after = "HAQ.DI.End.of.trial.change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"HAQ.DI.End.of.trial.change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"HAQ.DI.End.of.trial.change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"HAQ.DI.End.of.trial.change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "HAQ.DI.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.SD, .after = "HAQ.DI.End.of.trial.change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$HAQ.DI.End.of.trial.change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"HAQ.DI.End.of.trial.change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"HAQ.DI.End.of.trial.change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"HAQ.DI.End.of.trial.change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$HAQ.DI.End.of.trial.change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"HAQ.DI.End.of.trial.change.from.baseline.SD_edit"]<-((temp_sd[index7,"HAQ.DI.End.of.trial.change.from.baseline.Upper.CI"]-temp_sd[index7,"HAQ.DI.End.of.trial.change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"HAQ.DI.End.of.trial.change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"HAQ.DI.End.of.trial.change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$HAQ.DI.End.of.trial.change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.N"][which.max(temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.N"][which.min(temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.p.value"]))
      temp_sd2[index9,"HAQ.DI.End.of.trial.change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}

options(warn = defaultW)
##Replacing processed DAS28 outcome-EndOfTrial
HAQ_End_of_trial_continuous_outcome<-cbind(actuals,changeFbaseline)

####Final Continuous outcome dataset-HAQ
final_mtxIRDataset_HAQ<-cbind(final_mtx_irDataset_actual[,c(1:27,30:37,40:43,47:57)],
                                 HAQ_12months_continuous_outcome,
                                 HAQ_3months_continuous_outcome,
                                 HAQ_6months_continuous_outcome,
                                 HAQBaseline_continuous_outcome,
                                 HAQ_End_of_trial_continuous_outcome,
                                 DoseLevel01=final_mtx_irDataset_actual[,c("DoseLevel01")])
final_mtxIRDataset_HAQ[,c(51:141)] <- lapply(final_mtxIRDataset_HAQ[,c(51:141)], function(x) as.numeric(as.character(x)))#converting character numeric to numeric

View(final_mtxIRDataset_HAQ)

write.csv(final_mtxIRDataset_HAQ,"HAQ MTX-IR_AUG2026.csv",row.names = F)


##fill time points of no interest with NA
final_mtxIRDataset_HAQ2 <- final_mtxIRDataset_HAQ

# End of trial columns
eot_cols <- grep(
  "^HAQ\\.DI\\.End\\.of\\.trial",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)

# 12-month columns
m12_cols <- grep(
  "^HAQ\\.DI\\.12\\.months",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)

# 6-month columns
m6_cols <- grep(
  "^HAQ\\.DI\\.6\\.months",
  names(final_biologicDataset_HAQ2),
  value = TRUE
)

# 3-month columns
m3_cols <- grep(
  "^HAQ\\.DI\\.3\\.months",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)

# Apply same rules as dichotomous outcome
final_mtxIRDataset_HAQ2[
  final_mtxIRDataset_HAQ2$`Trial.duration..in.months` %in%
    c(3,4,6,9,11,12,13,18),
  eot_cols
] <- NA

final_mtxIRDataset_HAQ2[
  final_mtxIRDataset_HAQ2$`Trial.duration..in.months` %in%
    c(3,4,6,23,24,63),
  m12_cols
] <- NA

final_mtxIRDataset_HAQ2[
  final_mtxIRDataset_HAQ2$`Trial.duration..in.months` %in%
    c(3,4,9,11,12,13,18,23,24,63),
  m6_cols
] <- NA

final_mtxIRDataset_HAQ2[
  final_mtxIRDataset_HAQ2$`Trial.duration..in.months` %in%
    c(6,9,11,12,13,18,23,24,63),
  m3_cols
] <- NA

##Emery 2006
months6_cols<-grep(
  "^HAQ\\.DI\\.6\\.months",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Emery 2006",months6_cols]<-NA

#Tanaka 2011
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Tanaka 2011",c("HAQ.DI.Baseline.SD","HAQ.DI.Baseline.SD_edit")]<-c(0.37,0.27)

##Emery 2006, no baseline
months6_cols<-grep(
  "^HAQ\\.DI\\.6\\.months",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Emery 2010",months6_cols]<-NA

##Bi 2019, no baseline
months6_cols<-grep(
  "^HAQ\\.DI\\.6\\.months",
  names(final_mtxIRDataset_HAQ2),
  value = TRUE
)
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Bi 2019",months6_cols]<-NA

#Fisher 2016
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Fisher 2016","HAQ.DI.3.months...change.from.baseline.Mean_edit"]<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Fisher 2016","HAQ.DI.3.months...change.from.baseline.Median"]

#Kim 2000
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Kim 2000",c("HAQ.DI.Baseline.SD","HAQ.DI.Baseline.SD_edit")]<-c(0.30,0.30)

#Lan 2004
final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Lan 2004",c("HAQ.DI.3.months...actual.value.SD","HAQ.DI.3.months...actual.value.SD_edit")]<-exactp.sd(0.021, 29, 29, 0.99, 0.34)[[2]]

#Keystone 2009

final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Keystone 2009",c("HAQ.DI.Baseline.SD_edit")]<-c(c(iqr2sd(0.750, 1.750)),c(iqr2sd(1.000, 1.875)))


final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Keystone 2009",c("HAQ.DI.6.months...change.from.baseline.Mean_edit")]<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study.Identifier=="Keystone 2009",c("HAQ.DI.6.months...change.from.baseline.Median")]

temp_data3<-list()
empty_entries3<-list()

tttt<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study_no==16,c(51:140)]
tttt_studyInfo<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study_no==16,c(1:3)]
tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
tttt
as.vector(which(colSums(!is.na(tttt)) > 0))
colnames(tttt)

getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
getting_actuals2<-getting_actuals[, grepl("actual.value",names(getting_actuals))]
View(getting_actuals2)
dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
# get the last 3 columns
last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]



for(i in 1:length(unique(final_mtxIRDataset_HAQ2$Study_no))){
  tttt<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study_no==i,c(51:140)]
  tttt_studyInfo<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study_no==i,c(1:3)]
  tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
  
  if(sum(which(colSums(!is.na(tttt)) > 0))==0){
    empty_entries3[[i]]<-final_mtxIRDataset_HAQ2[final_mtxIRDataset_HAQ2$Study_no==i,c(1:3,52:140)]
  }else if(identical(c(85:89),as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##Only baseline value
    eee<-cbind(tttt_studyInfo,tttt[colSums(!is.na(tttt)) > 0]%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(tttt[colSums(!is.na(tttt)) > 0])[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data3[[i]]<- eee
  }else if(any(c(5:9,25:29,45:49,75:79)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##actual only but need to choose month
    getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
    getting_actuals2<-getting_actuals[, grepl("actual.value",names(getting_actuals))]
    if(ncol(getting_actuals2)==5){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[3])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==6){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[4])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==10){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==11){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==15){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[13])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==16){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[15])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==18){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[17])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==20){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[18])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }
  }else if(any(c(11,31,51,71)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##Baseline value and change in baseline value
    n_col<-ncol(tttt[colSums(!is.na(tttt)) > 0])
    #baseline<-tttt[colSums(!is.na(tttt)) > 0][,sort(n_col:((n_col-5)+1))]
    baseline<-tttt[colSums(!is.na(tttt)) > 0][,grep("^HAQ\\.DI\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE)]
    #change_mean_values<-tttt[colSums(!is.na(tttt)) > 0][,-sort(n_col:((n_col-5)+1))]
    change_mean_values<-tttt[colSums(!is.na(tttt)) > 0]%>%
      select(!grep("^HAQ\\.DI\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE))
    change_mean_values2<-as.data.frame(change_mean_values[ , grepl( "change.from.baseline.Mean_edit" , names(change_mean_values) ) ])
    actual_mean<-baseline[,2]+change_mean_values2[,ncol(change_mean_values2)]
    baseline[,"HAQ.DI.Baseline.Mean_edit"]<-actual_mean
    replaced_baseline_mean_with_actual<-baseline
    eee<-cbind(tttt_studyInfo,replaced_baseline_mean_with_actual%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(replaced_baseline_mean_with_actual)[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data3[[i]]<- eee
  }
  
  
}


HAQ_dataset_analysis<-do.call(rbind,temp_data3)
HAQ_dataset_analysis_mtxIR<-HAQ_dataset_analysis


HAQ_dataset_analysis_mtxIR$SD<-round(as.numeric(HAQ_dataset_analysis_mtxIR$SD),2)
HAQ_dataset_analysis_mtxIR$Mean<-round(as.numeric(HAQ_dataset_analysis_mtxIR$Mean),2)

##Drop furst study due to high ROB
HAQ_dataset_analysis_mtxIR<- HAQ_dataset_analysis_mtxIR %>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                           "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                           "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                           "Kameda 2011","Li 2016","Tam 2012"))

HAQ_dataset_analysis_mtxIR$Study_no<-renumber(HAQ_dataset_analysis_mtxIR,"Study_no")


#HAQ_dataset_analysis1<-HAQ_dataset_analysis1[-16,]
#####Formatting the month column
HAQ_dataset_analysis_mtxIR$Timepoint <- gsub(
  "HAQ\\.DI\\.|\\.months\\.{2,3}actual\\.value\\.(N|Mean_edit|Median|Lower\\.CI)|\\.(Mean_edit|Median)",
  "",
  HAQ_dataset_analysis_mtxIR$Timepoint
)


HAQ_dataset_analysis_mtxIR<- HAQ_dataset_analysis_mtxIR%>%
  dplyr::rename(treatment = Intervention)


library(dplyr)
for(i in 1:nrow(HAQ_dataset_analysis_mtxIR)){
  ppp<- master.medtable %>% filter(HAQ_dataset_analysis_mtxIR$treatment[i] ==master.medtable$medname )
  HAQ_dataset_analysis_mtxIR$treatment[i]<-ppp[1,1]
}


##Creating initial legend
mtxIR.HAQ.medlegend <- new.medtable(HAQ_dataset_analysis_mtxIR,'treatment')
mtxIR.HAQ.medlegend

##NMA plot reason for narrowing of network
nma.plot(HAQ_dataset_analysis_mtxIR[,1],HAQ_dataset_analysis_mtxIR[,3],data = HAQ_dataset_analysis_mtxIR, trtname = as.vector(mtxIR.HAQ.medlegend$medname))


################################################################################################################
#trying the narrowing a network function
#############################################################################################
HAQ_dataset_analysis_mtxIR
colnames(HAQ_dataset_analysis_mtxIR)[c(1,2)]<-c("study","study_name")

#Narrow 1
mtxIR_HAQ_narrow<-narrow.all(HAQ_dataset_analysis_mtxIR)
mtxIR_HAQ_narrow

#sub network plot 1
mtxIR.HAQ.d<-mtxIR_HAQ_narrow[[1]]
length(unique(mtxIR.HAQ.d$treatment))
mtxIR_HAQ_subl<-new.medtable(mtxIR.HAQ.d,'treatment')

# 1. Open jpeg file
jpeg("MTX-IR-Health Assessment Questionnaire (HAQ) Network plot with named interventions_AUG2026.jpg", width = 3480, height = 2050, res=300)
nma.plot(mtxIR.HAQ.d[,1],mtxIR.HAQ.d[,3],data = mtxIR.HAQ.d, trtname = as.vector(mtxIR_HAQ_subl$medname),title = "Health Assessment Questionnaire (HAQ) Network plot")
# Close the pdf file
dev.off() 


mtxIR.HAQ.d1<-mtxIR.HAQ.d
mtxIR.HAQ.d1$treatment<-renumber(mtxIR.HAQ.d1,"treatment")

# 1. Open jpeg file
jpeg("MTX-IR-Health Assessment Questionnaire (HAQ) Network plot without named interventions_AUG2026.jpg", width = 3480, height = 2050, res=300)
nma.plot(mtxIR.HAQ.d1[,1],mtxIR.HAQ.d1[,3],data = mtxIR.HAQ.d1,title = "Health Assessment Questionnaire (HAQ) Network plot")
# Close the pdf file
dev.off() 


##############################################################################
#Performing sub network plot 3
library(coda)
library(gemtc)
library(netmeta)
library(tidyr)
View(mtxIR.HAQ.d)
mtxIR_HAQ_subl
col_names<-colnames(mtxIR.HAQ.d)
network_data_HAQ_mtxIR_1<-mtxIR.HAQ.d %>%
  select(col_names)
network_data_HAQ_mtxIR_1$treatment<-renumber(network_data_HAQ_mtxIR_1,"treatment")
network_data_HAQ_mtxIR_1$study<-renumber(network_data_HAQ_mtxIR_1,"study")
network_data_HAQ_mtxIR_1<-network_data_HAQ_mtxIR_1[order(network_data_HAQ_mtxIR_1$study),]
row.names(network_data_HAQ_mtxIR_1)<-NULL

##To download HAQ data
network_data_HAQ_mtxIR_data_download<-network_data_HAQ_mtxIR_1
network_data_HAQ_mtxIR_data_download$treatment<-selected_treatment_name(network_data_HAQ_mtxIR_data_download$treatment,mtxIR_HAQ_subl)
write.csv(network_data_HAQ_mtxIR_data_download,"HAQ MTX-IR dataset_AUG2026.csv",row.names = F)

nma.plot(network_data_HAQ_mtxIR_1$study,network_data_HAQ_mtxIR_1$treatment,data = network_data_HAQ_mtxIR_1, trtname = as.vector(mtxIR_HAQ_subl$medname),title = "Health Assessment Questionnaire (HAQ) Network plot")

network_data_HAQ_mtxIR_1 <- network_data_HAQ_mtxIR_1 %>%
  dplyr:: rename(std.dev = SD, mean = Mean)


p_HAQ_MD<-pairwise(treat=treatment,mean=mean,sd=std.dev,data=network_data_HAQ_mtxIR_1,
                   studlab = study,n=sampleSize,sm="MD")
#colnames(p_HAQ_MD)[c(1,4:5)]<-c("study","diff", "std.err")
p_HAQ_MD <- p_HAQ_MD %>%
  dplyr:: rename("diff"=TE, "std.err"=seTE)

library(magrittr) #for set_colnames
p_HAQ_MD %>%
  dplyr::select(c(study_name,diff, std.err,treat1,treat2)) %>%
  pivot_longer(-study_name,
               names_to = c(".value"),
               names_pattern = "(..)") %>% 
  set_colnames(c("study", "diff", 
                 "std.err", "treatment")) -> data2


clean_and_impute_se <- function(dat) {
  
  dat %>%
    group_by(study) %>%
    group_modify(~ {
      
      df <- .x
      
      # ---------------------------------------------------------
      # STEP 1: Clean multi-arm studies
      # ---------------------------------------------------------
      
      if (nrow(df) > 2) {
        
        # Multi-arm data should consist of two-row comparisons
        if (nrow(df) %% 2 != 0) {
          warning(
            paste0(
              "Study '", .y$study,
              "' has an odd number of rows and was left unchanged."
            )
          )
          
          return(df)
        }
        
        # Count rows where both effect estimate and SE are missing
        na_counts <- df %>%
          group_by(treatment) %>%
          summarise(
            n_na = sum(is.na(diff) & is.na(std.err)),
            .groups = "drop"
          )
        
        # Treatment with the highest number of missing rows
        trt_keep <- na_counts$treatment[
          which.max(na_counts$n_na)
        ]
        
        # Assign every two consecutive rows to a comparison
        df <- df %>%
          mutate(
            .row_order = row_number(),
            .comparison = rep(
              seq_len(nrow(df) / 2),
              each = 2
            )
          )
        
        # Keep only comparisons containing trt_keep
        df <- df %>%
          group_by(.comparison) %>%
          filter(any(treatment == trt_keep)) %>%
          ungroup()
        
        # Keep only the final occurrence of trt_keep
        last_keep <- max(
          df$.row_order[df$treatment == trt_keep]
        )
        
        df <- df %>%
          filter(
            treatment != trt_keep |
              .row_order == last_keep
          ) %>%
          arrange(.row_order) %>%
          select(-.row_order, -.comparison)
      }
      
      # ---------------------------------------------------------
      # STEP 2: Impute missing SE within each study
      # ---------------------------------------------------------
      
      se_obs <- df$std.err[!is.na(df$std.err)]
      
      # Impute when at least one SE is missing and at least two
      # observed SEs are available
      if (
        any(is.na(df$std.err)) &&
        length(se_obs) >= 2
      ) {
        
        se_imp <- sqrt(
          se_obs[1] * se_obs[2] * 0.6
        )
        
        df$std.err[is.na(df$std.err)] <- se_imp
      }
      
      df
    }) %>%
    ungroup()
}

data2 <- clean_and_impute_se(data2)

network_data_HAQ_mtxIR_1<-mtc.network(data.re=data2)
summary(network_data_HAQ_mtxIR_1)
View(network_data_HAQ_mtxIR_1$data.re)

##RUN the model
#network_data_DAS28_mtxIR_1.model  <- mtc.model(network_data_DAS28_mtxIR_1 , likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
set.seed(58)
network_data_HAQ_mtxIR_1.model <- mtc.model(network_data_HAQ_mtxIR_1, likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
#network_data_HAQ_mtxIR_1.model <- mtc.model(network_data_HAQ_mtxIR_1, likelihood='normal', link='identity')

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_HAQ_mtxIR_1.model$inits <- mapply(c, network_data_HAQ_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_HAQ<-mtc.run(network_data_HAQ_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_HAQ)

View(as.data.frame(relative.effect.table(results_mtxIR_HAQ)))
View(as.data.frame(relative.effect.table.2DP2(results_mtxIR_HAQ, comparator = "1")))
#refME_mtxIR_DAS28<-as.data.frame(relative.effect.table(results_mtxIR_DAS28))
refME_mtxIR_HAQ<-as.data.frame(relative.effect.table.2DP2(results_mtxIR_HAQ, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refME_mtxIR_HAQ)<-mapvalues(row.names(refME_mtxIR_HAQ),
                                    from = c(1:(nrow(mtxIR_HAQ_subl))),
                                    to = as.vector(mtxIR_HAQ_subl$medname))
colnames(refME_mtxIR_HAQ)<-mapvalues(colnames(refME_mtxIR_HAQ),
                                   from = c(1:(nrow(mtxIR_HAQ_subl))),
                                   to = as.vector(mtxIR_HAQ_subl$medname))
View(refME_mtxIR_HAQ)

index_diagonal5<-data.frame(index1=1:nrow(refME_mtxIR_HAQ),index2=1:nrow(refME_mtxIR_HAQ))
options(warn=-1)
for(i in 1:nrow(refME_mtxIR_HAQ)){
  refME_mtxIR_HAQ[index_diagonal5$index1[i],index_diagonal5$index2[i]]<-selected_treatment_name(as.numeric(refME_mtxIR_HAQ[row(refME_mtxIR_HAQ)==col(refME_mtxIR_HAQ)])[i],mtxIR_HAQ_subl)
}
options(warn=0)
View(refME_mtxIR_HAQ)
write.csv(refME_mtxIR_HAQ,"HAQ MTX-IR pairwise comparisons_AUG2026.csv", row.names = T)

#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "HAQ-DI")
writeData(wb, "HAQ-DI", refME_mtxIR_HAQ, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)

##GRADE Assessment
nma_estimate_mtxIR_HAQ_dataframe<-nma.estimate.table.MD(results_mtxIR_HAQ)#use results
View(nma_estimate_mtxIR_HAQ_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_HAQ_mtxIR_1.ns<-mtc.nodesplit(network_data_HAQ_mtxIR_1, thin=50,likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 

names(network_data_HAQ_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_HAQ_mtxIR_1.ns)

# Overall summary and plot
summary.network_data_HAQ_mtxIR_1.ns <- summary(network_data_HAQ_mtxIR_1.ns)

#Direct effects
direct.data<-summary.network_data_HAQ_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_HAQ_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_HAQ_mtxIR_1.ns$cons.effect

#P-values
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
P_values_nodesplit_haqmtxIR<-summary.network_data_HAQ_mtxIR_1.ns$p.value
P_values_nodesplit_haqmtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_haqmtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_haqmtxIR2<-P_values_nodesplit_haqmtxIR
P_values_nodesplit_haqmtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_haqmtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_haqmtxIR2$comb<-paste0(P_values_nodesplit_haqmtxIR2$t1,":",P_values_nodesplit_haqmtxIR2$t2)

direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(direct.data[,c(3,4,5)],2),round(indirect.data[,c(3,4,5)],2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                -round(direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)
##Acquiring both direct and indirect estimate from gemtc
grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_HAQ_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)
###Adding ROB data
HAQ_mtxIR<-unique(mtxIR.HAQ.d$study_name)
data_qq_mtxIR<-mtxIR.HAQ.d
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%HAQ_mtxIR,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<- current_assessment_for_studies_mtxIR %>%
  select(c(col_names))

list_appraisal_mtxIR3<-list()
for(i in 1:length(HAQ_mtxIR)){
  temp11<-data_qq_mtxIR[data_qq_mtxIR$study_name==HAQ_mtxIR[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==HAQ_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal_mtxIR3[[i]]<-temp13
}
HAQ_mtxIR_withROB<-do.call(rbind,list_appraisal_mtxIR3)
View(HAQ_mtxIR_withROB)
HAQ_mtxIR_withROB2<-HAQ_mtxIR_withROB%>%
  select(c("study","study_name", "treatment", "SD", "sampleSize",
           "Mean","Timepoint","OVERALL..ALL.OTHER.OUTCOMES"))
View(HAQ_mtxIR_withROB2)


###Adding Indirectness data
study_unique<-unique(HAQ_mtxIR_withROB2$study_name)
data_nn<-HAQ_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp[!is.na(resolved_indirectness_consensus_temp$`Final Rating`),]
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]
row.names(current_indirect_assessment_for_studies)<-NULL
# current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:11),],current_indirect_assessment_for_studies[11,],current_indirect_assessment_for_studies[c(12:33),])
# current_indirect_assessment_for_studies$intervention[12]<-"Baricitinib (4 mg)"
row.names(current_indirect_assessment_for_studies)<-NULL
list_appraisal3<-list()
for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal3[[i]]<-temp13
}
HAQ_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal3)
HAQ_mtxIR_withROB2_indrectness <- HAQ_mtxIR_withROB2_indrectness %>%
  dplyr::rename(
    ROB = OVERALL..ALL.OTHER.OUTCOMES,
    Indirectness = `Final Rating`
  )


HAQ_mtxIR_withROB2_indrectness$Indirectness<-tolower(HAQ_mtxIR_withROB2_indrectness$Indirectness)

View(HAQ_mtxIR_withROB2_indrectness)

HAQ_mtxIR_data_new_withROB2_indrectness2<-HAQ_mtxIR_withROB2_indrectness
HAQ_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(HAQ_mtxIR_data_new_withROB2_indrectness2,"treatment")
HAQ_mtxIR_data_new_withROB2_indrectness2$ROB<-ifelse(HAQ_mtxIR_data_new_withROB2_indrectness2$ROB=="unclear","some concerns",HAQ_mtxIR_data_new_withROB2_indrectness2$ROB)

p_HAQ_MD2<-p_HAQ_MD
#colnames(p_HAQ_MD2)[c(1,4:5)]<-c("studlab","TE","seTE")
p_HAQ_MD2<- p_HAQ_MD2 %>%
  dplyr:: rename(TE="diff", seTE ="std.err")

net1 <- netmeta(p_HAQ_MD2)
dat.net<-netsplit(net1)
View(dat.net$indirect.common)
View(cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)],dat.net$compare.common[,c(2:5)]))
View(dat.net)
dat.net.contrib<-netcontrib(net1, method = "r")
View(dat.net.contrib$random)
View(dat.net.contrib$common)

#direct_evidence_new<-HAQ_mtxIR_data_new_withROB2_indrectness2_details$hatMatrix$side
direct_evidence_new<-cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)])
View(direct_evidence_new)
View(HAQ_mtxIR_data_new_withROB2_indrectness2_details$contributionMatrix)
#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")


#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]

##Getting direct evidence for intervention in only multi-arm trials
zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
colnames(zzz)<-c("treat1", "treat2")
grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]


grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_HAQ_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_HAQ_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))

#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
HAQ_output_df2_list<-list()
I_square_list<-list()



for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-HAQ_mtxIR_data_new_withROB2_indrectness2[HAQ_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select(study, Mean.1, SD.1, Mean.2, SD.2 )
  temp11.1[,c("SD.1","SD.2")]<-temp11.1[,c("SD.1","SD.2")]/(sqrt(temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c("sampleSize.1" , "sampleSize.2")]))
  
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","se.t","x.c","se.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    model.file<-"Cont_direct_onestudy.txt"
    n.chains=3
    set.seed(12451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("diff")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(MD=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    model.file<-"Cont_direct_multiplestudy_RE.txt"
    n.chains=3
    set.seed(11451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c, "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("diff")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_dataHAQ_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select(study, sampleSize.1, Mean.1, SD.1, sampleSize.2, Mean.2, SD.2 )
  #nma_dataDAS28_3[,c("SD.1","SD.2")]<-nma_dataDAS28_3[,c("SD.1","SD.2")]/(sqrt(nma_dataDAS28_3[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(2,5)]))
  nma_dataHAQ_4<-nma_dataHAQ_3
  colnames(nma_dataHAQ_4)<-c("NS","n.t","mean.t","sd.t","n.c","mean.c","sd.c")
  row.names(nma_dataHAQ_4)<-NULL
  m.cont<-metacont(n.e = n.t, 
                   mean.e = mean.t, 
                   sd.e = sd.t, 
                   n.c = n.c, 
                   mean.c = mean.c, 
                   sd.c =sd.c,
                   data = nma_dataHAQ_4,
                   sm = "SMD",
                   method.tau = "PM",
                   random = TRUE,
                   method.random.ci = TRUE,
                   title = "HAQ"
                   
  )
  
  summary_m.cont<-summary(m.cont)
  contribution_weight<-as.character(round((summary_m.cont$w.random/sum(summary_m.cont$w.random))*100,2))
  pooled_estimate2<-data.frame(Mean=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.cont$tau2,4),
                          Lower=round(summary_m.cont$lower.tau2,4),
                          Upper=round(summary_m.cont$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.cont$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.cont$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.cont$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.cont$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  Mean_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1,temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(5,12)])
  base_data <- tibble(mean  = Mean_individual_studies$MD, 
                      lower = Mean_individual_studies$Lower,
                      upper = Mean_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      treatment1 = paste0(as.character(temp11.2$x.c),"      ","(",as.character(temp11.2$sampleSize.2),")"),
                      treatment2 = paste0(as.character(temp11.2$x.t),"      ","(",as.character(temp11.2$sampleSize.1),")"),
                      MD = paste(format(round(Mean_individual_studies$MD,2),nsmall=2),paste("(",format(round(Mean_individual_studies$Lower,2),nsmall=2),", ",format(round(Mean_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$Mean, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   MD = paste(format(round(pooled_estimate2$Mean,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("HAQ (final)", "Study"),
                   treatment1 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.2),mtxIR_HAQ_subl)),
                   treatment2 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.1),mtxIR_HAQ_subl)),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   MD = c("", "MD"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  HAQ_output_df <- bind_rows(header,
                             base_data,
                             empty_row,
                             pooled)
  
  
  
  
  HAQ_output_df2<-rbind(HAQ_output_df,HAQ_output_df[c(nrow(HAQ_output_df)),])
  HAQ_output_df2[nrow(HAQ_output_df),1]<-I2
  HAQ_output_df2[nrow(HAQ_output_df),2]<-NA
  HAQ_output_df2[nrow(HAQ_output_df),c(6:11)]<-c(NA,NA,NA,NA,NA,NA,NA)
  HAQ_output_df2_list[[i]]<-HAQ_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-HAQ_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, MD,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               clip = c(-1,1), 
               xlog = F, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_HAQ_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_HAQ_subl))
    )
  
  
  
}#END


View(HAQ_output_df2_list[[19]])
n_studies<-c()
n_participants<-c()

extract_num<-function(data){
  return(str_extract(data,"(?<=\\()\\d+(?=\\))"))
}
countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}




for(i in 1:length(graph_list)){
  df<-HAQ_output_df2_list[[i]][-c(1,2,countdown(length(HAQ_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df,1,extract_num))))
}

length(graph_list)
graph_list[[18]]
Grade_assessment_dataset_direct_HAQ<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_HAQ_subl),
                                                selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_HAQ_subl),
                                                Outcome="HAQ",
                                                `Effect measure`="MD",
                                                `Number of studies`=n_studies,
                                                `Number of participant`=n_participants,
                                                ROB="",Indirectness="",Imprecision="",
                                                Inconsistency="",
                                                Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_HAQ)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                       "Intevention",
                                                                       "Effect measure",
                                                                       "Number of studies",
                                                                       "Number of participant",
                                                                       "Risk of Bias",
                                                                       "Publication Bias",
                                                                       "Certainty of evidence (ignores imprecision)")

View(Grade_assessment_dataset_direct_HAQ)
#Saving grade assessment details for HAQ BIOLOGICS as an object of lists
grade_details_HAQ_mtxIR<-list(graph_list,HAQ_output_df2_list,Grade_assessment_dataset_direct_HAQ,I_square_list)
#saveRDS(grade_details_HAQ_bio,file="List of HAQ Biologics data components for GRADE_OCT31.rds")
saveRDS(grade_details_HAQ_mtxIR,file="List of HAQ MTX-IR data components for GRADE_AUG2026.rds")
#
haq_output_one_study_list <- list()
haq_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

n_studies<-c()
for(i in 1:length(graph_list)){
  df<-HAQ_output_df2_list[[i]][-c(1,2,countdown(length(HAQ_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
}

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
haq_output_one_study_list  <- HAQ_output_df2_list[one_idx]
haq_output_more_than_2_list <- HAQ_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_HAQ[
    Grade_assessment_dataset_direct_HAQ$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_HAQ[
    Grade_assessment_dataset_direct_HAQ$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_haq_mtxIR_one_study <- list(
  graph_list_one_study_list,
  haq_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_haq_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  haq_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_haq_mtxIR_one_study,
  file = "List of HAQ MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_haq_mtxIR_more_than_2,
  file = "List of HAQ MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)


#####################################
#####################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_HAQ_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct_HAQ
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_HAQ$Treatment, 
                                                Grade_assessment_dataset_direct_HAQ$Intevention), 
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(HAQ_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","MD"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 

glen_direct_evidence_rating_haq<-read.csv("Glen_2026-09-03 22_05_HAQ-MTXIR_complete.csv")
glen_direct_evidence_rating_haq2<-glen_direct_evidence_rating_haq[,-c(13:15)]
colnames(glen_direct_evidence_rating_haq2)[1:2]<-colnames(grade_table4)[1:2]

glen_direct_evidence_rating_haq2
View(merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_haq2,
           by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T))
glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_haq2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"HAQ"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"MD"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR_HAQ_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR_HAQ_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}

##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for HAQ_MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for HAQ_MTX-IR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR HAQ OUTCOME
initial_gradeTable_haq<-Temp_GRADE5

estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)
flip_to_smaller_to_larger <- function(vector) {
  # Extract values before and after ":" and split them into two vectors
  before_colon <- as.numeric(sub(":.*", "", vector))
  after_colon <- as.numeric(sub(".*:", "", vector))
  
  # Determine which values are smaller
  smaller_values <- before_colon > after_colon
  
  # Create a new vector with flipped values
  flipped_vector <- ifelse(smaller_values, paste(after_colon, ":", before_colon, sep = ""), vector)
  
  return(flipped_vector)
}

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_haq2<-initial_gradeTable_haq[match(all_nma_estimates,initial_gradeTable_haq$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()

#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_haq2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_haq2.1<-cbind(initial_gradeTable_haq2[,c(1:14)],first_order_loop_contri,
                                 initial_gradeTable_haq2[,-c(1:14)])

initial_gradeTable_haq2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_haq2.1,
                                                            all_nma_estimates,
                                                            estimate_contribution_matrix,
                                                            glen_direct_evidence_rating6_onlyDirect,
                                                            top_n = 2))
ncol(initial_gradeTable_haq2.1)
View(initial_gradeTable_haq2.1)

# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_haq2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_haq2.1<-result_df
initial_gradeTable_haq2.1$combI<-str_replace(initial_gradeTable_haq2.1$comb, ":", " vs ")

initial_gradeTable_haq2.1<-process_dataset(initial_gradeTable_haq2.1, 2)
initial_gradeTable_haq3<-initial_gradeTable_haq2.1[match(initial_gradeTable_haq$comb,initial_gradeTable_haq2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for HAQ_MTX-IR outcome.csv",row.names = T)


###########Assessing Intransitivity
##Here, since all the indirectness assessments are low, intransitivity will be not serious for all.
##Else manually evaluate intransitivity of each direct evidence contributing to the indirect evidence.
check_vector <- function(vector) {
  intransitivity_rating <- ifelse("Serious" %in% vector, "Serious",
                                  ifelse("Some Concerns" %in% vector, "Some Concerns", "Not Serious"))
  degree <- ifelse((intransitivity_rating == "Not Serious")|(intransitivity_rating == "Some Concerns"), 0, -1)
  
  return(c(degree, intransitivity_rating))
}


sub_data_intransitivity <- initial_gradeTable_haq3[initial_gradeTable_haq3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_haq3[initial_gradeTable_haq3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_haq3.1<-cbind(initial_gradeTable_haq3[,c(1:23)],
                                 Intransitivity="",
                                 `Certainty of evidence (ignores imprecision)-Indirect`="",
                                 initial_gradeTable_haq3[,-c(1:23)])
initial_gradeTable_haq3.1[initial_gradeTable_haq3.1$Indirect!="-",c("Intransitivity",
                                                                    "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                          "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_haq4<-initial_gradeTable_haq3.1

# ###Applying the P-values, when no p-values is available because no both direct and indirect evidence available
# #P_values_nodesplit_haqmtxIR2
# grade_temp<-cbind(initial_gradeTable_haq4[,c(1:25)],DirectPerContri="-",IndirectPerContri="-",initial_gradeTable_haq4[,-c(1:25)])
# grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])
# 
# ##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
# grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
# #When only direct or only indirect estimate is present
# grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
#                                  grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
#                                  ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
#                                         grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))

###Applying the P-values
P_values_nodesplit_haqmtxIR2
grade_temp<-cbind(initial_gradeTable_haq4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_haq4[,-c(1:25)])
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])
drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_haqmtxIR3 <- P_values_nodesplit_haqmtxIR2[match(drt$comb, P_values_nodesplit_haqmtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_haqmtxIR3$p,2)
#grade_temp2<-grade_temp

#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)

subHAQSOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatHAQ<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subHAQSOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(nrow(subcontriMatHAQ))
IndirPer <- numeric(nrow(subcontriMatHAQ))

# Loop to calculate dirPer and IndirPer
for (i in 1:nrow(subcontriMatHAQ)) {
  dirPer[i] <- subcontriMatHAQ[i, colnames(subcontriMatHAQ) == row.names(subcontriMatHAQ)[i]]
  IndirPer[i] <- 1-dirPer[i]
  #IndirPer[i] <- sum(subcontriMatHAQ[i, colnames(subcontriMatHAQ) != row.names(subcontriMatHAQ)[i]])
}

# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")

grade_temp2[match(row.names(subcontriMatHAQ), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(row.names(subcontriMatHAQ), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)

##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
#When both direct and indirect are present

fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg



###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="-",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho

initial_gradeTable_haq5<-grade_temp2.2
initial_gradeTable_haq5<-initial_gradeTable_haq5[,-c(30,33)]
colnames(initial_gradeTable_haq5)[c(29,30)]<-c("Initial NMA evidence rating",
                                               "P-Value")
###After incoherence evaluation
initial_gradeTable_haq5.1<-cbind(initial_gradeTable_haq5[,c(1:31)],
                                 qoeAfterIncohenrence="",
                                 initial_gradeTable_haq5[,-c(1:31)])
initial_gradeTable_haq5.1[(initial_gradeTable_haq5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_haq5.1[(initial_gradeTable_haq5.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_haq5.1[initial_gradeTable_haq5.1$Incoherence!="-",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_haq5.1[(initial_gradeTable_haq5.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence


#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_haq5.1[(initial_gradeTable_haq5.1$Direct!="-")&(initial_gradeTable_haq5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_haq5.1[(initial_gradeTable_haq5.1$Direct!="-")&(initial_gradeTable_haq5.1$Indirect=="-"),"Direct"]

###Rating NMA after imprecision
initial_gradeTable_haq6<-cbind(initial_gradeTable_haq5.1,
                               p.estimate=apply(as.data.frame(initial_gradeTable_haq5.1$NMA),1,function(x){
                                 as.numeric(sub(" .*", "", x))
                               }),
                               t(apply(as.data.frame(initial_gradeTable_haq5.1$NMA),1,function(x){
                                 as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                 
                               })))

colnames(initial_gradeTable_haq6)[c(36,37)]<-c("l.bound","u.bound")

to_higher_lower_text <- function(e1) {
  sign_val <- sign(e1)
  abs_val <- abs(round(e1, 2))
  text <- ifelse(sign_val < 0, "lower", "higher")
  return(paste(abs_val, text))
}


####Running Continuous NMA model
# load data for MTC
##Getting baseline risk using the placebo intervention
# HAQ_bio_withROB2_indrectness$treatment<-renumber(HAQ_bio_withROB2_indrectness,"treatment")
# HAQ_bio_withROB2_indrectness$study<-renumber(HAQ_bio_withROB2_indrectness,"study")
HAQ_data<-HAQ_mtxIR_data_new_withROB2_indrectness2
HAQ_data<-HAQ_data[order(HAQ_data$study),c(1:6)]


HAQ_data_wide2<-reshape(number_of_arms(HAQ_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(HAQ_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(HAQ_data_wide2[,c(3,7,11)],1,xcv)

HAQ_data_wide3<-cbind(HAQ_data_wide2,na=numberOfArms)
View(HAQ_data_wide3)
#finding standard errors
HAQ_data_wide3$SD.1<-HAQ_data_wide3$SD.1/sqrt(HAQ_data_wide3$sampleSize.1)
HAQ_data_wide3$SD.2<-HAQ_data_wide3$SD.2/sqrt(HAQ_data_wide3$sampleSize.2)
HAQ_data_wide3$SD.3<-HAQ_data_wide3$SD.3/sqrt(HAQ_data_wide3$sampleSize.3)
#HAQ_data_wide3$SD.4<-HAQ_data_wide3$SD.4/sqrt(HAQ_data_wide3$sampleSize.4)
colnames(HAQ_data_wide3)<-c("study", "study_name", "t.1","se.1","n.1", "x.1",
                            "t.2","se.2","n.2","x.2","t.3","se.3","n.3","x.3", "na")
                            #"t.4","se.4","n.4","x.4","na" )
MTCData_cont <- as.data.frame(HAQ_data_wide3)
x <- data.matrix(MTCData_cont[,c(6,10,14)])
n <- data.matrix(MTCData_cont[,c(5,9,13)])
t <- data.matrix(apply(MTCData_cont[,c("t.1","t.2","t.3")],c(1,2),as.numeric))
se <- data.matrix(MTCData_cont[,c(4,8,12)])
na <- data.matrix(MTCData_cont[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- x[t[,'t.1']==1,"x.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
sePBO <- se[t[,'t.1']==1,"se.1"]

# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
#model.file4<-"Cont_multiarm_RE.txt"
model.file4<-"Cont_new_multiarm.txt"
n.chains=3
model.inits <- initv1 

cont.model3<-jags.model(file=model.file4,
                        data= list("x"=x,"t"=t,"se"=se,"na" = na[,1], "NT" = NT, "NS" = NS,
                                   "NTpbo"= NTpbo,"xPBO"=xPBO,"sePBO"=sePBO),
                        n.chains=n.chains,
                        inits = model.inits,
                        n.adapt=5000
)

update(object=cont.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("diff","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( cont.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]

#Converting estimates relative to baseline to absolutes
estimates_relative_to_baseline<-initial_gradeTable_haq6[1:(NT-1),c(1,2,35:37)]
typeof(estimates_relative_to_baseline)

absolutes_from_baseline<-cbind(baseline_risk =round(baseline_risk,2),round(baseline_risk,2)+estimates_relative_to_baseline[,c(3:5)])
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(round(baseline_risk,2),estimates_relative_to_baseline_combine[,7]))
colnames(all_baseline)<-"base"

com.list<-list()
unique_treat2<-unique(initial_gradeTable_haq6$treat2)
d.sub<-initial_gradeTable_haq6[initial_gradeTable_haq6$treat2==unique_treat2[1],c(1,2,35:37)]
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_haq6[initial_gradeTable_haq6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,base+d.sub[,c(3:5)])
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe<-com.dataframe[c(1:nrow(initial_gradeTable_haq6)),]
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9])

higher_lower_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"lower"
  }else{
    text<-"higher"
  }
  return(paste(abs(round(e1-e2,2)),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,13])))


prefix.assumedrisk <- ' Mean '
prefix.corrrisk<-' Mean '
prefix.md<- " MD "
suffix.assumedrisk <- ''
suffix.corrrisk <- ''
time<-1


result_com.dataframe_withtext2<-cbind(com.dataframe_withtext[,1:2],
                                      paste(prefix.assumedrisk,com.dataframe_withtext[,6], suffix.assumedrisk, sep=''),
                                      paste(prefix.corrrisk, com.dataframe_withtext[,7], suffix.corrrisk, ' (', com.dataframe_withtext[,8], ' to ',com.dataframe_withtext[,9],')', sep=''),
                                      paste(prefix.md, com.dataframe_withtext[,10], ' (', com.dataframe_withtext[,11], ' to ',com.dataframe_withtext[,12],')', sep=''),
                                      # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                      com.dataframe_withtext[,c(3:5,10:12)]
)

colnames(result_com.dataframe_withtext2)[c(3:11)]<-c("Assumed mean (comparator)",	"Corresponding mean (intervention)",
                                                     "Mean difference","p.estimate","l.bound","u.bound","","","")

View(cbind(initial_gradeTable_haq6[,-c(35:37)],
           result_com.dataframe_withtext2[,-c(1,2)]))

initial_gradeTable_haq7<-cbind(initial_gradeTable_haq6[,-c(35:37)],
                               result_com.dataframe_withtext2[,-c(1,2)],pair="")
#colnames(initial_gradeTable_haq7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_haq8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_haq7$treat1,mtxIR_HAQ_subl)),
                               `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_haq7$treat2,mtxIR_HAQ_subl)),
                               initial_gradeTable_haq7)
colnames(initial_gradeTable_haq8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_haq8<-initial_gradeTable_haq8[,-c(46)]
write.csv(initial_gradeTable_haq8,"HAQ MTX-IR SoF table without imprecision assessment AUG2026.csv",row.names = F)



effect.size.rating_list<-lapply(split(initial_gradeTable_haq8, seq_len(nrow(initial_gradeTable_haq8))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 0.21, large.effect = 0.5, framework = T)
})



as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_haq9<-initial_gradeTable_haq8

initial_gradeTable_haq9<-cbind(initial_gradeTable_haq9[,c(1:34,37:45,35:36)])

initial_gradeTable_haq9.1<-cbind(initial_gradeTable_haq9[,c(1:37)],
                                 Small.effect.threshold= 0.2,
                                 Large.effect.threshold= 0.5,
                                 Framework = "NULL",
                                 Imprecision.1 = "",
                                 Certainty.of.evidence = "",
                                 initial_gradeTable_haq9[,c(44,45)])



initial_gradeTable_haq9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_HAQ_mtxIR<-initial_gradeTable_haq9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_HAQ_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_HAQ_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_HAQ_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_HAQ_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_HAQ_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_HAQ_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_HAQ_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_HAQ_mtxIR[is.na(detailed_SOF_HAQ_mtxIR$Certainty.of.evidence),])
detailed_SOF_HAQ_mtxIR[is.na(detailed_SOF_HAQ_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

detailed_SOF_HAQ_mtxIR$comb<-initial_gradeTable_haq7$comb

write.csv(detailed_SOF_HAQ_mtxIR,"HAQ_MTX-IR_Detail_SoF_AUG2026.csv", row.names = F)


################################################################################################
################################################################################################
##############PAIN Analysis###################################
###################################################
############################
#######Graphical data for Pain and errors


dataset_Pain<-final_mtx_irDataset_actual%>% 
  select(matches("^Pain\\..*")) %>%
  select(-contains("VAS"))%>%
  select(-contains("Reported"))

colnames(dataset_Pain)
dataset_Pain[dataset_Pain=="NA"]<-NA
View(dataset_Pain)

#12 months
actuals<-dataset_Pain%>%
  select(matches("^Pain\\.12\\.months..actual*"))
View(actuals)
colnames(actuals)
library(tibble)# library to add column
actuals<-add_column(actuals, "Pain.12.months..actual.value.Mean_edit" = actuals$Pain.12.months..actual.value.Mean, .after = "Pain.12.months..actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
index<-which(is.na(actuals$Pain.12.months..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){#evaluate when entire dataset is NA
  actuals<-add_column(actuals, "Pain.12.months..actual.value.SD_edit" = actuals$Pain.12.months..actual.value.SD, .after = "Pain.12.months..actual.value.SD")
  changeFbaseline<-dataset_Pain%>%
    select(matches("^Pain\\.12.months\\..change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Pain.12.months..change.from.baseline.Mean_edit" =changeFbaseline$Pain.12.months..change.from.baseline.Mean, .after = "Pain.12.months..change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Pain.12.months..change.from.baseline.SD_edit" = changeFbaseline$Pain.12.months..change.from.baseline.SD, .after = "Pain.12.months..change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Pain.12.months..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Pain.12.months..actual.value.Mean_edit"]<-sub_actuals[index2,"Pain.12.months..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  actuals<-add_column(actuals, "Pain.12.months..actual.value.SD_edit" = actuals$Pain.12.months..actual.value.SD, .after = "Pain.12.months..actual.value.SD")
  
  #for missing variability
  index3<-which(is.na(actuals$Pain.12.months..actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Pain.12.months..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Pain.12.months..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Pain.12.months..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Pain.12.months..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Pain.12.months..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Pain.12.months..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Pain.12.months..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.12.months..actual.value.Upper.CI))
      temp_sd[index7,"Pain.12.months..actual.value.SD_edit"]<-((temp_sd[index7,"Pain.12.months..actual.value.Upper.CI"]-temp_sd[index7,"Pain.12.months..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.12.months..actual.value.N"])
      sub_actuals_sd[index6,] <-temp_sd
      actuals[index3,]<-sub_actuals_sd
    } 
    if((length(which(is.na(sub_actuals_sd[,"Pain.12.months..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Pain.12.months..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.12.months..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.12.months..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.12.months..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.12.months..actual.value.N"][which.max(temp_sd2[index9,"Pain.12.months..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.12.months..actual.value.N"][which.min(temp_sd2[index9,"Pain.12.months..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.12.months..actual.value.p.value"]))
      values<-as.numeric(c(p.val, Nrx, Npbo, meanhigh, meanlow))
      temp_sd2[index9,"Pain.12.months..actual.value.SD_edit"]<- exactp.sd(p.val=values[1], Nrx=values[2], Npbo=values[3], meanhigh=values[4], meanlow=values[5])[[2]]
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  
  ####12 months Change from baseline
  changeFbaseline<-dataset_Pain%>%
    select(matches("^Pain\\.12.months\\..change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Pain.12.months..change.from.baseline.Mean_edit" =changeFbaseline$Pain.12.months..change.from.baseline.Mean , .after = "Pain.12.months..change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Pain.12.months..change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Pain.12.months..change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Pain.12.months..change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Pain.12.months..change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Pain.12.months..change.from.baseline.SD_edit" = changeFbaseline$Pain.12.months..change.from.baseline.SD, .after = "Pain.12.months..change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Pain.12.months..change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Pain.12.months..change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Pain.12.months..change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Pain.12.months..change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.12.months..change.from.baseline.Upper.CI))
      temp_sd[index7,"Pain.12.months..change.from.baseline.SD_edit"]<-((as.numeric(temp_sd[index7,"Pain.12.months..change.from.baseline.Upper.CI"])-as.numeric(temp_sd[index7,"Pain.12.months..change.from.baseline.Lower.CI"]))/3.92)*sqrt(as.numeric(temp_sd[index7,"Pain.12.months..change.from.baseline.N"]))
      sub_changeFbaseline_sd[index6,] <-temp_sd
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Pain.12.months..change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.12.months..change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.12.months..change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.12.months..change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.12.months..change.from.baseline.N"][which.max(temp_sd2[index9,"Pain.12.months..change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.12.months..change.from.baseline.N"][which.min(temp_sd2[index9,"Pain.12.months..change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.12.months..change.from.baseline.p.value"]))
      values<-as.numeric(c(p.val, Nrx, Npbo, meanhigh, meanlow))
      temp_sd2[index9,"Pain.12.months..change.from.baseline.SD_edit"]<- exactp.sd(p.val=values[1], Nrx=values[2], Npbo=values[3], meanhigh=values[4], meanlow=values[5])[[2]]
      # temp_sd2[index9,"DAS.12.months..change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.12.months..change.from.baseline.p.value"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.Mean"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}


##Replacing processed Pain outcome-12months
Pain_12months_continuous_outcome<-cbind(actuals,changeFbaseline)

####3 months Actuals
dataset_Pain_3months<-dataset_Pain%>%
  select(matches("^Pain\\.3\\.months*"))
dataset_Pain_3months[dataset_Pain_3months=="NA"]<-NA
View(dataset_Pain_3months)
actuals<-dataset_Pain_3months%>%
  select(matches("^Pain\\.3\\.months...actual*"))
View(actuals)
colnames(actuals)
#for missing means
actuals<-add_column(actuals, "Pain.3.months....actual.value.Mean_edit" = actuals$Pain.3.months...actual.value.Mean, .after = "Pain.3.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Pain.3.months...actual.value.Mean))

if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Pain.3.months...actual.value.SD_edit" = actuals$Pain.3.months...actual.value.SD, .after = "Pain.3.months...actual.value.SD")
  changeFbaseline<-dataset_Pain_3months%>%
    select(matches("^Pain\\.3.months\\...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Pain.3.months...change.from.baseline.Mean_edit" =changeFbaseline$Pain.3.months...change.from.baseline.Mean, .after = "Pain.3.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Pain.3.months...change.from.baseline.SD_edit" = changeFbaseline$Pain.3.months...change.from.baseline.SD, .after = "Pain.3.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Pain.3.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Pain.3.months...actual.value.Mean_edit"]<-sub_actuals[index2,"Pain.3.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Pain.3.months...actual.value.SD_edit" = actuals$Pain.3.months...actual.value.SD, .after = "Pain.3.months...actual.value.SD")
  index3<-which(is.na(actuals$Pain.3.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Pain.3.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Pain.3.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Pain.3.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Pain.3.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Pain.3.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Pain.3.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Pain.3.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.3.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.3.months...actual.value.SD_edit"]<-((temp_sd[index7,"Pain.3.months...actual.value.Upper.CI"]-temp_sd[index7,"Pain.3.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.3.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"Pain.3.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Pain.3.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.3.months...actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.3.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.3.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.3.months...actual.value.N"][which.max(temp_sd2[index9,"Pain.3.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.3.months...actual.value.N"][which.min(temp_sd2[index9,"Pain.3.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.3.months...actual.value.p.value"]))
      temp_sd2[index9,"Pain.3.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####3 months Change from baseline
  changeFbaseline<-dataset_Pain_3months%>%
    select(matches("^Pain\\.3.months\\...change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Pain.3.months...change.from.baseline.Mean_edit" =changeFbaseline$Pain.3.months...change.from.baseline.Mean, .after = "Pain.3.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Pain.3.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Pain.3.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Pain.3.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Pain.3.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Pain.3.months...change.from.baseline.SD_edit" = changeFbaseline$Pain.3.months...change.from.baseline.SD, .after = "Pain.3.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Pain.3.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Pain.3.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Pain.3.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Pain.3.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.3.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.3.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"Pain.3.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"Pain.3.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.3.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Pain.3.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.3.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.3.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.3.months...change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.3.months...change.from.baseline.N"][which.max(temp_sd2[index9,"Pain.3.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.3.months...change.from.baseline.N"][which.min(temp_sd2[index9,"Pain.3.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.3.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"Pain.3.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}




##Replacing processed Pain outcome-3months
Pain_3months_continuous_outcome<-cbind(actuals,changeFbaseline)


#################Pain 6 months
####6 months Actuals
dataset_Pain_6months<-dataset_Pain%>%
  select(matches("^Pain\\.6\\.months*"))
dataset_Pain_6months[dataset_Pain_6months=="NA"]<-NA
View(dataset_Pain_6months)
actuals<-dataset_Pain_6months%>%
  select(matches("^Pain\\.6\\.months...actual*"))
View(actuals)
colnames(actuals)

#for missing means
actuals<-add_column(actuals, "Pain.6.months...actual.value.Mean_edit" = actuals$Pain.6.months...actual.value.Mean, .after = "Pain.6.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Pain.6.months...actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Pain.6.months...actual.value.SD_edit" = actuals$Pain.6.months...actual.value.SD, .after = "Pain.6.months...actual.value.SD")
  changeFbaseline<-dataset_Pain_6months%>%
    select(matches("^Pain\\.6.months\\...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Pain.6.months...change.from.baseline.Mean_edit" =changeFbaseline$Pain.6.months...change.from.baseline.Mean, .after = "Pain.6.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Pain.6.months...change.from.baseline.SD_edit" = changeFbaseline$Pain.6.months...change.from.baseline.SD, .after = "Pain.6.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Pain.6.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Pain.6.months...actual.value.Mean_edit"]<-sub_actuals[index2,"Pain.6.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Pain.6.months...actual.value.SD_edit" = actuals$Pain.6.months...actual.value.SD, .after = "Pain.6.months...actual.value.SD")
  index3<-which(is.na(actuals$Pain.6.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Pain.6.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Pain.6.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Pain.6.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Pain.6.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Pain.6.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Pain.6.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Pain.6.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.6.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.6.months...actual.value.SD_edit"]<-((temp_sd[index7,"Pain.6.months...actual.value.Upper.CI"]-temp_sd[index7,"Pain.6.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.6.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"Pain.6.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Pain.6.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.6.months...actual.value.p.valu))
      meanhigh<-max(temp_sd2[index9,"Pain.6.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.6.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.6.months...actual.value.N"][which.max(temp_sd2[index9,"Pain.6.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.6.months...actual.value.N"][which.min(temp_sd2[index9,"Pain.6.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.6.months...actual.value.p.value"]))
      temp_sd2[index9,"Pain.6.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####6 months Change from baseline
  changeFbaseline<-dataset_Pain_6months%>%
    select(matches("^Pain\\.6.months\\...change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #View(changeFbaseline)
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Pain.6.months...change.from.baseline.Mean_edit" =changeFbaseline$Pain.6.months...change.from.baseline.Mean, .after = "Pain.6.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Pain.6.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Pain.6.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Pain.6.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Pain.6.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Pain.6.months...change.from.baseline.SD_edit" = changeFbaseline$Pain.6.months...change.from.baseline.SD, .after = "Pain.6.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Pain.6.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Pain.6.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Pain.6.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Pain.6.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.6.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.6.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"Pain.6.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"Pain.6.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.6.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Pain.6.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.6.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.6.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.6.months...change.from.baseline.Mean.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.6.months...change.from.baseline.N"][which.max(temp_sd2[index9,"Pain.6.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.6.months...change.from.baseline.N"][which.min(temp_sd2[index9,"Pain.6.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.6.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"Pain.6.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
  
  
}

##Replacing processed Pain outcome-6months
Pain_6months_continuous_outcome<-cbind(actuals,changeFbaseline)


options(warn = defaultW)
################Pain baseline
####Pain baseline Actuals
dataset_PainBaseline<-dataset_Pain%>%
  select(matches("^Pain\\.Baseline*"))
dataset_PainBaseline[dataset_PainBaseline=="NA"]<-NA
View(dataset_PainBaseline)
actuals<-dataset_PainBaseline
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "Pain.Baseline.Mean_edit" = actuals$Pain.Baseline.Mean, .after = "Pain.Baseline.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Pain.Baseline.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Pain.Baseline.SD_edit" = actuals$Pain.Baseline.SD, .after = "Pain.Baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Pain.Baseline.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Pain.Baseline.Mean_edit"]<-sub_actuals[index2,"Pain.Baseline.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Pain.Baseline.SD_edit" = actuals$Pain.Baseline.SD, .after = "Pain.Baseline.SD")
  index3<-which(is.na(actuals$Pain.Baseline.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Pain.Baseline.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Pain.Baseline.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Pain.Baseline.SD_edit"]<-(sub_actuals_sd[index4,"Pain.Baseline.SE"])*sqrt(sub_actuals_sd[index4,"Pain.Baseline.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Pain.Baseline.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Pain.Baseline.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.Baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.Baseline.SD_edit"]<-((temp_sd[index7,"Pain.Baseline.Upper.CI"]-temp_sd[index7,"Pain.Baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.Baseline.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"Pain.Baseline.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Pain.Baseline.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.Baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.Baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.Baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.Baseline.N"][which.max(temp_sd2[index9,"Pain.Baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.Baseline.N"][which.min(temp_sd2[index9,"Pain.Baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.Baseline.p.value"]))
      temp_sd2[index9,"Pain.Baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
}


options(warn = defaultW)
##Replacing processed Pain outcome-Baseline
PainBaseline_continuous_outcome<-actuals

#################Pain End of trial
####End of trial Actuals
dataset_Pain_Endoftrial<-dataset_Pain%>%
  select(matches("^Pain\\.End.of.trial*"))
dataset_Pain_Endoftrial[dataset_Pain_Endoftrial=="NA"]<-NA
View(dataset_Pain_Endoftrial)
actuals<-dataset_Pain_Endoftrial%>%
  select(matches("^Pain\\.End.of.trial\\..actual*"))
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "Pain.End.of.trial..actual.value.Mean_edit" = actuals$Pain.End.of.trial..actual.value.Mean, .after = "Pain.End.of.trial..actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
index<-which(is.na(actuals$Pain.End.of.trial..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Pain.End.of.trial..actual.value.SD_edit" = actuals$Pain.End.of.trial..actual.value.SD, .after = "Pain.End.of.trial..actual.value.SD")
  changeFbaseline<-dataset_Pain_Endoftrial%>%
    select(matches("^Pain\\.End.of.trial\\.change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Pain.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$Pain.End.of.trial.change.from.baseline.Mean, .after = "Pain.End.of.trial.change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Pain.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$Pain.End.of.trial.change.from.baseline.SD, .after = "Pain.End.of.trial.change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Pain.End.of.trial..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Pain.End.of.trial..actual.value.Mean_edit"]<-sub_actuals[index2,"Pain.End.of.trial..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Pain.End.of.trial..actual.value.SD_edit" = actuals$Pain.End.of.trial..actual.value.SD, .after = "Pain.End.of.trial..actual.value.SD")
  index3<-which(is.na(actuals$Pain.End.of.trial..actual.value.SD))
  
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Pain.End.of.trial..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Pain.End.of.trial..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Pain.End.of.trial..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.End.of.trial..actual.value.Upper.CI))
      
      if(length(index7)!=0){
        temp_sd[index7,"Pain.End.of.trial..actual.value.SD_edit"]<-((temp_sd[index7,"Pain.End.of.trial..actual.value.Upper.CI"]-temp_sd[index7,"Pain.End.of.trial..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.End.of.trial..actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Pain.End.of.trial..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.End.of.trial..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.End.of.trial..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.End.of.trial..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.End.of.trial..actual.value.N"][which.max(temp_sd2[index9,"Pain.End.of.trial..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.End.of.trial..actual.value.N"][which.min(temp_sd2[index9,"Pain.End.of.trial..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.End.of.trial..actual.value.p.value"]))
      temp_sd2[index9,"Pain.End.of.trial..actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  options(warn = defaultW)
  
  
  ####End of trial Change from baseline
  changeFbaseline<-dataset_Pain_Endoftrial%>%
    select(matches("^Pain\\.End.of.trial\\.change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  #colnames(changeFbaseline)
  defaultW <- getOption("warn") 
  options(warn = -1) 
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Pain.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$Pain.End.of.trial.change.from.baseline.Mean, .after = "Pain.End.of.trial.change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Pain.End.of.trial.change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Pain.End.of.trial.change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Pain.End.of.trial.change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Pain.End.of.trial.change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Pain.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$Pain.End.of.trial.change.from.baseline.SD, .after = "Pain.End.of.trial.change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Pain.End.of.trial.change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Pain.End.of.trial.change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Pain.End.of.trial.change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Pain.End.of.trial.change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Pain.End.of.trial.change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Pain.End.of.trial.change.from.baseline.SD_edit"]<-((temp_sd[index7,"Pain.End.of.trial.change.from.baseline.Upper.CI"]-temp_sd[index7,"Pain.End.of.trial.change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Pain.End.of.trial.change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Pain.End.of.trial.change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Pain.End.of.trial.change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.N"][which.max(temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.N"][which.min(temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.p.value"]))
      temp_sd2[index9,"Pain.End.of.trial.change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}

options(warn = defaultW)
##Replacing processed Pain outcome-EndOfTrial
Pain_End_of_trial_continuous_outcome<-cbind(actuals,changeFbaseline)


####Final Continuous outcome dataset-Pain
final_mtx_irDataset_actual_Pain<-cbind(final_mtx_irDataset_actual[,c(1:27,30:37,40:43,47:57)],
                                  Pain_12months_continuous_outcome,
                                  Pain_3months_continuous_outcome,
                                  Pain_6months_continuous_outcome,
                                  PainBaseline_continuous_outcome,
                                  Pain_End_of_trial_continuous_outcome,
                                  DoseLevel01=final_mtx_irDataset_actual[,"DoseLevel01"])
final_mtx_irDataset_actual_Pain[,c(51:141)] <- lapply(final_mtx_irDataset_actual_Pain[,c(51:141)], function(x) as.numeric(as.character(x)))#converting character numeric to numeric

View(final_mtx_irDataset_actual_Pain)

write.csv(final_mtx_irDataset_actual_Pain,"Pain MTX-IR AUG2026.csv",row.names = F)

##fill time points of no interest with NA
final_mtx_irDataset_actual_Pain2 <- final_mtx_irDataset_actual_Pain

# End of trial columns
eot_cols <- grep(
  "^Pain\\.End\\.of\\.trial",
  names(final_mtx_irDataset_actual_Pain2 ),
  value = TRUE
)

# 12-month columns
m12_cols <- grep(
  "^Pain\\.12\\.months",
  names(final_mtx_irDataset_actual_Pain2 ),
  value = TRUE
)

# 6-month columns
m6_cols <- grep(
  "^Pain\\.6\\.months",
  names(final_mtx_irDataset_actual_Pain2 ),
  value = TRUE
)

# 3-month columns
m3_cols <- grep(
  "^Pain\\.3\\.months",
  names(final_mtx_irDataset_actual_Pain2 ),
  value = TRUE
)

# Apply same rules as dichotomous outcome
final_mtx_irDataset_actual_Pain2[
  final_mtx_irDataset_actual_Pain2$`Trial.duration..in.months` %in%
    c(3,4,6,9,11,12,13,18),
  eot_cols
] <- NA

final_mtx_irDataset_actual_Pain2[
  final_mtx_irDataset_actual_Pain2$`Trial.duration..in.months` %in%
    c(3,4,6,23,24,63),
  m12_cols
] <- NA

final_mtx_irDataset_actual_Pain2[
  final_mtx_irDataset_actual_Pain2$`Trial.duration..in.months` %in%
    c(3,4,9,11,12,13,18,23,24,63),
  m6_cols
] <- NA

final_mtx_irDataset_actual_Pain2[
  final_mtx_irDataset_actual_Pain2$`Trial.duration..in.months` %in%
    c(6,9,11,12,13,18,23,24,63),
  m3_cols
] <- NA

tttt<-final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study_no==97,c(51:140)]
tttt_studyInfo<-final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study_no==97,c(1:3)]
tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
tttt[colSums(!is.na(tttt)) > 0]
as.vector(which(colSums(!is.na(tttt)) > 0))

#Tanaka 2011, no baseline values
months3_change_cols<-grep(
  "^Pain\\.3\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Tanaka 2011",c(months3_change_cols)]<-NA

#Emery 2010, no baseline values
months6_change_cols<-grep(
  "^Pain\\.6\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Emery 2010",c(months6_change_cols)]<-NA


#Tanaka 2019.a, no baseline values
months6_change_cols<-grep(
  "^Pain\\.6\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Tanaka 2019.a",c(months6_change_cols)]<-NA


#Cohen 2020, no baseline values
months3_change_cols<-grep(
  "^Pain\\.3\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Cohen 2020",c(months3_change_cols)]<-NA

#Keystone 2015
months3_change_cols<-grep(
  "^Pain\\.3\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Keystone 2015",c(months3_change_cols)]<-NA

#Kim 2000
final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Kim 2000",c("Pain.3.months...change.from.baseline.SD_edit","Pain.Baseline.SD_edit")]<- c((100 - 11) /
                                                                                                                                                                                (2 * qnorm((20 - 0.375) / (20 + 0.25))),(81 - 10) /
                                                                                                                                                                                (2 * qnorm((20 - 0.375) / (20 + 0.25))))

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Kim 2000",c("Pain.3.months...change.from.baseline.Mean_edit")]<- final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Kim 2000",c("Pain.3.months...change.from.baseline.Median")]

#Karanikolas 2006, no baseline values
months12_change_cols<-grep(
  "^Pain\\.12\\.months\\..change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Karanikolas 2006",c(months12_change_cols)]<-NA


#Dougados 2013, no baseline values
months6_change_cols<-grep(
  "^Pain\\.6\\.months\\...change",
  names(final_mtx_irDataset_actual_Pain2),
  value = TRUE
)

final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study.Identifier=="Dougados 2013",c(months6_change_cols)]<-NA



temp_data4<-list()
empty_entries4<-list()


for(i in 1:length(unique(final_mtx_irDataset_actual_Pain2$Study_no))){
  tttt<-final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study_no==i,c(51:140)]
  tttt_studyInfo<-final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study_no==i,c(1:3)]
  tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
  
  if(sum(which(colSums(!is.na(tttt)) > 0))==0){
    empty_entries4[[i]]<-final_mtx_irDataset_actual_Pain2[final_mtx_irDataset_actual_Pain2$Study_no==i,c(1:3,52:140)]
  }else if(identical(c(85:89),as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##Only baseline value
    eee<-cbind(tttt_studyInfo,tttt[colSums(!is.na(tttt)) > 0]%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(tttt[colSums(!is.na(tttt)) > 0])[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data4[[i]]<- eee
  }else if(any(c(5:9,25:29,45:49,75:79)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##actual only but need to choose month
    getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
    getting_actuals2<-getting_actuals[, grepl("actual.value",names(getting_actuals))]
    if(ncol(getting_actuals2)==5){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[3])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==6){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[4])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==10){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==11){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==15){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[13])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==16){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[15])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==18){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[17])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }else if(ncol(getting_actuals2)==20){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[18])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data4[[i]]<- eee
    }
  }else if(any(c(11,31,51,71)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##Baseline value and change in baseline value
    n_col<-ncol(tttt[colSums(!is.na(tttt)) > 0])
    #baseline<-tttt[colSums(!is.na(tttt)) > 0][,sort(n_col:((n_col-5)+1))]
    baseline<-tttt[colSums(!is.na(tttt)) > 0][,grep("^Pain\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE)]
    #change_mean_values<-tttt[colSums(!is.na(tttt)) > 0][,-sort(n_col:((n_col-5)+1))]
    change_mean_values<-tttt[colSums(!is.na(tttt)) > 0]%>%
      select(!grep("^Pain\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE))
    change_mean_values2<-as.data.frame(change_mean_values[ , grepl( "change.from.baseline.Mean_edit" , names(change_mean_values) ) ])
    actual_mean<-baseline[,2]+change_mean_values2[,ncol(change_mean_values2)]
    baseline[,"Pain.Baseline.Mean_edit"]<-actual_mean
    replaced_baseline_mean_with_actual<-baseline
    eee<-cbind(tttt_studyInfo,replaced_baseline_mean_with_actual%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(replaced_baseline_mean_with_actual)[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data4[[i]]<- eee
  }
  
  
}

View(do.call(rbind,temp_data4))
View(do.call(rbind,empty_entries4))
Pain_dataset_analysis_mtxIR<-do.call(rbind,temp_data4)
Pain_dataset_analysis_mtxIR1<-Pain_dataset_analysis_mtxIR
row.names(Pain_dataset_analysis_mtxIR1)<-NULL

#####Formatting the month column
Pain_dataset_analysis_mtxIR1$Timepoint <- gsub(
  "Pain\\.|\\.months\\.{2,3}actual\\.value\\.(N|Mean_edit|Median|Lower|SD|SD_edit|_edit\\.CI)|\\.(Mean_edit|Median)",
  "",
  Pain_dataset_analysis_mtxIR1$Timepoint
)

Pain_dataset_analysis_mtxIR1<- Pain_dataset_analysis_mtxIR1[order(Pain_dataset_analysis_mtxIR1$Study_no),]

row.names(Pain_dataset_analysis_mtxIR1)<-NULL

Pain_dataset_analysis_mtxIR1<- Pain_dataset_analysis_mtxIR1 %>%
  dplyr::rename(treatment = Intervention)

##Drop furst study due to high ROB
Pain_dataset_analysis_mtxIR1<- Pain_dataset_analysis_mtxIR1 %>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012"))

library(dplyr)
for(i in 1:nrow(Pain_dataset_analysis_mtxIR1)){
  ppp<- master.medtable%>% filter(Pain_dataset_analysis_mtxIR1$treatment[i] ==master.medtable$medname )
  Pain_dataset_analysis_mtxIR1$treatment[i]<-ppp[1,1]
}

##Creating initial legend
PAIN.mtxIR.medlegend <- new.medtable(Pain_dataset_analysis_mtxIR1,'treatment')
PAIN.mtxIR.medlegend

##NMA plot reason for narrowing of network
nma.plot(Pain_dataset_analysis_mtxIR1[,1],Pain_dataset_analysis_mtxIR1[,3],data = Pain_dataset_analysis_mtxIR1, trtname = as.vector(PAIN.mtxIR.medlegend$medname))

################################################################################################################
#trying the narrowing a network function
#############################################################################################
Pain_dataset_analysis_mtxIR1
Pain_dataset_analysis_mtxIR1$Study_no<-renumber(Pain_dataset_analysis_mtxIR1,"Study_no")
Pain_dataset_analysis_mtxIR1_old<-Pain_dataset_analysis_mtxIR1

##Converting those with 1-10 scale to 1-100 VAS
Pain_dataset_analysis_mtxIR1$SD <- ifelse(Pain_dataset_analysis_mtxIR1$Mean > 10, Pain_dataset_analysis_mtxIR1$SD, Pain_dataset_analysis_mtxIR1$SD * 10)
Pain_dataset_analysis_mtxIR1$Mean <- ifelse(Pain_dataset_analysis_mtxIR1$Mean > 10, Pain_dataset_analysis_mtxIR1$Mean, Pain_dataset_analysis_mtxIR1$Mean * 10)

##force mtx+ada arm to be 9.6 in the Breedveld 2006 study
#Pain_dataset_analysis_mtxIR1[Pain_dataset_analysis_mtxIR1$Study.Identifier=="Breedveld 2006",c("SD","Mean")][3,]<-Pain_dataset_analysis_mtxIR1[Pain_dataset_analysis_mtxIR1$Study.Identifier=="Breedveld 2006",c("SD","Mean")][3,]/10

Pain_dataset_analysis_mtxIR1 <- Pain_dataset_analysis_mtxIR1 %>%
  dplyr::rename("study" = Study_no,"study_name" = Study.Identifier)

#Narrow 1
mtxIR_pain_narrow<-narrow.all(Pain_dataset_analysis_mtxIR1)
mtxIR_pain_narrow


#sub network plot 1
mtxIR.pain.d<-mtxIR_pain_narrow[[1]]
length(unique(mtxIR.pain.d$treatment))
mtxIR_pain_subl<-new.medtable(mtxIR.pain.d,'treatment')
mtxIR.pain.d_download<-mtxIR.pain.d
mtxIR.pain.d_download$treatment<-renumber(mtxIR.pain.d_download,"treatment")
mtxIR.pain.d_download$treatment<-selected_treatment_name(mtxIR.pain.d_download$treatment,mtxIR_pain_subl)
write.csv(mtxIR.pain.d_download[,-c(8,9)],"Pain MTX-IR dataset_AUG2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("Pain MTX-IR Network plot with named interventions.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.pain.d[,1],mtxIR.pain.d[,3],data = mtxIR.pain.d, trtname = as.vector(mtxIR_pain_subl$medname),title = "PAIN MTX-IR Network plot")
# Close the pdf file
dev.off() 
mtxIR.pain.d1<-mtxIR.pain.d
mtxIR.pain.d1$treatment<-renumber(mtxIR.pain.d1,"treatment")
# 2. Open jpeg file
jpeg("Pain MTX-IR Network plot without named interventions.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.pain.d1[,1],mtxIR.pain.d1[,3],data = mtxIR.pain.d1,title = "PAIN MTX-IR Network plot")
# Close the pdf file
dev.off() 


##############################################################################
#Performing sub network plot 3
library(coda)
library(gemtc)
library(netmeta)
library(tidyr)
View(mtxIR.pain.d)
mtxIR_pain_subl

network_data_pain_mtxIR_1<-mtxIR.pain.d%>% 
  select(study, study_name, treatment, SD, sampleSize, Mean ,Timepoint)
network_data_pain_mtxIR_1$treatment<-renumber(network_data_pain_mtxIR_1,"treatment")
network_data_pain_mtxIR_1$study<-renumber(network_data_pain_mtxIR_1,"study")
network_data_pain_mtxIR_1<-network_data_pain_mtxIR_1[order(network_data_pain_mtxIR_1$study),]
row.names(network_data_pain_mtxIR_1)<-NULL
nma.plot(network_data_pain_mtxIR_1$study,network_data_pain_mtxIR_1$treatment,data = network_data_pain_mtxIR_1, trtname = as.vector(mtxIR_pain_subl$medname),title = "PAIN MTX-IR Network plot")

network_data_pain_mtxIR_1 <- network_data_pain_mtxIR_1 %>%
  dplyr::rename("std.dev" = SD, "mean"= Mean)

p_pain_MD<-pairwise(treat=treatment,mean=mean,sd=std.dev,data=network_data_pain_mtxIR_1,
                    studlab = study,n=sampleSize,sm="MD")
#colnames(p_pain_MD)[c(1,4:5)]<-c("study","diff", "std.err")
p_pain_MD <- p_pain_MD %>%
  dplyr::rename("diff" = TE, "std.err" = seTE)

p_pain_MD %>%
  dplyr::select(c(study_name,diff, std.err,treat1,treat2)) %>%
  pivot_longer(-study_name,
               names_to = c(".value"),
               names_pattern = "(..)") %>% 
  set_colnames(c("study", "diff", 
                 "std.err", "treatment")) -> data2

data2 <- clean_and_impute_se(data2)

network_data_pain_mtxIR_1<-mtc.network(data.re=data2)
summary(network_data_pain_mtxIR_1)
View(network_data_pain_mtxIR_1$data.re)

##RUN the model
#network_data_DAS28_bio_1.model  <- mtc.model(network_data_DAS28_bio_1 , likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
set.seed(5811)
network_data_pain_mtxIR_1.model <- mtc.model(network_data_pain_mtxIR_1, likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
#network_data_HAQ_bio_1.model <- mtc.model(network_data_HAQ_bio_1, likelihood='normal', link='identity')

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_pain_mtxIR_1.model$inits <- mapply(c, network_data_pain_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_pain<-mtc.run(network_data_pain_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_pain)

View(as.data.frame(relative.effect.table(results_mtxIR_pain)))
View(as.data.frame(relative.effect.table.2DP2(results_mtxIR_pain, comparator = "1")))

refME_mtxIR_pain<-as.data.frame(relative.effect.table.2DP2(results_mtxIR_pain, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refME_mtxIR_pain)<-mapvalues(row.names(refME_mtxIR_pain),
                                       from = c(1:(nrow(mtxIR_pain_subl))),
                                       to = as.vector(mtxIR_pain_subl$medname))
colnames(refME_mtxIR_pain)<-mapvalues(colnames(refME_mtxIR_pain),
                                      from = c(1:(nrow(mtxIR_pain_subl))),
                                      to = as.vector(mtxIR_pain_subl$medname))
View(refME_mtxIR_pain)

index_diagonal5<-data.frame(index1=1:nrow(refME_mtxIR_pain),index2=1:nrow(refME_mtxIR_pain))
options(warn=-1)
for(i in 1:nrow(refME_mtxIR_pain)){
  refME_mtxIR_pain[index_diagonal5$index1[i],index_diagonal5$index2[i]]<-selected_treatment_name(as.numeric(refME_mtxIR_pain[row(refME_mtxIR_pain)==col(refME_mtxIR_pain)])[i],mtxIR_pain_subl)
}
options(warn=0)
View(refME_mtxIR_pain)


write.csv(refME_mtxIR_pain,"Pain MTX-IR pairwise comparisons_AUG2026.csv", row.names = T)
#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "PAIN")
writeData(wb, "PAIN", refME_mtxIR_pain, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)


##GRADE Assessment
nma_estimate_mtxIR_Pain_dataframe<-nma.estimate.table.MD(results_mtxIR_pain)#use results
View(nma_estimate_mtxIR_Pain_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_pain_mtxIR_1.ns<-mtc.nodesplit(network_data_pain_mtxIR_1, thin=50,likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 
#network_data_HAQ_mtxIR_1.ns<-mtc.nodesplit(network_data_HAQ_mtxIR_1, thin=50,likelihood='normal', link='identity') 
mtxIR_pain_subl
names(network_data_pain_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_pain_mtxIR_1.ns)

# Overall summary and plot
summary.network_data_pain_mtxIR_1.ns<- summary(network_data_pain_mtxIR_1.ns)

#Direct effects
direct.data<-summary.network_data_pain_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_pain_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_pain_mtxIR_1.ns$cons.effect

#P-values
P_values_nodesplit_painmtxIR<-summary.network_data_pain_mtxIR_1.ns$p.value
P_values_nodesplit_painmtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_painmtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_painmtxIR2<-P_values_nodesplit_painmtxIR
P_values_nodesplit_painmtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_painmtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_painmtxIR2$comb<-paste0(P_values_nodesplit_painmtxIR2$t1,":",P_values_nodesplit_painmtxIR2$t2)


direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(direct.data[,c(3,4,5)],2),round(indirect.data[,c(3,4,5)],2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                -round(direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)
##Acquiring both direct and indirect estimate from gemtc
grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_Pain_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)

network_data_pain_mtxIR_1<-mtxIR.pain.d[,c(1:7)]
network_data_pain_mtxIR_1$treatment<-renumber(network_data_pain_mtxIR_1,"treatment")
network_data_pain_mtxIR_1$study<-renumber(network_data_pain_mtxIR_1,"study")
network_data_pain_mtxIR_1<-network_data_pain_mtxIR_1[order(network_data_pain_mtxIR_1$study),]
row.names(network_data_pain_mtxIR_1)<-NULL

#Adding ROB Data
Pain_mtxIR<-unique(network_data_pain_mtxIR_1$study_name)
data_qq_mtxIR<-network_data_pain_mtxIR_1
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%Pain_mtxIR,]

col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<- current_assessment_for_studies_mtxIR %>%
  select(c(col_names))


list_appraisal_mtxIR9<-list()
for(i in 1:length(Pain_mtxIR)){
  temp11<-data_qq_mtxIR[data_qq_mtxIR$study_name==Pain_mtxIR[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==Pain_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal_mtxIR9[[i]]<-temp13
}
Pain_mtxIR_withROB<-do.call(rbind,list_appraisal_mtxIR9)
View(Pain_mtxIR_withROB)
Pain_mtxIR_withROB2<-Pain_mtxIR_withROB%>%
  select(c("study","study_name", "treatment", "SD", "sampleSize",
           "Mean","Timepoint","OVERALL..ALL.OTHER.OUTCOMES"))
View(Pain_mtxIR_withROB2)

###Adding Indirectness data
study_unique<-unique(Pain_mtxIR_withROB2$study_name)
data_nn<-Pain_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp[!is.na(resolved_indirectness_consensus_temp$`Final Rating`),]
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]
row.names(current_indirect_assessment_for_studies)<-NULL
# current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:2),],current_indirect_assessment_for_studies[2,],current_indirect_assessment_for_studies[c(3,4),])
# current_indirect_assessment_for_studies$intervention[3]<-"Baricitinib (4 mg)"
# row.names(current_indirect_assessment_for_studies)<-NULL
list_appraisal10<-list()
for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal10[[i]]<-temp13
}
Pain_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal10)
Pain_mtxIR_withROB2_indrectness <- Pain_mtxIR_withROB2_indrectness%>%
  dplyr:: rename(ROB = OVERALL..ALL.OTHER.OUTCOMES , "Indirectness" = `Final Rating`)

Pain_mtxIR_withROB2_indrectness$Indirectness<-tolower(Pain_mtxIR_withROB2_indrectness$Indirectness)
#Pain_mtxIR_withROB2_indrectness[Pain_mtxIR_withROB2_indrectness$study_name=="Kremer 2016","Indirectness"]<-c("low","low")
View(Pain_mtxIR_withROB2_indrectness)

Pain_mtxIR_data_new_withROB2_indrectness2<-Pain_mtxIR_withROB2_indrectness
Pain_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(Pain_mtxIR_data_new_withROB2_indrectness2,"treatment")
Pain_mtxIR_data_new_withROB2_indrectness2$ROB<-ifelse(Pain_mtxIR_data_new_withROB2_indrectness2$ROB=="unclear","some concerns",Pain_mtxIR_data_new_withROB2_indrectness2$ROB)
# Pain_mtxIR_data_new_withROB2_indrectness2_details<-getContributionMatrix(Pain_mtxIR_data_new_withROB2_indrectness2,type = "long_binary",model="random",tau=NA,  sm="OR")
p_pain_MD1<-pairwise(treat=treatment,mean=Mean,sd=SD,data=network_data_pain_mtxIR_1,
                     studlab = study,n=sampleSize,sm="MD")
p_pain_MD2<-p_pain_MD1
colnames(p_pain_MD2)[c(1,4:5)]<-c("studlab","TE","seTE")
net1 <- netmeta(p_pain_MD2)
dat.net<-netsplit(net1)
dat.net.contrib<-netcontrib(net1, method = "r")

direct_evidence_new<-cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)])
View(direct_evidence_new)
#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")

# grade_table1<-cbind(nma_estimate_mtxIR_Pain_dataframe[,c(1,2)],
#                     as.data.frame(matrix(NA,nrow(nma_estimate_mtxIR_Pain_dataframe),6)),
#                     nma_estimate_mtxIR_Pain_dataframe[,-c(1,2)])
# colnames(grade_table1)[3:8]<-c("MD1","low1","high1","MD2","low2","high2")

#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]


##Getting direct evidence for intervention in only multi-arm trials
zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
colnames(zzz)<-c("treat1", "treat2")
grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]


grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_pain_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_pain_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))


#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)

graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
Pain_output_df2_list<-list()
I_square_list<-list()

for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-Pain_mtxIR_data_new_withROB2_indrectness2[Pain_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select(study, Mean.1, SD.1, Mean.2, SD.2 )
  temp11.1[,c("SD.1","SD.2")]<-temp11.1[,c("SD.1","SD.2")]/(sqrt(temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c("sampleSize.1" , "sampleSize.2")]))
  
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","se.t","x.c","se.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    #model.file<-"Cont_direct_SMD_onestudy.txt"
    model.file<-"Cont_direct_onestudy.txt"
    n.chains=3
    set.seed(66451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c#,
                                   # "pooled.sd.c"=sqrt(sum(temp11[c(4,11)]^2)),
                                   # "pooled.sd.t"=sqrt(sum(temp11[c(4,11)]^2))
                         ),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("diff")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(MD=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    model.file<-"Cont_direct_multiplestudy_RE.txt" #Check this again  #"Cont_direct_SMD_multiplestudy_RE.txt"
    n.chains=3
    set.seed(77451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c, "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("diff")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_dataPain_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select(study, sampleSize.1, Mean.1, SD.1, sampleSize.2, Mean.2, SD.2 )
  #nma_dataDAS28_3[,c("SD.1","SD.2")]<-nma_dataDAS28_3[,c("SD.1","SD.2")]/(sqrt(nma_dataDAS28_3[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(2,5)]))
  nma_dataPain_4<-nma_dataPain_3
  colnames(nma_dataPain_4)<-c("NS","n.t","mean.t","sd.t","n.c","mean.c","sd.c")
  row.names(nma_dataPain_4)<-NULL
  m.cont<-metacont(n.e = n.t, 
                   mean.e = mean.t, 
                   sd.e = sd.t, 
                   n.c = n.c, 
                   mean.c = mean.c, 
                   sd.c =sd.c,
                   data = nma_dataPain_4,
                   sm = "SMD",
                   method.tau = "PM",
                   random = TRUE,
                   method.random.ci = TRUE,
                   title = "Pain (final)"
                   
  )
  
  summary_m.cont<-summary(m.cont)
  contribution_weight<-as.character(round((summary_m.cont$w.random/sum(summary_m.cont$w.random))*100,2))
  pooled_estimate2<-data.frame(Mean=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.cont$tau2,4),
                          Lower=round(summary_m.cont$lower.tau2,4),
                          Upper=round(summary_m.cont$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.cont$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.cont$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.cont$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.cont$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  Mean_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1,temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(5,12)])
  base_data <- tibble(mean  = Mean_individual_studies$MD, 
                      lower = Mean_individual_studies$Lower,
                      upper = Mean_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      treatment1 = paste0(as.character(temp11.2$x.c),"      ","(",as.character(temp11.2$sampleSize.2),")"),
                      treatment2 = paste0(as.character(temp11.2$x.t),"      ","(",as.character(temp11.2$sampleSize.1),")"),
                      MD = paste(format(round(Mean_individual_studies$MD,2),nsmall=2),paste("(",format(round(Mean_individual_studies$Lower,2),nsmall=2),", ",format(round(Mean_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$Mean, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   MD = paste(format(round(pooled_estimate2$Mean,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("Pain (final)", "Study"),
                   treatment1 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.2),mtxIR_pain_subl)),
                   treatment2 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.1),mtxIR_pain_subl)),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   MD = c("", "MD"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  Pain_output_df <- bind_rows(header,
                              base_data,
                              empty_row,
                              pooled)
  
  
  
  
  Pain_output_df2<-rbind(Pain_output_df,Pain_output_df[c(nrow(Pain_output_df)),])
  Pain_output_df2[nrow(Pain_output_df),1]<-I2
  Pain_output_df2[nrow(Pain_output_df),2]<-NA
  Pain_output_df2[nrow(Pain_output_df),c(6:11)]<-c(NA,NA,NA,NA,NA,NA,NA)
  Pain_output_df2_list[[i]]<-Pain_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-Pain_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, MD,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(-1, 1), 
               xlog = F, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_pain_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_pain_subl))
    )
  
  
  
}#END

n_studies<-c()
n_participants<-c()

extract_num<-function(data){
  return(str_extract(data,"(?<=\\()\\d+(?=\\))"))
}
countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}




for(i in 1:length(graph_list)){
  df<-Pain_output_df2_list[[i]][-c(1,2,countdown(length(Pain_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df,1,extract_num))))
}

graph_list[[1]]

Grade_assessment_dataset_direct_Pain<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_pain_subl),
                                                 selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_pain_subl),
                                                 Outcome="PAIN",
                                                 `Effect measure`="MD",
                                                 `Number of studies`=n_studies,
                                                 `Number of participant`=n_participants,
                                                 ROB="",Indirectness="",Imprecision="",
                                                 Inconsistency="",
                                                 Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_Pain)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                        "Intevention",
                                                                        "Effect measure",
                                                                        "Number of studies",
                                                                        "Number of participant",
                                                                        "Risk of Bias",
                                                                        "Publication Bias",
                                                                        "Certainty of evidence (ignores imprecision)")


View(Grade_assessment_dataset_direct_Pain)

#Saving grade assessment details for Pain mtxIRLOGICS as an object of lists
grade_details_Pain_mtxIR_new<-list(graph_list,Pain_output_df2_list,Grade_assessment_dataset_direct_Pain,I_square_list)
saveRDS(grade_details_Pain_mtxIR_new,file="List of Pain MTX-IR data components for GRADE_AUG2026.rds")

pain_output_one_study_list <- list()
pain_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

n_studies<-c()
for(i in 1:length(graph_list)){
  df<-Pain_output_df2_list[[i]][-c(1,2,countdown(length(Pain_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
}

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
pain_output_one_study_list  <- Pain_output_df2_list[one_idx]
pain_output_more_than_2_list <- Pain_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_Pain[
    Grade_assessment_dataset_direct_Pain$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_Pain[
    Grade_assessment_dataset_direct_Pain$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_pain_mtxIR_one_study <- list(
  graph_list_one_study_list,
  pain_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_pain_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  pain_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_pain_mtxIR_one_study,
  file = "List of Pain MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_pain_mtxIR_more_than_2,
  file = "List of Pain MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)


#####################################
#####################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_pain_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")

get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_Pain$Treatment,
                                                Grade_assessment_dataset_direct_Pain$Intevention),
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(Pain_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","MD"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 
glen_direct_evidence_rating_pain<-read.csv("Glen_2026-09-04 16_20_PAIN-MTXIR_complete.csv")
glen_direct_evidence_rating_pain2<-glen_direct_evidence_rating_pain[,-c(13:15)]
colnames(glen_direct_evidence_rating_pain2)[1:2]<-colnames(grade_table4)[1:2]


glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_pain2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"PAIN"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"MD"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR_pain_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR_pain_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}


##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for Pain_MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for Pain_MTX-IR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR Pain OUTCOME
#initial_gradeTable_haq<-read.csv("Initial GRADE table for HAQ modified.csv")
#initial_gradeTable_haq<-read.csv("Initial GRADE table for HAQ_OCT31 modified.csv")
initial_gradeTable_pain<-Temp_GRADE5
estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))



initial_gradeTable_pain2<-initial_gradeTable_pain[match(all_nma_estimates,initial_gradeTable_pain$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()

#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_pain2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_pain2.1<-cbind(initial_gradeTable_pain2[,c(1:14)],first_order_loop_contri,
                                  initial_gradeTable_pain2[,-c(1:14)])

initial_gradeTable_pain2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_pain2.1,
                                                             all_nma_estimates,
                                                             estimate_contribution_matrix,
                                                             glen_direct_evidence_rating6_onlyDirect,
                                                             top_n = 2))

# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_pain2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_pain2.1<-result_df
initial_gradeTable_pain2.1$combI<-str_replace(initial_gradeTable_pain2.1$comb, ":", " vs ")

initial_gradeTable_pain2.1<-process_dataset(initial_gradeTable_pain2.1, 2)
initial_gradeTable_pain3<-initial_gradeTable_pain2.1[match(initial_gradeTable_pain$comb,initial_gradeTable_pain2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for PAIN_MTX-IR outcome_2026.csv",row.names = T)
#
###########Assessing Intransitivity
sub_data_intransitivity <- initial_gradeTable_pain3[initial_gradeTable_pain3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_pain3[initial_gradeTable_pain3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_pain3.1<-cbind(initial_gradeTable_pain3[,c(1:23)],
                                  Intransitivity="",
                                  `Certainty of evidence (ignores imprecision)-Indirect`="",
                                  initial_gradeTable_pain3[,-c(1:23)])
initial_gradeTable_pain3.1[initial_gradeTable_pain3.1$Indirect!="-",c("Intransitivity",
                                                                      "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                            "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_pain4<-initial_gradeTable_pain3.1


#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)
grade_temp<-cbind(initial_gradeTable_pain4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_pain4[,-c(1:25)])

###Applying the P-values
P_values_nodesplit_painmtxIR2
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_painmtxIR3 <- P_values_nodesplit_painmtxIR2[match(drt$comb, P_values_nodesplit_painmtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_painmtxIR3$p,2)


subPAINSOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatPAIN<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subPAINSOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(1)
IndirPer <- numeric(1)

# # Loop to calculate dirPer and IndirPer
# for (i in 1:nrow(subcontriMatPAIN)) {
#   dirPer[i] <- subcontriMatPAIN[i, colnames(subcontriMatPAIN) == row.names(subcontriMatPAIN)[i]]
#   IndirPer[i] <-  1- dirPer[i]
# }

#Because it has only one p-value
dirPer <- subcontriMatPAIN[names(subcontriMatPAIN)==subPAINSOF$comb]
IndirPer <-  1- dirPer
# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")


grade_temp2[match(subPAINSOF$comb, grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(subPAINSOF$comb, grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)

##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))

#When both direct and indirect are present
fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg


####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho

initial_gradeTable_pain5<-grade_temp2.2
initial_gradeTable_pain5<-initial_gradeTable_pain5[,-c(30,33)]
colnames(initial_gradeTable_pain5)[c(29,30)]<-c("Initial NMA evidence rating",
                                                "P-Value")

###After incoherence evaluation
initial_gradeTable_pain5.1<-cbind(initial_gradeTable_pain5[,c(1:31)],
                                  qoeAfterIncohenrence="",
                                  initial_gradeTable_pain5[,-c(1:31)])
initial_gradeTable_pain5.1[(initial_gradeTable_pain5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_pain5.1[(initial_gradeTable_pain5.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_pain5.1[initial_gradeTable_pain5.1$Incoherence!="",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_pain5.1[(initial_gradeTable_pain5.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence

#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_pain5.1[(initial_gradeTable_pain5.1$Direct!="-")&(initial_gradeTable_pain5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_pain5.1[(initial_gradeTable_pain5.1$Direct!="-")&(initial_gradeTable_pain5.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_pain6<-cbind(initial_gradeTable_pain5.1,
                                p.estimate=apply(as.data.frame(initial_gradeTable_pain5.1$NMA),1,function(x){
                                  as.numeric(sub(" .*", "", x))
                                }),
                                t(apply(as.data.frame(initial_gradeTable_pain5.1$NMA),1,function(x){
                                  as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                  
                                })))

colnames(initial_gradeTable_pain6)[c(36,37)]<-c("l.bound","u.bound")

to_higher_lower_text <- function(e1) {
  sign_val <- sign(e1)
  abs_val <- abs(round(e1, 2))
  text <- ifelse(sign_val < 0, "lower", "higher")
  return(paste(abs_val, text))
}

####Running Continuous NMA model
# load data for MTC
##Getting baseline risk using the placebo intervention

Pain_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(Pain_mtxIR_data_new_withROB2_indrectness2,"treatment")
Pain_mtxIR_data_new_withROB2_indrectness2$study<-renumber(Pain_mtxIR_data_new_withROB2_indrectness2,"study")
Pain_data<-Pain_mtxIR_data_new_withROB2_indrectness2
Pain_data<-Pain_data[order(Pain_data$study),c(1:6)]


Pain_data_wide2<-reshape(number_of_arms(Pain_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(Pain_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(Pain_data_wide2[,c(3,7,11)],1,xcv)

Pain_data_wide3<-cbind(Pain_data_wide2,na=numberOfArms)
View(Pain_data_wide3)
#finding standard errors
Pain_data_wide3$SD.1<-Pain_data_wide3$SD.1/sqrt(Pain_data_wide3$sampleSize.1)
Pain_data_wide3$SD.2<-Pain_data_wide3$SD.2/sqrt(Pain_data_wide3$sampleSize.2)
Pain_data_wide3$SD.3<-Pain_data_wide3$SD.3/sqrt(Pain_data_wide3$sampleSize.3)
colnames(Pain_data_wide3)<-c("study", "study_name", "t.1","se.1","n.1", "x.1",
                             "t.2","se.2","n.2","x.2","t.3","se.3","n.3","x.3", "na")
                             #"t.4","se.4","n.4","x.4","na" )
MTCData_cont <- as.data.frame(Pain_data_wide3)
x <- data.matrix(MTCData_cont[,c(6,10,14)])
n <- data.matrix(MTCData_cont[,c(5,9,13)])
t <- data.matrix(apply(MTCData_cont[,c("t.1","t.2","t.3")],c(1,2),as.numeric))
se <- data.matrix(MTCData_cont[,c(4,8,12)])
na <- data.matrix(MTCData_cont[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- x[t[,'t.1']==1,"x.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
sePBO <- se[t[,'t.1']==1,"se.1"]

# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
#model.file4<-"Cont_multiarm_RE.txt"
model.file4<-"Cont_new_multiarm.txt"
n.chains=3
model.inits <- initv1 

cont.model3<-jags.model(file=model.file4,
                        data= list("x"=x,"t"=t,"se"=se,"na" = na[,1], "NT" = NT, "NS" = NS,
                                   "NTpbo"= NTpbo,"xPBO"=xPBO,"sePBO"=sePBO),
                        n.chains=n.chains,
                        inits = model.inits,
                        n.adapt=5000
)

update(object=cont.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("diff","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( cont.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]

#Converting estimates relative to baseline to absolutes
estimates_relative_to_baseline<-initial_gradeTable_pain6[1:(NT-1),c(1,2,35:37)]
typeof(estimates_relative_to_baseline)

absolutes_from_baseline<-cbind(baseline_risk =round(baseline_risk,2),round(baseline_risk,2)+estimates_relative_to_baseline[,c(3:5)])
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(round(baseline_risk,2),estimates_relative_to_baseline_combine[,7]))
colnames(all_baseline)<-"base"

com.list<-list()
unique_treat2<-unique(initial_gradeTable_pain6$treat2)
d.sub<-initial_gradeTable_pain6[initial_gradeTable_pain6$treat2==unique_treat2[1],c(1,2,35:37)]
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_pain6[initial_gradeTable_pain6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,base+d.sub[,c(3:5)])
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe<-com.dataframe[c(1:nrow(initial_gradeTable_pain6)),]
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9])

higher_lower_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"lower"
  }else{
    text<-"higher"
  }
  return(paste(abs(round(e1-e2,2)),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,13])))


prefix.assumedrisk <- ' Mean '
prefix.corrrisk<-' Mean '
prefix.md<- " MD "
suffix.assumedrisk <- ''
suffix.corrrisk <- ''
time<-1


result_com.dataframe_withtext2<-cbind(com.dataframe_withtext[,1:2],
                                      paste(prefix.assumedrisk,com.dataframe_withtext[,6], suffix.assumedrisk, sep=''),
                                      paste(prefix.corrrisk, com.dataframe_withtext[,7], suffix.corrrisk, ' (', com.dataframe_withtext[,8], ' to ',com.dataframe_withtext[,9],')', sep=''),
                                      paste(prefix.md, com.dataframe_withtext[,10], ' (', com.dataframe_withtext[,11], ' to ',com.dataframe_withtext[,12],')', sep=''),
                                      # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                      com.dataframe_withtext[,c(3:5,10:12)]
)

colnames(result_com.dataframe_withtext2)[c(3:11)]<-c("Assumed mean (comparator)",	"Corresponding mean (intervention)",
                                                     "Mean difference","p.estimate","l.bound","u.bound","","","")


initial_gradeTable_pain7<-cbind(initial_gradeTable_pain6[,-c(35:37)],
                                result_com.dataframe_withtext2[,-c(1,2)],pair="")
#colnames(initial_gradeTable_pain7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_pain8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_pain7$treat1,mtxIR_pain_subl)),
                                `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_pain7$treat2,mtxIR_pain_subl)),
                                initial_gradeTable_pain7)
colnames(initial_gradeTable_pain8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_pain8<-initial_gradeTable_pain8[,-c(46)]
write.csv(initial_gradeTable_pain8,"Pain MTX-IR SoF table without imprecision assessment AUG2026.csv",row.names = F)



effect.size.rating_list<-lapply(split(initial_gradeTable_pain8, seq_len(nrow(initial_gradeTable_pain8))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 1, large.effect = 3, framework = T)
})



as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)


###Example dataset rather than using imprecision application
initial_gradeTable_pain9<-initial_gradeTable_pain8
initial_gradeTable_pain9<-cbind(initial_gradeTable_pain9[,c(1:34,37:45,35:36)])
initial_gradeTable_pain9.1<-cbind(initial_gradeTable_pain9[,c(1:37)],
                                  Small.effect.threshold=1,
                                  Large.effect.threshold= 3,
                                  Framework = "NULL",
                                  Imprecision.1 = "",
                                  Certainty.of.evidence = "",
                                  initial_gradeTable_pain9[,c(44,45)])


initial_gradeTable_pain9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_Pain_mtxIR<-initial_gradeTable_pain9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_Pain_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_Pain_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_Pain_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_Pain_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_Pain_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_Pain_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_Pain_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_Pain_mtxIR[is.na(detailed_SOF_Pain_mtxIR$Certainty.of.evidence),])
detailed_SOF_Pain_mtxIR[is.na(detailed_SOF_Pain_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

detailed_SOF_Pain_mtxIR$comb<-initial_gradeTable_pain7$comb

write.csv(detailed_SOF_Pain_mtxIR,"Pain_MTX-IR_Detail_SoF_AUG2026.csv", row.names = F)

###############################################################################
###############################################################################
##############
##DAS28 Remission

final_mtx_irDataset_actual2 <- final_mtx_irDataset_actual%>%
  mutate(
    `DAS.28.Remission.End.of.trial.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
        NA,
        `DAS.28.Remission.End.of.trial.n`
      ),
    
    `DAS.28.Remission.End.of.trial.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,9,11,12,13,18),
        NA,
        `DAS.28.Remission.End.of.trial.N`
      ),
    
    `DAS.28.Remission.12.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,6,23,24,63),
        NA,
        `DAS.28.Remission.12.months.n`
      ),
    
    `DAS.28.Remission.12.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,23,24,63),
        NA,
        `DAS.28.Remission.12.months.N`
      ),
    
    `DAS.28.Remission.6.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24,63),
        NA,
        `DAS.28.Remission.6.months.n`
      ),
    
    `DAS.28.Remission.6.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(3,4,9,11,12,13,18,23,24,63),
        NA,
        `DAS.28.Remission.6.months.N`
      ),
    
    `DAS.28.Remission.3.months.n` =
      ifelse(
        `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24,63),
        NA,
        `DAS.28.Remission.3.months.n`
      ),
    `DAS.28.Remission.3.months.N` =
      ifelse(
        `Trial.duration..in.months` %in% c(6,9,11,12,13,18,23,24,63),
        NA,
        `DAS.28.Remission.3.months.N`
      )
  )


df_dasR <- final_mtx_irDataset_actual2%>% 
  select(matches("^DAS.28.Remission\\..*\\.(n|N)$")) %>%
  select(-contains("CRP"))%>%
  select(-contains(c("4.months","22.Weeks")))
colnames(df_dasR)

DAS28remission_dataset_mtxIR<-cbind(final_mtx_irDataset_actual2[,c(1:3)], df_dasR)

DAS28remission_dataset_mtxIR2<-as.data.frame(DAS28remission_dataset_mtxIR)

DAS28remission_datalist_mtxIR <-list()
n_study1<-unique(DAS28remission_dataset_mtxIR2$Study_no)

for(i in 1:length(n_study1) ){
  eee<-DAS28remission_dataset_mtxIR2[DAS28remission_dataset_mtxIR2$Study_no==n_study1[i],]
  eee[eee=="NA"]<-NA
  ggg<-eee[colSums(!is.na(eee)) > 0]
  if(ncol(ggg)==11){
    ggg<-ggg[,-c(4:9)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    DAS28remission_datalist_mtxIR[[i]]<-ggg
  } else if(ncol(ggg)==9){
    ggg<-ggg[,-c(4:7)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    DAS28remission_datalist_mtxIR[[i]]<-ggg
  }else if(ncol(ggg)==7){
    ggg<-ggg[,-c(4:5)]
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    DAS28remission_datalist_mtxIR[[i]]<-ggg
  }else if(ncol(ggg)==5){
    ggg<-cbind(ggg,rep(colnames(ggg)[4],nrow(ggg)))
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    DAS28remission_datalist_mtxIR[[i]]<-ggg
  }else if((ncol(ggg)<=4) || (ncol(ggg)==6) || (ncol(ggg)==8) || (ncol(ggg)==10)){
    ggg<-ggg[,c(1,2,3)]
    aaaa<-matrix(NA,nrow(ggg),3)
    ggg<-cbind(ggg,aaaa)
    colnames(ggg)<- c('study','study_name','treatment','responders', 'sampleSize','month')
    DAS28remission_datalist_mtxIR[[i]]<-ggg
  }
  
  
}


DAS28remission_datalist_mtxIR
DAS28remission_big_data_mtxIR = do.call(rbind, DAS28remission_datalist_mtxIR)
DAS28remission_big_data_mtxIR
DAS28remission_big_data_mtxIR<-cbind(DAS28remission_big_data_mtxIR,DoseLevel01=final_mtx_irDataset_actual2$DoseLevel01,TrialDuration=final_mtx_irDataset_actual2$Trial.duration..in.months)
View(DAS28remission_big_data_mtxIR)
dummy <- NULL
dummy[is.na(DAS28remission_big_data_mtxIR$month)] = 0
dummy[DAS28remission_big_data_mtxIR$month=="DAS.28.Remission.3.months.n"] = 1
#dummy[ACR50_big_data_biologic$month=="ACR.50.4.months.n"] = 2
dummy[DAS28remission_big_data_mtxIR$month=="DAS.28.Remission.6.months.n"] = 2
dummy[DAS28remission_big_data_mtxIR$month=="DAS.28.Remission.12.months.n"] = 3
dummy[DAS28remission_big_data_mtxIR$month=="DAS.28.Remission.End.of.trial.n"] = 4
DAS28remission_big_data_mtxIR<-cbind(DAS28remission_big_data_mtxIR,dummy)
colnames(DAS28remission_big_data_mtxIR)[9]<-"Time_point"
DAS28remission_big_data_mtxIR

##getting the actual contributed studies for DAS28remission
DAS28remission_big_data_mtxIR2<-DAS28remission_big_data_mtxIR[DAS28remission_big_data_mtxIR$Time_point!=0,]

##Solving error in Burnmester 2013 placebo arm
#DAS28remission_big_data_mtxIR2[DAS28remission_big_data_mtxIR2$study_name=="Burmester 2013",][,c(4,5)]<-c(c(2,8),c(120,119))
#DAS28remission_big_data_mtxIR2[DAS28remission_big_data_mtxIR2$study_name=="Burmester 2013",][,c(6,8,9)]<-c(rep("DAS.28.Remission.3.months.n",2),c(3,3),c(1,1))
#DAS28remission_big_data_mtxIR2<-na.omit(DAS28remission_big_data_mtxIR2)
DAS28remission_big_data_mtxIR2$month <- sub("^DAS.28.Remission\\.(.*)\\.(n|N)$", "\\1", DAS28remission_big_data_mtxIR2$month)
DAS28remission_big_data_mtxIR2$responders<-as.numeric(DAS28remission_big_data_mtxIR2$responders)
DAS28remission_big_data_mtxIR2$sampleSize<-as.numeric(DAS28remission_big_data_mtxIR2$sampleSize)
DAS28remission_big_data_mtxIR2$DoseLevel01<-as.numeric(DAS28remission_big_data_mtxIR2$DoseLevel01)
DAS28remission_big_data_mtxIR2$DoseLevel01<-as.numeric(DAS28remission_big_data_mtxIR2$Time_point)

DAS28remission_big_data_mtxIR3<- DAS28remission_big_data_mtxIR2 %>%
  mutate(across(c(responders, sampleSize,DoseLevel01, Time_point), as.numeric))

##Drop furst study due to high ROB
DAS28remission_big_data_mtxIR3 <- DAS28remission_big_data_mtxIR3%>%
  filter(!study_name %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012"))
DAS28remission_big_data_mtxIR3$study<-renumber(DAS28remission_big_data_mtxIR3,"study")
row.names(DAS28remission_big_data_mtxIR3)<-NULL
View(DAS28remission_big_data_mtxIR3)

##Data corrections
DAS28remission_big_data_mtxIR3<- DAS28remission_big_data_mtxIR3%>%
  filter(!study_name %in%c("Choy 2012", "Maini 2006", "Maini 2006", "Fleischmann 2023.a", "Liu 2025"))


#Mean and median time frame
mean(as.numeric(DAS28remission_big_data_mtxIR3$TrialDuration))
median(as.numeric(DAS28remission_big_data_mtxIR3$TrialDuration))

library(dplyr)
for(i in 1:nrow(DAS28remission_big_data_mtxIR3)){
  ppp<- master.medtable %>% filter(DAS28remission_big_data_mtxIR3$treatment[i] ==master.medtable$medname )
  DAS28remission_big_data_mtxIR3$treatment[i]<-ppp[1,1]
}


##Creating initial legend
mtxIR.medlegend.DAS28remission <- new.medtable(DAS28remission_big_data_mtxIR3,'treatment')
mtxIR.medlegend.DAS28remission

##NMA plot reason for narrowing of network
nma.plot(DAS28remission_big_data_mtxIR3[,1],DAS28remission_big_data_mtxIR3[,3],data = DAS28remission_big_data_mtxIR3, trtname = as.vector(mtxIR.medlegend.DAS28remission$medname))

################################################################################################################
#trying the narrowing a network function
#############################################################################################
DAS28remission_big_data_mtxIR3

#Narrow 1
mtxIR_narrow<-narrow.all(DAS28remission_big_data_mtxIR3)
mtxIR_narrow


#sub network plot 1
mtxIR.d.DAS28remission<-mtxIR_narrow[[1]]
length(unique(mtxIR.d.DAS28remission$treatment))
mtxIR_subl_DAS28remission<-new.medtable(mtxIR.d.DAS28remission,'treatment')
nma.plot(mtxIR.d.DAS28remission[,1],mtxIR.d.DAS28remission[,3],data = mtxIR.d.DAS28remission, trtname = as.vector(mtxIR_subl_DAS28remission$medname), title = "DAS28 Remission Network plot")


#sub network plot 2
mtxIR.d.DAS28remission2<-mtxIR_narrow[[2]]
length(unique(mtxIR.d.DAS28remission2$treatment))
mtxIR_subl_DAS28remission2<-new.medtable(mtxIR.d.DAS28remission2,'treatment')
nma.plot(mtxIR.d.DAS28remission2[,1],mtxIR.d.DAS28remission2[,3],data = mtxIR.d.DAS28remission2, trtname = as.vector(mtxIR_subl_DAS28remission2$medname), title = "DAS28 Remission Network plot")

##############################################################################
##############################################################################
#Performing sub network plot 1
library(rjags)
library(coda)
library(gemtc)
View(mtxIR.d.DAS28remission)
mtxIR_subl_DAS28remission

# network_data_DAS28remission_mtxIR_1<-mtxIR.d.DAS28remission%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
#            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
network_data_DAS28remission_mtxIR_1<-mtxIR.d.DAS28remission%>%
  select(c("study", "study_name","treatment","responders", "sampleSize", "month"))
network_data_DAS28remission_mtxIR_1$treatment<-renumber(network_data_DAS28remission_mtxIR_1,"treatment")
##To download DAS28remission data
network_data_DAS28remission_mtxIR_data_download<-network_data_DAS28remission_mtxIR_1
network_data_DAS28remission_mtxIR_1$responders<-as.numeric(network_data_DAS28remission_mtxIR_1$responders)
network_data_DAS28remission_mtxIR_1$sampleSize<-as.numeric(network_data_DAS28remission_mtxIR_1$sampleSize)
network_data_DAS28remission_mtxIR_data_download$treatment<-selected_treatment_name(network_data_DAS28remission_mtxIR_data_download$treatment,mtxIR_subl_DAS28remission)
write.csv(network_data_DAS28remission_mtxIR_data_download,"DAS28remission MTX-IR dataset_AUG2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("Disease Activity Score 28 Remission (DAS28 Remission) MTX-IR-Network plot with named interventions_AUG2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(network_data_DAS28remission_mtxIR_1$study,network_data_DAS28remission_mtxIR_1$treatment,data = network_data_DAS28remission_mtxIR_1, trtname = as.vector(mtxIR_subl_DAS28remission$medname),title = "Disease Activity Score 28 MTX-IR Remission Network plot")
# Close the pdf file
dev.off() 
# 1. Open jpeg file
jpeg("Disease Activity Score 28 Remission (DAS28 Remission) MTX-IR-Network plot without named interventions_JUNE2026_V2.jpg", width = 3050, height = 2050, res=300)
nma.plot(network_data_DAS28remission_mtxIR_1$study,network_data_DAS28remission_mtxIR_1$treatment,data = network_data_DAS28remission_mtxIR_1,title = "Disease Activity Score 28 MTX-IR Remission Network plot")
# Close the pdf file
dev.off() 


network_data_DAS28remission_mtxIR_1<-mtc.network(network_data_DAS28remission_mtxIR_1)

##RUN the model

set.seed(5014)
network_data_DAS28remission_mtxIR_1.model <- mtc.model(network_data_DAS28remission_mtxIR_1, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_DAS28remission_mtxIR_1.model$inits<- mapply(c, network_data_DAS28remission_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_DAS28remission<-mtc.run(network_data_DAS28remission_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_DAS28remission)

View(relative.effect.table.OR(results_mtxIR_DAS28remission))
refOR_mtxIR_DAS28remission<-as.data.frame(relative.effect.table.OR2(results_mtxIR_DAS28remission, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refOR_mtxIR_DAS28remission)<-mapvalues(row.names(refOR_mtxIR_DAS28remission),
                                               from = c(1:(nrow(mtxIR_subl_DAS28remission))),
                                               to = as.vector(mtxIR_subl_DAS28remission$medname))
colnames(refOR_mtxIR_DAS28remission)<-mapvalues(colnames(refOR_mtxIR_DAS28remission),
                                              from = c(1:(nrow(mtxIR_subl_DAS28remission))),
                                              to = as.vector(mtxIR_subl_DAS28remission$medname))
View(refOR_mtxIR_DAS28remission)
# which(row(refOR)==col(refOR))
# refOR[row(refOR)==col(refOR)]
index_diagonal2<-data.frame(index1=1:nrow(refOR_mtxIR_DAS28remission),index2=1:nrow(refOR_mtxIR_DAS28remission))
options(warn=-1)
for(i in 1:nrow(refOR_mtxIR_DAS28remission)){
  refOR_mtxIR_DAS28remission[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_mtxIR_DAS28remission[row(refOR_mtxIR_DAS28remission)==col(refOR_mtxIR_DAS28remission)])[i],mtxIR_subl_DAS28remission)
}
options(warn=0)
View(refOR_mtxIR_DAS28remission)
#library(xlsx)
#write.xlsx(refOR_mtxIR_DAS28remission,"Results for mtxIRlogics_OCT31.xlsx",sheetName = "_DAS28remission",row.names = T)
#write.xlsx(refOR_mtxIR,"Results for mtxIRlogics Informative Prior.xlsx",sheetName = "ACR50",row.names = T)
write.csv(refOR_mtxIR_DAS28remission,"DAS28 Remission response pairwise comparison MTX-IR_AUG2026.csv", row.names = T)

#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "DAS28 Remission")
writeData(wb, "DAS28 Remission", refOR_mtxIR_DAS28remission, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)

#Converting NMA estimates to OR in a dataframe format
nma_estimate_mtxIR_DAS28remission_dataframe<-nma.estimate.table.OR(results_mtxIR_DAS28remission)
View(nma_estimate_mtxIR_DAS28remission_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_DAS28remission_mtxIR_1.ns<-mtc.nodesplit(network_data_DAS28remission_mtxIR_1, thin=50,likelihood="binom", link="logit",hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 

names(network_data_DAS28remission_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_DAS28remission_mtxIR_1)

# Overall summary and plot
summary.network_data_DAS28remission_mtxIR_1.ns <- summary(network_data_DAS28remission_mtxIR_1.ns)

#Direct effects
direct.data<-summary.network_data_DAS28remission_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_DAS28remission_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_DAS28remission_mtxIR_1.ns$cons.effect

#P-values
P_values_nodesplit_das28Remission_mtxIR<-summary.network_data_DAS28remission_mtxIR_1.ns$p.value
P_values_nodesplit_das28Remission_mtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_das28Remission_mtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_das28Remission_mtxIR2<-P_values_nodesplit_das28Remission_mtxIR
P_values_nodesplit_das28Remission_mtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_das28Remission_mtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_das28Remission_mtxIR2$comb<-paste0(P_values_nodesplit_das28Remission_mtxIR2$t1,":",P_values_nodesplit_das28Remission_mtxIR2$t2)



direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(exp(direct.data[,c(3:5)]),2),round(exp(indirect.data[,c(3:5)]),2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","OR1","low1","high1","OR2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                round(1/direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","OR1","low1","high1","OR2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)
grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_DAS28remission_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)

#Getting the quality assessment for each study
###Adding ROB data
mtxIR.d.DAS28remission[,c(1:11)]
ffff_mtxIR<-unique(mtxIR.d.DAS28remission$study_name)
data_qq_mtxIR<-mtxIR.d.DAS28remission[,-c(10:11)]
View(ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff_mtxIR,])
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff_mtxIR,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<-current_assessment_for_studies_mtxIR%>%
  select(c(col_names))

list_appraisal_mtxIR<-list()
for(i in 1:length(ffff_mtxIR)){
  temp11<-data_qq_mtxIR[data_qq_mtxIR$study_name==ffff_mtxIR[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==ffff_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal_mtxIR[[i]]<-temp13
}
DAS28remission_mtxIR_withROB<-do.call(rbind,list_appraisal_mtxIR)
View(DAS28remission_mtxIR_withROB)
#ACR50_mtxIR_withROB2<-ACR50_mtxIR_withROB[,-c(10:22)]
DAS28remission_mtxIR_withROB2<-DAS28remission_mtxIR_withROB%>%
  select(c("study","study_name", "treatment", "responders", "sampleSize",
           "month", "DoseLevel01", "TrialDuration","Time_point","OVERALL..ALL.OTHER.OUTCOMES"))
View(DAS28remission_mtxIR_withROB2)

###Adding Indirectness data
#ACR50_mtxIR_withROB2_indrectness<-cbind(ACR50_mtxIR_withROB2,Indirectness="low")
study_unique<-unique(DAS28remission_mtxIR_withROB2$study_name)
data_nn<-DAS28remission_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp[!is.na(resolved_indirectness_consensus_temp$`Final Rating`),]
View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]
row.names(current_indirect_assessment_for_studies)<-NULL
##because of bari 4mg
# current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:17),],current_indirect_assessment_for_studies[c(17),],current_indirect_assessment_for_studies[c(18:28),])
# current_indirect_assessment_for_studies$intervention[18]<-"Baricitinib (4 mg)"
#current_indirect_assessment_for_studies<-rbind(current_indirect_assessment_for_studies[c(1:19),],current_indirect_assessment_for_studies[c(19),],current_indirect_assessment_for_studies[c(20:30),])
#current_indirect_assessment_for_studies$intervention[20]<-"Baricitinib (4 mg)"
list_appraisal2<-list()
for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier`==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
DAS28remission_try_data_new_withROB2_indrectness<-do.call(rbind,list_appraisal2)
colnames(DAS28remission_try_data_new_withROB2_indrectness)[11]<-"Indirectness"

DAS28remission_try_data_new_withROB2_indrectness$Indirectness<-tolower(DAS28remission_try_data_new_withROB2_indrectness$Indirectness)
View(DAS28remission_try_data_new_withROB2_indrectness)


DAS28remission_mtxIR_withROB2_indrectness<-DAS28remission_try_data_new_withROB2_indrectness
#####Formatting the month column
DAS28remission_mtxIR_withROB2_indrectness$month<-sub(".months","",DAS28remission_mtxIR_withROB2_indrectness$month)
DAS28remission_mtxIR_withROB2_indrectness$month<-as.numeric(ifelse(DAS28remission_mtxIR_withROB2_indrectness$month=="End.of.trial.n",DAS28remission_mtxIR_withROB2_indrectness$TrialDuration,DAS28remission_mtxIR_withROB2_indrectness$month))

##Taking Huffstutter 2017 out, it does not allow the model to converge because of zero cells
#DAS28remission_mtxIR_withROB2_indrectness<-DAS28remission_mtxIR_withROB2_indrectness[!DAS28remission_mtxIR_withROB2_indrectness$study_name%in%c("Huffstutter 2017","Schiff 2014","Cohen 2006"),]
View(DAS28remission_mtxIR_withROB2_indrectness)


DAS28remission_mtxIR_data_new_withROB2_indrectness2<-DAS28remission_mtxIR_withROB2_indrectness
DAS28remission_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(DAS28remission_mtxIR_data_new_withROB2_indrectness2,"treatment")
#ACR50_data_for_analysis2$OVERALL..ALL.OTHER.OUTCOMES<-ifelse(ACR50_data_for_analysis2$OVERALL..ALL.OTHER.OUTCOMES=="unclear","some concerns",ACR50_data_for_analysis2$OVERALL..ALL.OTHER.OUTCOMES)
DAS28remission_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES<-ifelse(DAS28remission_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES=="unclear","some concerns",DAS28remission_mtxIR_data_new_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES)
pairwise_DAS28remission_OR<-pairwise(treat=treatment,event = responders,data=network_data_DAS28remission_mtxIR_1$data.ab,
                                     studlab = study,n=sampleSize,sm="OR")

netPairwiseDAS28remission<-netmeta(pairwise_DAS28remission_OR)
netPairwiseDAS28remission_split<-netsplit(netPairwiseDAS28remission)
View(netPairwiseDAS28remission_split$direct.common)
View(netPairwiseDAS28remission_split$indirect.common)
netPairwiseDAS28remission_contribmat<-netcontrib(netPairwiseDAS28remission, method = "r")
View(netPairwiseDAS28remission_contribmat$random)
direct_evidence_new<-cbind(netPairwiseDAS28remission_split$direct.common[,c(1:5)],netPairwiseDAS28remission_split$indirect.common[,c(2:5)])


View(direct_evidence_new)
#getting only direct estimates
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]
library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")
direct_only4%in%paste(nma_estimate_mtxIR_DAS28remission_dataframe$treat1,nma_estimate_mtxIR_DAS28remission_dataframe$treat2,sep = ":")

#####Renaming the number treatment to named treatment
selected_treatment_name<-function(dataset,sub.medtable){
  suppressWarnings(library(dplyr))
  for(i in 1:length(dataset)){
    selected_treatment_name<- sub.medtable %>% filter(dataset[i] ==sub.medtable$new_mednum )
    dataset[i]<-selected_treatment_name[1,2]
  }
  return(dataset)
}

#Setting up initial grade table
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]
##Getting direct evidence for intervention in only multi-arm trials that fails to get direct estimate
all_actual_direct_evidence<-t(apply(apply(str_split_fixed(colnames(netPairwiseDAS28remission_contribmat$random), ":", 2),2,as.numeric),1,switch_func))
colnames(all_actual_direct_evidence)<-c("treat1", "treat2")
#Not present in dmard mtxIR datatset for acr50
inthreeArmtrial<-anti_join(as.data.frame(all_actual_direct_evidence), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]



grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_subl_DAS28remission)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_subl_DAS28remission)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))


#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
DAS28remission_output_df2_list<-list()
I_square_list<-list()

DAS28remission_mtxIR_data_new_withROB2_indrectness2$responders<-as.numeric(DAS28remission_mtxIR_data_new_withROB2_indrectness2$responders)
DAS28remission_mtxIR_data_new_withROB2_indrectness2$sampleSize<-as.numeric(DAS28remission_mtxIR_data_new_withROB2_indrectness2$sampleSize)
for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-DAS28remission_mtxIR_data_new_withROB2_indrectness2[DAS28remission_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:5)
  
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select("study", "responders.1", "sampleSize.1","responders.2", "sampleSize.2" )
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"OVERALL..ALL.OTHER.OUTCOMES"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","n.t","x.c","n.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    model.file<-"Dichot_direct_onestudy.txt"
    n.chains=3
    set.seed(77453)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(#p.c=c(0.5),p.t=c(0.5),
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "n.t"=temp12$n.t,"x.c"=temp12$x.c,"n.c"=temp12$n.c),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("OR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(OR=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
  }
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    model.file<-"Dichot_direct_multiplestudy_RE.txt"
    n.chains=3
    #model.inits <- list(p.c=c(rep(0.5,length(temp12$NS))),p.t=c(rep(0.5,length(temp12$NS))))
    set.seed(433663)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "n.t"=temp12$n.t,"x.c"=temp12$x.c,"n.c"=temp12$n.c, "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("OR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select("study", "responders.1", "sampleSize.1","responders.2", "sampleSize.2" )
  #nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(Studies_contributing_to_direct_evidence_wide$study),c(1,4,5,13,14)]
  nma_data50_4<-nma_data50_3
  colnames(nma_data50_4)<-c("NS","x.t","n.t","x.c","n.c")
  row.names(nma_data50_4)<-NULL
  
  m.bin<-metabin(event.e = x.t, 
                 n.e = n.t,
                 event.c = x.c,
                 n.c = n.c,
                 # studlab = study_name,
                 data = nma_data50_4,
                 sm = "OR",
                 method = "MH",
                 MH.exact = TRUE,
                 common = FALSE,
                 random = TRUE,
                 method.tau = "PM",
                 method.random.ci= TRUE,
                 title = "DAS28Remission"
                 
  )
  
  summary_m.bin<-summary(m.bin)
  contribution_weight<-as.character(round((summary_m.bin$w.random/sum(summary_m.bin$w.random))*100,2))
  pooled_estimate2<-data.frame(OR=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.bin$tau2,4)
                          ,
                          Lower=round(summary_m.bin$lower.tau2,4),
                          Upper=round(summary_m.bin$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.bin$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.bin$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.bin$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.bin$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  
  
  #}
  
  OR_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1)
  base_data <- tibble(mean  = OR_individual_studies$OR, 
                      lower = OR_individual_studies$Lower,
                      upper = OR_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      # treatment1 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      # treatment2 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      treatment1 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      treatment2 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      OR = paste(format(round(OR_individual_studies$OR,2),nsmall=2),paste("(",format(round(OR_individual_studies$Lower,2),nsmall=2),", ",format(round(OR_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$OR, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   OR = paste(format(round(pooled_estimate2$OR,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("DAS28 Remission", "Study"),
                   # treatment1 = c("",selected_treatment_name(unique(temp11$treatment.1),subl)),
                   # treatment2 = c("",selected_treatment_name(unique(temp11$treatment.2),subl)),
                   treatment1 = c("",selected_treatment_name(unique(temp11$treatment.2),mtxIR_subl_DAS28remission)),
                   treatment2 = c("",selected_treatment_name(unique(temp11$treatment.1),mtxIR_subl_DAS28remission)),
                   #ROB=c(paste("Favours",selected_treatment_name(unique(temp11$treatment.1),subl)),"ROB"),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   # OR = c(paste("Favours",selected_treatment_name(unique(temp11$treatment.2),subl)), "OR"),
                   OR = c("", "OR"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  DAS28remission_output_df <- bind_rows(header,
                                        base_data,
                                        empty_row,
                                        pooled)
  
  
  
  
  DAS28remission_output_df2<-rbind(DAS28remission_output_df,DAS28remission_output_df[c(nrow(DAS28remission_output_df)),])
  DAS28remission_output_df2[nrow(DAS28remission_output_df),1]<-I2
  DAS28remission_output_df2[nrow(DAS28remission_output_df),2]<-NA
  DAS28remission_output_df2[nrow(DAS28remission_output_df),c(6:10)]<-c(NA,NA,NA,NA,NA)
  
  DAS28remission_output_df2_list[[i]]<-DAS28remission_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-DAS28remission_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, OR,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(0.1, 2.4), 
               xlog = TRUE, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               # grid = structure(c(1.2), 
               #                  gp = gpar(lty = 2, col = "#F73831")),
               #txt_gp = fpTxtGp(cex=0.85),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               # xlab = paste("Outcome: ACR50-Direct Comparison",selected_treatment_name(unique(temp11$treatment.1),subl),"and",selected_treatment_name(unique(temp11$treatment.2),subl))
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_subl_DAS28remission),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_subl_DAS28remission))
    )
  
}

length(graph_list)
graph_list[[3]]
graph_list[[13]]
DAS28remission_output_df2_list[[7]]
n_studies<-c()
n_participants<-c()
for(i in 1:length(graph_list)){
  df<-DAS28remission_output_df2_list[[i]][grep("\\d+/\\d+", DAS28remission_output_df2_list[[i]]$treatment1),c(2,3) ]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Select denominator values
  df$col1_after_slash <- sapply(strsplit(df$treatment1, "/"), "[[", 2)
  df$col2_after_slash <- sapply(strsplit(df$treatment2, "/"), "[[", 2)
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(unlist(df[,c(3,4)]))))
  
}

Grade_assessment_dataset_direct_DAS28remission_mtxIR<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_subl_DAS28remission),
                                                               selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_subl_DAS28remission),
                                                               Outcome="DAS28 Remission",
                                                               `Effect measure`="OR",
                                                               `Number of studies`=n_studies,
                                                               `Number of participant`=n_participants,
                                                               ROB="",Indirectness="",Imprecision="",
                                                               Inconsistency="",
                                                               Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_DAS28remission_mtxIR)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                                      "Intevention",
                                                                                      "Effect measure",
                                                                                      "Number of studies",
                                                                                      "Number of participant",
                                                                                      "Risk of Bias",
                                                                                      "Publication Bias",
                                                                                      "Certainty of evidence (ignores imprecision)")
View(Grade_assessment_dataset_direct_DAS28remission_mtxIR)
#Saving grade assessment details for DAS28remission mtxIR as an object of lists
grade_details_DAS28remission_mtxIR<-list(graph_list,DAS28remission_output_df2_list,Grade_assessment_dataset_direct_DAS28remission_mtxIR,I_square_list)
saveRDS(grade_details_DAS28remission_mtxIR,file="List of DAS28remission MTX-IR data components for GRADE_AUG2026.rds")

#########
DAS28Remission_output_one_study_list <- list()
DAS28Remission_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

# Count studies for each comparison
n_studies <- sapply(DAS28remission_output_df2_list, function(x) {
  
  df <- x[
    grep("\\d+/\\d+", x$treatment1),
    c(2,3)
  ]
  
  nrow(df)
})

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
DAS28Remission_output_one_study_list  <- DAS28remission_output_df2_list[one_idx]
DAS28Remission_output_more_than_2_list <- DAS28remission_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_DAS28remission_mtxIR[
    Grade_assessment_dataset_direct_DAS28remission_mtxIR$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_DAS28remission_mtxIR[
    Grade_assessment_dataset_direct_DAS28remission_mtxIR$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_DAS28Remission_mtxIR_one_study <- list(
  graph_list_one_study_list,
  DAS28Remission_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_DAS28Remission_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  DAS28Remission_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_DAS28Remission_mtxIR_one_study,
  file = "List of DAS28Remission MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_DAS28Remission_mtxIR_more_than_2,
  file = "List of DAS28Remission MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)

#####################################
#####################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_subl_DAS28remission)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct_DAS28remission_mtxIR
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_DAS28remission_mtxIR$Treatment, 
                                                Grade_assessment_dataset_direct_DAS28remission_mtxIR$Intevention), 
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
#Getting the actual direct estimate for missing three arm trial
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(DAS28remission_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","OR"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"

#############
##Creating indirect assessment 
####Get direct evidence report 
# setwd("/Users/mujaab/Documents/OneDrive - University of Calgary/Manucripts for PhD papers/Creating various datasets for GRADE evaluation")
# kamso_direct_evidence_rating_acr50<-read.csv("mtxIR-IR_ACR50_directAssessment_April2023.csv")
kamso_direct_evidence_rating_DAS28remission<-read.csv("Glen_2026-09-03 22_14_DAS28Rem-MTXIR_complete.csv")
kamso_direct_evidence_rating_DAS28remission2<-kamso_direct_evidence_rating_DAS28remission[,-c(13:15)]
colnames(kamso_direct_evidence_rating_DAS28remission2)[1:2]<-colnames(grade_table4)[1:2]

kamso_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],kamso_direct_evidence_rating_DAS28remission2,
                                     by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)


kamso_direct_evidence_rating4<-cbind(kamso_direct_evidence_rating3[,c(1,2,4:7)],
                                     kamso_direct_evidence_rating3[,c(3,8:13)])
#View(kamso_direct_evidence_rating4[order(kamso_direct_evidence_rating4$Outcome),])##Perhaps might use this in the long run
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Outcome),"Outcome"]<-"DAS28 Remission"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"OR"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
kamso_direct_evidence_rating4[is.na(kamso_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"

kamso_direct_evidence_rating5<-cbind(kamso_direct_evidence_rating4,
                                     num1=as.numeric(changing_treatment_name2num(kamso_direct_evidence_rating4$treat1,mtxIR_subl_DAS28remission)),
                                     num2=as.numeric(changing_treatment_name2num(kamso_direct_evidence_rating4$treat2,mtxIR_subl_DAS28remission)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}

##ratings with comb
kamso_direct_evidence_rating5[,14:15]<-t(apply(kamso_direct_evidence_rating5[,14:15],1,switch_func2))
kamso_direct_evidence_rating6<-cbind(kamso_direct_evidence_rating5,
                                     comb=paste(kamso_direct_evidence_rating5[,14],
                                                kamso_direct_evidence_rating5[,15],sep = ":"))
kamso_direct_evidence_rating6_onlyDirect<-kamso_direct_evidence_rating6[!is.na(kamso_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
kamso_direct_evidence_rating6$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
##Getting indirect and nma estimates
Temp_GRADE<-merge(kamso_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(kamso_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("3:13","3:5"),"Indirect"]<-"NE"
#write.csv(Temp_GRADE4,"Initial GRADE table for ACR50_APR24.csv",row.names = F)
write.csv(Temp_GRADE4,"Initial GRADE table for DAS28remission_MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for DAS28remission_MTX-IR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR DAS28remission OUTCOME
initial_gradeTable_DAS28remission<-Temp_GRADE5

estimate_contribution_matrix<-netPairwiseDAS28remission_contribmat$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_DAS28remission2<-initial_gradeTable_DAS28remission[match(all_nma_estimates,initial_gradeTable_DAS28remission$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()
#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_DAS28remission2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames

initial_gradeTable_DAS28remission2.1<-cbind(initial_gradeTable_DAS28remission2[,c(1:14)],first_order_loop_contri,
                                            initial_gradeTable_DAS28remission2[,-c(1:14)])

initial_gradeTable_DAS28remission2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_DAS28remission2.1,
                                                                       all_nma_estimates,
                                                                       estimate_contribution_matrix,
                                                                       kamso_direct_evidence_rating6_onlyDirect,
                                                                       top_n = 2))

View(initial_gradeTable_DAS28remission2.1)
# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_DAS28remission2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_DAS28remission2.1<-result_df
initial_gradeTable_DAS28remission2.1$combI<-str_replace(initial_gradeTable_DAS28remission2.1$comb, ":", " vs ")

initial_gradeTable_DAS28remission2.1<-process_dataset(initial_gradeTable_DAS28remission2.1, 2)


initial_gradeTable_DAS28remission3<-initial_gradeTable_DAS28remission2.1[match(initial_gradeTable_DAS28remission$comb,initial_gradeTable_DAS28remission2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for DAS28Remission_MTX-IR outcome_AUG2026.csv",row.names = T)

###########Assessing Intransitivity
sub_data_intransitivity<-initial_gradeTable_DAS28remission3[initial_gradeTable_DAS28remission3$Indirect!="-",]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_DAS28remission3[initial_gradeTable_DAS28remission3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_DAS28remission3.1<-cbind(initial_gradeTable_DAS28remission3[,c(1:23)],
                                            Intransitivity="",
                                            `Certainty of evidence (ignores imprecision)-Indirect`="",
                                            initial_gradeTable_DAS28remission3[,-c(1:23)])
initial_gradeTable_DAS28remission3.1[initial_gradeTable_DAS28remission3.1$Indirect!="-",c("Intransitivity",
                                                                                          "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                                                "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_DAS28remission4<-initial_gradeTable_DAS28remission3.1

###Applying the P-values
P_values_nodesplit_das28Remission_mtxIR2
grade_temp<-cbind(initial_gradeTable_DAS28remission4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_DAS28remission4[,-c(1:25)])
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_das28Remission_mtxIR3 <- P_values_nodesplit_das28Remission_mtxIR2[match(drt$comb, P_values_nodesplit_das28Remission_mtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_das28Remission_mtxIR3$p,2)

#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)

subDAS28RemissionSOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatDAS28Remission<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subDAS28RemissionSOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(nrow(subcontriMatDAS28Remission))
IndirPer <- numeric(nrow(subcontriMatDAS28Remission))

# Loop to calculate dirPer and IndirPer
for (i in 1:nrow(subcontriMatDAS28Remission)) {
  dirPer[i] <- subcontriMatDAS28Remission[i, colnames(subcontriMatDAS28Remission) == row.names(subcontriMatDAS28Remission)[i]]
  IndirPer[i] <- 1 - dirPer[i]
  #IndirPer[i] <- sum(subcontriMatDAS28Remission[i, colnames(subcontriMatDAS28Remission) != row.names(subcontriMatDAS28Remission)[i]][-1])
}

# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")


grade_temp2[match(row.names(subcontriMatDAS28Remission), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(row.names(subcontriMatDAS28Remission), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)


##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
#When both direct and indirect are present

fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg

####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho


initial_gradeTable_DAS28remission5<-grade_temp2.2

initial_gradeTable_DAS28remission5<-initial_gradeTable_DAS28remission5[,-c(30,33)]
colnames(initial_gradeTable_DAS28remission5)[c(29,30)]<-c("Initial NMA evidence rating",
                                                          "P-Value")
###After incoherence evaluation
initial_gradeTable_DAS28remission5.1<-cbind(initial_gradeTable_DAS28remission5[,c(1:31)],
                                            qoeAfterIncohenrence="",
                                            initial_gradeTable_DAS28remission5[,-c(1:31)])
initial_gradeTable_DAS28remission5.1[(initial_gradeTable_DAS28remission5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_DAS28remission5.1[(initial_gradeTable_DAS28remission5.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_DAS28remission5.1[initial_gradeTable_DAS28remission5.1$Incoherence!="",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_DAS28remission5.1[(initial_gradeTable_DAS28remission5.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence

#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_DAS28remission5.1[(initial_gradeTable_DAS28remission5.1$Direct!="-")&(initial_gradeTable_DAS28remission5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_DAS28remission5.1[(initial_gradeTable_DAS28remission5.1$Direct!="-")&(initial_gradeTable_DAS28remission5.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_DAS28remission6<-cbind(initial_gradeTable_DAS28remission5.1,
                                          p.estimate=apply(as.data.frame(initial_gradeTable_DAS28remission5.1$NMA),1,function(x){
                                            as.numeric(sub(" .*", "", x))
                                          }),
                                          t(apply(as.data.frame(initial_gradeTable_DAS28remission5.1$NMA),1,function(x){
                                            as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                            
                                          })))

colnames(initial_gradeTable_DAS28remission6)[c(36,37)]<-c("l.bound","u.bound")


##Getting baseline risk using the placebo intervention
DAS28remission_data<-network_data_DAS28remission_mtxIR_1$data.ab%>%
  select(c("study","study_name", "treatment", "responders", "sampleSize"))
DAS28remission_data_wide2<-reshape(number_of_arms(DAS28remission_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(DAS28remission_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(DAS28remission_data_wide2[,c(3,7,11)],1,xcv)

DAS28remission_data_wide3<-cbind(DAS28remission_data_wide2,na=numberOfArms)
View(DAS28remission_data_wide3)
colnames(DAS28remission_data_wide3)<-c("study", "study_name", "t.1","x.1","n.1",
                                       "t.2","x.2","n.2","t.3","x.3","n.3",
                                       "na" )
#######
MTCData_dichot <- as.data.frame(DAS28remission_data_wide3)
x <- data.matrix(MTCData_dichot[,c(4,7,10)])
n <- data.matrix(MTCData_dichot[,c(5,8,11)])
t <- data.matrix(MTCData_dichot[,c(3,6,9)])
na <- data.matrix(MTCData_dichot[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <-length(t[which(t==1)])
xPBO<-x[which(t==1)]
nPBO <- n[which(t==1)]


# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Wichmann-Hill", .RNG.seed=1)

#loading rate model txt file
model.file4<-"Dichot_multiarm_RE.txt"
n.chains=4
model.inits <- initv1 

dichot.model3<-jags.model(file=model.file4,
                          data= list("x"=x,"t"=t,"na" = na[,1], "NT" = NT, "NS" = NS, "n"=n,
                                     "NTpbo"= NTpbo,"xPBO"=xPBO,"nPBO"=nPBO ),
                          n.chains=n.chains,
                          inits = model.inits,
                          n.adapt=10000
)

update(object=dichot.model3,n.iter=20000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("OR","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( dichot.model3, params4, n.iter=50000 )
summary_samps<-summary(samps4)
View(summary_samps$statistics)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]
typeof(estimates_relative_to_baseline)

#Converting estimates relative to baseline to absolutes
#nma_estimate_mtxIR_dataframe3<-as.data.frame(apply(t(do.call(rbind,nma_estimate_mtxIR_dataframe2)),2,as.numeric))
estimates_relative_to_baseline<-initial_gradeTable_DAS28remission6[1:(NT-1),c(1,2,35:37)]

absolutes_from_baseline<-cbind(baseline_risk,(estimates_relative_to_baseline[,3:5]*baseline_risk)/((1-baseline_risk)+(estimates_relative_to_baseline[,3:5]*baseline_risk)))
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(baseline_risk,estimates_relative_to_baseline_combine[,7]))

colnames(all_baseline)<-"base"
com.list<-list()
unique_treat2<-unique(initial_gradeTable_DAS28remission6$treat2)
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_DAS28remission6[initial_gradeTable_DAS28remission6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,(d.sub[,3:5]*base)/((1-base)+(d.sub[,3:5]*base)))
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9]*1000)
more_fewer_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"fewer"
  }else{
    text<-"more"
  }
  return(paste(round(abs(e1-e2),0),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,13]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,13])))

prefix.assumedrisk <- ''
prefix.corrrisk<-''
suffix.assumedrisk <- ' per 1000'
suffix.corrrisk <- ' per 1000'
result_com.dataframe_withtext<-cbind(com.dataframe_withtext[,1:2],
                                     paste(prefix.assumedrisk,round(com.dataframe_withtext[,6]*1000,0), suffix.assumedrisk, sep=''),
                                     paste(prefix.corrrisk, round(com.dataframe_withtext[,7]*1000,0), suffix.corrrisk, ' (', round(com.dataframe_withtext[,8]*1000,0), ' to ',round(com.dataframe_withtext[,9]*1000,0),')', sep=''),
                                     paste(prefix.corrrisk, com.dataframe_withtext[,13], ' (', com.dataframe_withtext[,14], ' to ',com.dataframe_withtext[,15],')', sep=''),
                                     # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                     com.dataframe_withtext[,10:15]
)

colnames(result_com.dataframe_withtext)[c(3:11)]<-c("Assumed risk (comparator)",	"Corresponding risk (intervention)",
                                                    "Risk difference","p.estimate","l.bound","u.bound","","","")


initial_gradeTable_DAS28remission7<-cbind(initial_gradeTable_DAS28remission6[,-c(35:37)],
                                          result_com.dataframe_withtext[,-c(1,2)],pair="")
#colnames(initial_gradeTable_DAS28remission7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_DAS28remission8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_DAS28remission7$treat1,mtxIR_subl_DAS28remission)),
                                          `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_DAS28remission7$treat2,mtxIR_subl_DAS28remission)),
                                          initial_gradeTable_DAS28remission7)
colnames(initial_gradeTable_DAS28remission8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_DAS28remission8<-initial_gradeTable_DAS28remission8[,-c(46)]

#write.csv(initial_gradeTable_ACR508,"ACR50 SoF table without imprecision assessment May2023.csv",row.names = F)
write.csv(initial_gradeTable_DAS28remission8,"DAS28remission MTX-IR SoF table without imprecision assessment AUG2026.csv",row.names = F)

effect.size.rating_list<-lapply(split(initial_gradeTable_DAS28remission7, seq_len(nrow(initial_gradeTable_DAS28remission7))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 50, large.effect = 200, framework = T)
})

as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_DAS28remission9<-initial_gradeTable_DAS28remission8
initial_gradeTable_DAS28remission9<-cbind(initial_gradeTable_DAS28remission9[,c(1:34,37:45,35:36)])
initial_gradeTable_DAS28remission9.1<-cbind(initial_gradeTable_DAS28remission9[,c(1:37)],
                                            Small.effect.threshold= 50,
                                            Large.effect.threshold= 200,
                                            Framework = "NULL",
                                            Imprecision.1 = "",
                                            Certainty.of.evidence = "",
                                            initial_gradeTable_DAS28remission9[,c(44,45)])



initial_gradeTable_DAS28remission9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_DAS28Remission_mtxIR<-initial_gradeTable_DAS28remission9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_DAS28Remission_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_DAS28Remission_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_DAS28Remission_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_DAS28Remission_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_DAS28Remission_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_DAS28Remission_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_DAS28Remission_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_DAS28Remission_mtxIR[is.na(detailed_SOF_DAS28Remission_mtxIR$Certainty.of.evidence),])
detailed_SOF_DAS28Remission_mtxIR[is.na(detailed_SOF_DAS28Remission_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

#detailed_SOF_DAS28Remission_mtxIR$comb<-initial_gradeTable_DAS28remission7$comb

write.csv(detailed_SOF_DAS28Remission_mtxIR,"DAS28Remission_MTX-IR_Detail_SoF_AUG2026.csv", row.names = F)




#Performing sub network plot 2
View(mtxIR.d.DAS28remission2)
mtxIR_subl_DAS28remission2

# network_data_DAS28remission_mtxIR_1<-mtxIR.d.DAS28remission%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
#            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
network_data_DAS28remission_mtxIR_12<-mtxIR.d.DAS28remission2%>%
  select(c("study", "study_name","treatment","responders", "sampleSize", "month"))
network_data_DAS28remission_mtxIR_12$treatment<-renumber(network_data_DAS28remission_mtxIR_12,"treatment")
##To download DAS28remission data
network_data_DAS28remission_mtxIR_data_download2<-network_data_DAS28remission_mtxIR_12
network_data_DAS28remission_mtxIR_12$responders<-as.numeric(network_data_DAS28remission_mtxIR_12$responders)
network_data_DAS28remission_mtxIR_12$sampleSize<-as.numeric(network_data_DAS28remission_mtxIR_12$sampleSize)
network_data_DAS28remission_mtxIR_data_download2$treatment<-selected_treatment_name(network_data_DAS28remission_mtxIR_data_download2$treatment,mtxIR_subl_DAS28remission2)
write.csv(network_data_DAS28remission_mtxIR_data_download2,"DAS28remission MTX-IR sub-network dataset_AUG2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("Disease Activity Score 28 Remission (DAS28 Remission) MTX-IR-Network plot sub-network with named interventions_AUG2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(network_data_DAS28remission_mtxIR_12$study,network_data_DAS28remission_mtxIR_12$treatment,data = network_data_DAS28remission_mtxIR_12, trtname = as.vector(mtxIR_subl_DAS28remission2$medname),title = "Disease Activity Score 28 MTX-IR Remission Network plot")
# Close the pdf file
dev.off() 
# 1. Open jpeg file
jpeg("Disease Activity Score 28 Remission (DAS28 Remission) MTX-IR-Network plot sub-network without named interventions_AUG2026_V2.jpg", width = 3050, height = 2050, res=300)
nma.plot(network_data_DAS28remission_mtxIR_12$study,network_data_DAS28remission_mtxIR_12$treatment,data = network_data_DAS28remission_mtxIR_12,title = "Disease Activity Score 28 MTX-IR Remission Network plot")
# Close the pdf file
dev.off() 


network_data_DAS28remission_mtxIR_12<-mtc.network(network_data_DAS28remission_mtxIR_12)

##RUN the model

set.seed(5014)
network_data_DAS28remission_mtxIR_1.model2 <- mtc.model(network_data_DAS28remission_mtxIR_12, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_DAS28remission_mtxIR_1.model2$inits<- mapply(c, network_data_DAS28remission_mtxIR_1.model2$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_DAS28remission2<-mtc.run(network_data_DAS28remission_mtxIR_1.model2,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_DAS28remission2)

View(relative.effect.table.OR(results_mtxIR_DAS28remission2))
refOR_mtxIR_DAS28remission2<-as.data.frame(relative.effect.table.OR2(results_mtxIR_DAS28remission2, comparator = "3"))

#library for mapvalues
library(plyr)
row.names(refOR_mtxIR_DAS28remission2)<-mapvalues(row.names(refOR_mtxIR_DAS28remission2),
                                                 from = c(1:(nrow(mtxIR_subl_DAS28remission2))),
                                                 to = as.vector(mtxIR_subl_DAS28remission2$medname))
colnames(refOR_mtxIR_DAS28remission2)<-mapvalues(colnames(refOR_mtxIR_DAS28remission2),
                                                from = c(1:(nrow(mtxIR_subl_DAS28remission2))),
                                                to = as.vector(mtxIR_subl_DAS28remission2$medname))
View(refOR_mtxIR_DAS28remission2)
# which(row(refOR)==col(refOR))
# refOR[row(refOR)==col(refOR)]
index_diagonal2<-data.frame(index1=1:nrow(refOR_mtxIR_DAS28remission2),index2=1:nrow(refOR_mtxIR_DAS28remission2))
options(warn=-1)
for(i in 1:nrow(refOR_mtxIR_DAS28remission2)){
  refOR_mtxIR_DAS28remission2[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_mtxIR_DAS28remission2[row(refOR_mtxIR_DAS28remission2)==col(refOR_mtxIR_DAS28remission2)])[i],mtxIR_subl_DAS28remission2)
}
options(warn=0)
View(refOR_mtxIR_DAS28remission2)
write.csv(refOR_mtxIR_DAS28remission2,"DAS28 Remission response pairwise comparison MTX-IR-sub-network_AUG2026.csv", row.names = T)

#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "DAS28 Remission-sub-network")
writeData(wb, "DAS28 Remission-sub-network", refOR_mtxIR_DAS28remission2, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)

###############################################################################
###############################################################################
###############################################################################
#####Radiographic Progression
#RadiograhicP_dataset_mtxIR<-final_dmardDataset[,c(1:51,862:933)]
RadiograhicP_dataset_mtxIR2<-final_mtx_irDataset_actual%>% 
  select(matches("^Radiographic.Progression\\..*")) %>%
  select(-contains("Reported"))%>%
  select(-matches("\\.1$"))
colnames(RadiograhicP_dataset_mtxIR2)

RadiograhicP_dataset_mtxIR2[RadiograhicP_dataset_mtxIR2=="NA"]<-NA
View(RadiograhicP_dataset_mtxIR2)
actuals<-RadiograhicP_dataset_mtxIR2%>%
  select(matches("^Radiographic.Progression\\.12\\.months..actual*"))
View(actuals)
colnames(actuals)
library(tibble)# library to add column
actuals<-add_column(actuals, "Radiographic.Progression.12.months..actual.value.Mean_edit" = actuals$Radiographic.Progression.12.months..actual.value.Mean, .after = "Radiographic.Progression.12.months..actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
index<-which(is.na(actuals$Radiographic.Progression.12.months..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){#evaluate when entire dataset is NA
  actuals<-add_column(actuals, "Radiographic.Progression.12.months..actual.value.SD_edit" = actuals$Radiographic.Progression.12.months..actual.value.SD, .after = "Radiographic.Progression.12.months..actual.value.SD")
  changeFbaseline<-RadiograhicP_dataset_mtxIR2%>%
    select(matches("^Radiographic.Progression\\.12.months\\..change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.12.months..change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.Mean , .after = "Radiographic.Progression.12.months..change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.12.months..change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.SD, .after = "Radiographic.Progression.12.months..change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Radiographic.Progression.12.months..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Radiographic.Progression.12.months..actual.value.Mean_edit"]<-sub_actuals[index2,"Radiographic.Progression.12.months..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  actuals<-add_column(actuals, "Radiographic.Progression.12.months..actual.value.SD_edit" = actuals$Radiographic.Progression.12.months..actual.value.SD, .after = "Radiographic.Progression.12.months..actual.value.SD")
  
  #for missing variability
  index3<-which(is.na(actuals$DAS.12.months..actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Radiographic.Progression.12.months..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Radiographic.Progression.12.months..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Radiographic.Progression.12.months..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$DAS.12.months..actual.value.Upper.CI))
      temp_sd[index7,"Radiographic.Progression.12.months..actual.value.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.12.months..actual.value.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.12.months..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.12.months..actual.value.N"])
      sub_actuals_sd[index6,] <-temp_sd
      actuals[index3,]<-sub_actuals_sd
    } 
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.12.months..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.12.months..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.N"][which.max(temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.N"][which.min(temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.12.months..actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.12.months..actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.12.months..actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.12.months..actual.value.Mean"],
      #                                                        temp_sd2[index9,"DAS.12.months..actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  
  ####12 months Change from baseline
  changeFbaseline<-RadiograhicP_dataset_mtxIR2%>%
    select(matches("^Radiographic.Progression\\.12.months\\..change.from.baseline*"))
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.12.months..change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.Mean , .after = "Radiographic.Progression.12.months..change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Radiographic.Progression.12.months..change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Radiographic.Progression.12.months..change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Radiographic.Progression.12.months..change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.12.months..change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.SD, .after = "Radiographic.Progression.12.months..change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Radiographic.Progression.12.months..change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Radiographic.Progression.12.months..change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Radiographic.Progression.12.months..change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Radiographic.Progression.12.months..change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.12.months..change.from.baseline.Upper.CI))
      temp_sd[index7,"Radiographic.Progression.12.months..change.from.baseline.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.12.months..change.from.baseline.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.12.months..change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.12.months..change.from.baseline.N"])
      sub_changeFbaseline_sd[index6,] <-temp_sd
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.12.months..change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.12.months..change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.N"][which.max(temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.N"][which.min(temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.12.months..change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.12.months..change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.12.months..change.from.baseline.p.value"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.Mean"],
      #                                                         temp_sd2[index9,"DAS.12.months..change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}


##Replacing processed Radiographic Progression outcome-12months
RadiographicProgression_12months_continuous_outcome<-cbind(actuals,changeFbaseline)

####3 months Actuals
dataset_RadiographicProgression3months<-RadiograhicP_dataset_mtxIR2%>%
  select(matches("^Radiographic.Progression\\.3\\.months*"))

dataset_RadiographicProgression3months[dataset_RadiographicProgression3months=="NA"]<-NA
View(dataset_RadiographicProgression3months)
actuals<-dataset_RadiographicProgression3months%>%
  select(matches("^Radiographic.Progression\\.3\\.months...actual*"))
View(actuals)
colnames(actuals)
#for missing means
actuals<-add_column(actuals, "Radiographic.Progression.3.months...actual.value.Mean_edit" = actuals$Radiographic.Progression.3.months...actual.value.Mean, .after = "Radiographic.Progression.3.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Radiographic.Progression.3.months...actual.value.Mean))

if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Radiographic.Progression.3.months...actual.value.SD_edit" = actuals$Radiographic.Progression.3.months...actual.value.SD, .after = "Radiographic.Progression.3.months...actual.value.SD")
  changeFbaseline<-dataset_RadiographicProgression3months%>%
    select(matches("^Radiographic.Progression\\.3.months\\...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.3.months...change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.Mean , .after = "Radiographic.Progression.3.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.3.months...change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.SD, .after = "Radiographic.Progression.3.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Radiographic.Progression.3.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Radiographic.Progression.3.months...actual.value.Mean_edit"]<-sub_actuals[index2,"Radiographic.Progression.3.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Radiographic.Progression.3.months...actual.value.SD_edit" = actuals$Radiographic.Progression.3.months...actual.value.SD, .after = "Radiographic.Progression.3.months...actual.value.SD")
  index3<-which(is.na(actuals$Radiographic.Progression.3.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Radiographic.Progression.3.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Radiographic.Progression.3.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Radiographic.Progression.3.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.3.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.3.months...actual.value.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.3.months...actual.value.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.3.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.3.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.3.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.3.months...actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.N"][which.max(temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.N"][which.min(temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.3.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####3 months Change from baseline
  changeFbaseline<-dataset_RadiographicProgression3months%>%
    select(matches("^Radiographic.Progression\\.3.months\\...change.from.baseline*"))
  changeFbaseline[] <- lapply(changeFbaseline, function(x) as.numeric(as.character(x)))#converting character numeric to numeric
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.3.months...change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.Mean , .after = "Radiographic.Progression.3.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Radiographic.Progression.3.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Radiographic.Progression.3.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Radiographic.Progression.3.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.3.months...change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.SD, .after = "Radiographic.Progression.3.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Radiographic.Progression.3.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Radiographic.Progression.3.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Radiographic.Progression.3.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Radiographic.Progression.3.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.3.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.3.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.3.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.3.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.3.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.3.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.3.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.N"][which.max(temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.N"][which.min(temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.3.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}




##Replacing processed Radiographic Progression outcome-3months
RadiographicProgression_3months_continuous_outcome<-cbind(actuals,changeFbaseline)


#################Radiographic Progression 6 months
####6 months Actuals
dataset_RadiographicProgression6months<-RadiograhicP_dataset_mtxIR2%>%
  select(matches("^Radiographic.Progression\\.6\\.months*"))
dataset_RadiographicProgression6months[dataset_RadiographicProgression6months=="NA"]<-NA
View(dataset_RadiographicProgression6months)
actuals<-dataset_RadiographicProgression6months%>%
  select(matches("^Radiographic.Progression\\.6\\.months...actual*"))
View(actuals)
colnames(actuals)

#for missing means
actuals<-add_column(actuals, "Radiographic.Progression.6.months...actual.value.Mean_edit" = actuals$Radiographic.Progression.6.months...actual.value.Mean, .after = "Radiographic.Progression.6.months...actual.value.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Radiographic.Progression.6.months...actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Radiographic.Progression.6.months...actual.value.SD_edit" = actuals$Radiographic.Progression.6.months...actual.value.SD, .after = "Radiographic.Progression.6.months...actual.value.SD")
  changeFbaseline<-dataset_RadiographicProgression6months%>%
    select(matches("^Radiographic.Progression\\.6\\.months...change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.6.months...change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.Mean , .after = "Radiographic.Progression.6.months...change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.6.months...change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.SD, .after = "Radiographic.Progression.6.months...change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Radiographic.Progression.6.months...actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Radiographic.Progression.6.months...actual.value.Mean_edit"]<-sub_actuals[index2,"Radiographic.Progression.6.months...actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Radiographic.Progression.6.months...actual.value.SD_edit" = actuals$Radiographic.Progression.6.months...actual.value.SD, .after = "Radiographic.Progression.6.months...actual.value.SD")
  index3<-which(is.na(actuals$Radiographic.Progression.6.months...actual.value.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Radiographic.Progression.6.months...actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Radiographic.Progression.6.months...actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Radiographic.Progression.6.months...actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.6.months...actual.value.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.6.months...actual.value.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.6.months...actual.value.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.6.months...actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.6.months...actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
      
    } 
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.6.months...actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.6.months...actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.N"][which.max(temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.N"][which.min(temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.6.months...actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  
  ####6 months Change from baseline
  changeFbaseline<-dataset_RadiographicProgression6months%>%
    select(matches("^Radiographic.Progression\\.6\\.months...change.from.baseline*"))
  colnames(changeFbaseline)
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.6.months...change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.Mean , .after = "Radiographic.Progression.6.months...change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Radiographic.Progression.6.months...change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Radiographic.Progression.6.months...change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Radiographic.Progression.6.months...change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.6.months...change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.SD, .after = "Radiographic.Progression.6.months...change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Radiographic.Progression.6.months...change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Radiographic.Progression.6.months...change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Radiographic.Progression.6.months...change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Radiographic.Progression.6.months...change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.6.months...change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.6.months...change.from.baseline.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.6.months...change.from.baseline.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.6.months...change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.6.months...change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
      
      
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.6.months...change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.6.months...change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.N"][which.max(temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.N"][which.min(temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.6.months...change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
  
  
}

##Replacing processed Radiographic Progression outcome-6months
RadiographicProgression_6months_continuous_outcome<-cbind(actuals,changeFbaseline)


options(warn = defaultW)
################Radiographic Progression baseline
####Radiographic Progression baseline Actuals
dataset_RadiographicProgressionBaseline<-RadiograhicP_dataset_mtxIR2%>%
  select(matches("^Radiographic.Progression\\.Baseline*"))
dataset_RadiographicProgressionBaseline[dataset_RadiographicProgressionBaseline=="NA"]<-NA
View(dataset_DASBaseline)
actuals<-dataset_RadiographicProgressionBaseline
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "Radiographic.Progression.Baseline.Mean_edit" = actuals$Radiographic.Progression.Baseline.Mean, .after = "Radiographic.Progression.Baseline.Mean")
actuals[] <- lapply(actuals, function(x) as.numeric(as.character(x)))#converting character numeric to numeric

index<-which(is.na(actuals$Radiographic.Progression.Baseline.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Radiographic.Progression.Baseline.SD_edit" = actuals$Radiographic.Progression.Baseline.SD, .after = "Radiographic.Progression.Baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Radiographic.Progression.Baseline.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Radiographic.Progression.Baseline.Mean_edit"]<-sub_actuals[index2,"Radiographic.Progression.Baseline.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Radiographic.Progression.Baseline.SD_edit" = actuals$Radiographic.Progression.Baseline.SD, .after = "Radiographic.Progression.Baseline.SD")
  index3<-which(is.na(actuals$Radiographic.Progression.Baseline.SD))
  #View(actuals[index3,])
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Radiographic.Progression.Baseline.SD_edit"]<-(sub_actuals_sd[index4,"DAS.Baseline.SE"])*sqrt(sub_actuals_sd[index4,"Radiographic.Progression.Baseline.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.Baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.Baseline.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.Baseline.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.Baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.Baseline.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.Baseline.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.Baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.Baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.Baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.Baseline.N"][which.max(temp_sd2[index9,"Radiographic.Progression.Baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.Baseline.N"][which.min(temp_sd2[index9,"Radiographic.Progression.Baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.Baseline.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.Baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
}


options(warn = defaultW)
##Replacing processed DAS outcome-Baseline
RadiographicProgressionBaseline_continuous_outcome<-actuals

#################RadiographicProgression End of trial
####End of trial Actuals
dataset_RadiographicProgression_Endoftrial<-RadiograhicP_dataset_mtxIR2%>%
  select(matches("^Radiographic.Progression\\.End.of.trial*"))%>%
  select(-contains("month.24"))
dataset_RadiographicProgression_Endoftrial[dataset_RadiographicProgression_Endoftrial=="NA"]<-NA
View(dataset_RadiographicProgression_Endoftrial)
actuals<-dataset_RadiographicProgression_Endoftrial%>%
  select(matches("^Radiographic.Progression\\.End.of.trial\\..actual*"))
View(actuals)
colnames(actuals)
defaultW <- getOption("warn") 
options(warn = -1) 
#for missing means
actuals<-add_column(actuals, "Radiographic.Progression.End.of.trial..actual.value.Mean_edit" = actuals$Radiographic.Progression.End.of.trial..actual.value.Mean, .after = "Radiographic.Progression.End.of.trial..actual.value.Mean")
index<-which(is.na(actuals$Radiographic.Progression.End.of.trial..actual.value.Mean))
if(length(as.vector(which(colSums(is.na(actuals)) == nrow(actuals))))==ncol(actuals)){
  actuals<-add_column(actuals, "Radiographic.Progression.End.of.trial..actual.value.SD_edit" = actuals$Radiographic.Progression.End.of.trial..actual.value.SD, .after = "Radiographic.Progression.End.of.trial..actual.value.SD")
  changeFbaseline<-dataset_RadiographicProgression_Endoftrial%>%
    select(matches("^Radiographic.Progression\\.End.of.trial\\.change.from.baseline*"))
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.Mean , .after = "Radiographic.Progression.End.of.trial.change.from.baseline.Mean")
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.SD, .after = "Radiographic.Progression.End.of.trial.change.from.baseline.SD")
  
}else{
  if(sum(sapply(actuals[index,], function(x)all(is.na(x))))!=8){
    sub_actuals<-actuals[index,]
    index2<-which(!is.na(sub_actuals[,"Radiographic.Progression.End.of.trial..actual.value.Median"]))
    if(length(index2)!=0){
      sub_actuals[index2,"Radiographic.Progression.End.of.trial..actual.value.Mean_edit"]<-sub_actuals[index2,"Radiographic.Progression.End.of.trial..actual.value.Median"]
      actuals[index,]<-sub_actuals
    }
  }
  
  
  #for missing variability
  actuals<-add_column(actuals, "Radiographic.Progression.End.of.trial..actual.value.SD_edit" = actuals$Radiographic.Progression.End.of.trial..actual.value.SD, .after = "Radiographic.Progression.End.of.trial..actual.value.SD")
  index3<-which(is.na(actuals$Radiographic.Progression.End.of.trial..actual.value.SD))
  if(sum(sapply(actuals[index3,], function(x)all(is.na(x))))!=8){
    sub_actuals_sd<-actuals[index3,]
    index4<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.SE"]))
    index5<-which(!is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.Upper.CI"]))
    if(length(index4)!=0){
      sub_actuals_sd[index4,"Radiographic.Progression.End.of.trial..actual.value.SD_edit"]<-(sub_actuals_sd[index4,"Radiographic.Progression.End.of.trial..actual.value.SE"])*sqrt(sub_actuals_sd[index4,"Radiographic.Progression.End.of.trial..actual.value.N"])
      actuals[index3,]<-sub_actuals_sd
    }
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.SD"])))!=0)){
      index6<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.SD"]))
      temp_sd<-sub_actuals_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.End.of.trial..actual.value.Upper.CI))
      
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.End.of.trial..actual.value.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.End.of.trial..actual.value.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.End.of.trial..actual.value.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.End.of.trial..actual.value.N"])
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }else{
        sub_actuals_sd[index6,] <-temp_sd
        actuals[index3,]<-sub_actuals_sd
      }
    } 
    if((length(which(is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.SD"])))!=0)){
      index8<-which(is.na(sub_actuals_sd[,"Radiographic.Progression.End.of.trial..actual.value.SD"]))
      temp_sd2<-sub_actuals_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.End.of.trial..actual.value.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.N"][which.max(temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.N"][which.min(temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.End.of.trial..actual.value.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...actual.value.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...actual.value.p.value"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.Mean"],
      #                                                         temp_sd2[index9,"DAS.3.months...actual.value.N"])
      sub_actuals_sd[index8,] <-temp_sd2
      actuals[index3,]<-sub_actuals_sd
      
    }
  }
  options(warn = defaultW)
  
  
  ####End of trial Change from baseline
  changeFbaseline<-dataset_RadiographicProgression_Endoftrial%>%
    select(matches("^Radiographic.Progression\\.End.of.trial\\.change.from.baseline*"))
  colnames(changeFbaseline)
  defaultW <- getOption("warn") 
  options(warn = -1) 
  #for missing means
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit" =changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.Mean , .after = "Radiographic.Progression.End.of.trial.change.from.baseline.Mean")
  index<-which(is.na(changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.Mean))
  if(sum(sapply(changeFbaseline[index,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline<-changeFbaseline[index,]
    index2<-which(!is.na(changeFbaseline[,"Radiographic.Progression.End.of.trial.change.from.baseline.Median"]))
    if(length(index2)!=0){
      sub_changeFbaseline[index2,"Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit"]<-sub_changeFbaseline[index2,"Radiographic.Progression.End.of.trial.change.from.baseline.Median"]
      changeFbaseline[index,]<-sub_changeFbaseline
    }
  }
  #for missing variability
  changeFbaseline<-add_column(changeFbaseline, "Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit" = changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.SD, .after = "Radiographic.Progression.End.of.trial.change.from.baseline.SD")
  index3<-which(is.na(changeFbaseline$Radiographic.Progression.End.of.trial.change.from.baseline.SD))
  if(sum(sapply(changeFbaseline[index3,], function(x)all(is.na(x))))!=8){
    sub_changeFbaseline_sd<-changeFbaseline[index3,]
    index4<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.SE"]))
    index5<-which(!is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.Lower.CI"]))
    if(length(index4)!=0){
      sub_changeFbaseline_sd[index4,"Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit"]<-(sub_changeFbaseline_sd[index4,"Radiographic.Progression.End.of.trial.change.from.baseline.SE"])*sqrt(sub_changeFbaseline_sd[index4,"Radiographic.Progression.End.of.trial.change.from.baseline.N"])
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
    }
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.SD"])))!=0)){
      index6<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.SD"]))
      temp_sd<-sub_changeFbaseline_sd[index6,]
      index7<-which(!is.na(temp_sd$Radiographic.Progression.End.of.trial.change.from.baseline.Upper.CI))
      if(length(index7)!=0){
        temp_sd[index7,"Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit"]<-((temp_sd[index7,"Radiographic.Progression.End.of.trial.change.from.baseline.Upper.CI"]-temp_sd[index7,"Radiographic.Progression.End.of.trial.change.from.baseline.Lower.CI"])/3.92)*sqrt(temp_sd[index7,"Radiographic.Progression.End.of.trial.change.from.baseline.N"])
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }else{
        sub_changeFbaseline_sd[index6,] <-temp_sd
        changeFbaseline[index3,]<-sub_changeFbaseline_sd
      }
    } 
    if((length(which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.SD"])))!=0)){
      index8<-which(is.na(sub_changeFbaseline_sd[,"Radiographic.Progression.End.of.trial.change.from.baseline.SD"]))
      temp_sd2<-sub_changeFbaseline_sd[index8,]
      index9<-which(!is.na(temp_sd2$Radiographic.Progression.End.of.trial.change.from.baseline.p.value))
      meanhigh<-max(temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit"])
      meanlow<-min(temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit"])
      Nrx<-temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.N"][which.max(temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit"])]
      Npbo<-temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.N"][which.min(temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.Mean_edit"])]
      p.val<-na.omit(unique(temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.p.value"]))
      temp_sd2[index9,"Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit"]<- exactp.sd(p.val, Nrx, Npbo, meanhigh, meanlow)[[2]]
      # temp_sd2[index9,"DAS.3.months...change.from.baseline.SD"]<-p2sd(temp_sd2[index9,"DAS.3.months...change.from.baseline.p.value"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.Mean"],
      #                                                                 temp_sd2[index9,"DAS.3.months...change.from.baseline.N"])
      sub_changeFbaseline_sd[index8,] <-temp_sd2
      changeFbaseline[index3,]<-sub_changeFbaseline_sd
      
    }
  }
}

options(warn = defaultW)
##Replacing processed Radiographic Progression outcome-EndOfTrial
RadiographicProgression_End_of_trial_continuous_outcome<-cbind(actuals,changeFbaseline)


####Final Continuous outcome dataset-Radiographic progression
final_mtxIRDataset_radiographicProgression<-cbind(final_mtx_irDataset_actual[,c(1:51)],
                                                  RadiographicProgression_12months_continuous_outcome,
                                                  RadiographicProgression_3months_continuous_outcome,
                                                  RadiographicProgression_6months_continuous_outcome,
                                                  RadiographicProgressionBaseline_continuous_outcome,
                                                  RadiographicProgression_End_of_trial_continuous_outcome)
final_mtxIRDataset_radiographicProgression[,c(51:141)] <- lapply(final_mtxIRDataset_radiographicProgression[,c(51:141)], function(x) as.numeric(as.character(x)))#converting character numeric to numeric

View(final_mtxIRDataset_radiographicProgression)

##fill time points of no interest with NA
final_mtxIRDataset_radiographicProgression2 <- final_mtxIRDataset_radiographicProgression

# End of trial columns
eot_cols <- grep(
  "^Radiographic.Progression\\.End\\.of\\.trial",
  names(final_mtxIRDataset_radiographicProgression2),
  value = TRUE
)

# 12-month columns
m12_cols <- grep(
  "^Radiographic.Progression\\.12\\.months",
  names(final_mtxIRDataset_radiographicProgression2),
  value = TRUE
)

# 6-month columns
m6_cols <- grep(
  "^Radiographic.Progression\\.6\\.months",
  names(final_mtxIRDataset_radiographicProgression2),
  value = TRUE
)

# 3-month columns
m3_cols <- grep(
  "^Radiographic.Progression\\.3\\.months",
  names(final_mtxIRDataset_radiographicProgression2),
  value = TRUE
)

# Apply same rules as dichotomous outcome
final_mtxIRDataset_radiographicProgression2[
  final_mtxIRDataset_radiographicProgression2$`Trial.duration..in.months` %in%
    c(3,4,6,9,11,12,13,18),
  eot_cols
] <- NA

final_mtxIRDataset_radiographicProgression2[
  final_mtxIRDataset_radiographicProgression2$`Trial.duration..in.months` %in%
    c(3,4,23,24,63),
  m12_cols
] <- NA

final_mtxIRDataset_radiographicProgression2[
  final_mtxIRDataset_radiographicProgression2$`Trial.duration..in.months` %in%
    c(3,4,9,11,12,13,18,23,24,63),
  m6_cols
] <- NA

final_mtxIRDataset_radiographicProgression2[
  final_mtxIRDataset_radiographicProgression2$`Trial.duration..in.months` %in%
    c(6,9,11,12,13,18,23,24,63),
  m3_cols
] <- NA


#############
###getting the radio actual dataset
radiog_actual<-final_mtxIRDataset_radiographicProgression2

#Conaghan 2013
radiog_actual[radiog_actual$Study.Identifier=="Conaghan 2013","Radiographic.Progression.Baseline.Median"]<-NA

#Kremer 2011
radiog_actual[radiog_actual$Study.Identifier=="Kremer 2011","Radiographic.Progression.Baseline.SD_edit"]<-c(190.5,171.4,178.7)/4

#Machado 2014
radiog_actual[radiog_actual$Study.Identifier=="Machado 2014","Radiographic.Progression.Baseline.SD_edit"]<-c(4.06*sqrt(247),5.35*sqrt(119))



#radiog_actual<-final_mtxIRDataset_radiographicProgression[final_mtxIRDataset_radiographicProgression$Study_no%in%c(8,15,17,21,40,49,58,60,65,66,72,75),]
View(radiog_actual[radiog_actual$Study_no==75,])
tttt<-radiog_actual[radiog_actual$Study_no==98,c(52:141)]
tttt_studyInfo<-radiog_actual[radiog_actual$Study_no==98,c(1:3)]
tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
tttt[colSums(!is.na(tttt)) > 0]
as.integer(as.vector(which(colSums(!is.na(tttt)) > 0)))
identical(as.integer(c(81:83,85,86)),as.integer(as.vector(which(colSums(!is.na(tttt)) > 0))))

getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
getting_actuals2<-getting_actuals[, grepl("actual.value",names(getting_actuals))]
View(getting_actuals2)

temp_data3<-list()
error_check_naive<-list()
empty_entries3 <- list()

for(i in 1:length(unique(radiog_actual$Study_no))){
  tttt<-radiog_actual[radiog_actual$Study_no==i,c(52:141)]
  tttt_studyInfo<-radiog_actual[radiog_actual$Study_no==i,c(1:3)]
  tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
  
  if(sum(which(colSums(!is.na(tttt)) > 0))==0){
    empty_entries3[[i]]<-radiog_actual[radiog_actual$Study_no==i,c(1:3,52:141)]
  }else if(identical(as.integer(c(81:83,85,86)),as.integer(as.vector(which(colSums(!is.na(tttt)) > 0))))){
    ##Only baseline value
    eee<-cbind(tttt_studyInfo,tttt[colSums(!is.na(tttt)) > 0]%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(tttt[colSums(!is.na(tttt)) > 0])[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data3[[i]]<- eee
  }else if(any(c(5:9,25:29,45:49,78:79, 85)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##actual only but need to choose month
    getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
    getting_actuals2<-getting_actuals[, grepl("actual.value",names(getting_actuals))]
    if(ncol(getting_actuals2)==5){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[3])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==6){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[4])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==7){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[4])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==10){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==11){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==15){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[13])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==16){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[15])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==18){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[17])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(ncol(getting_actuals2)==20){
      dd<-getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))
      # get the last 3 columns
      last3 <- names(dd)[(ncol(dd)-2):ncol(dd)]
      # order them by suffix
      ordered_last3 <- last3[order(match(
        str_extract(last3, "(SD_edit|N|Mean_edit)$"),
        c("SD_edit","N","Mean_edit")
      ))]
      eee<-cbind(tttt_studyInfo,dd %>%select(all_of(ordered_last3)),timepoint=colnames(getting_actuals2)[18])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data3[[i]]<- eee
    }else if(!ncol(getting_actuals2)%in%c(5,6,10,11,15,16,18,20)){
      error_check_naive[[i]]<-radiog_actual[radiog_actual$Study_no==i,]
    }
  }else if(any(c(11,31,51,71, 81)%in%as.vector(which(colSums(!is.na(tttt)) > 0)))){
    ##Baseline value and change in baseline value
    n_col<-ncol(tttt[colSums(!is.na(tttt)) > 0])
    #baseline<-tttt[colSums(!is.na(tttt)) > 0][,sort(n_col:((n_col-5)+1))]
    baseline<-tttt[colSums(!is.na(tttt)) > 0][,grep("^Radiographic.Progression\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE)]
    #change_mean_values<-tttt[colSums(!is.na(tttt)) > 0][,-sort(n_col:((n_col-5)+1))]
    change_mean_values<-tttt[colSums(!is.na(tttt)) > 0]%>%
      select(!grep("^Radiographic.Progression\\.Baseline",names(tttt[colSums(!is.na(tttt)) > 0]),value = TRUE))
    change_mean_values2<-as.data.frame(change_mean_values[ , grepl( "change.from.baseline.Mean_edit" , names(change_mean_values) ) ])
    actual_mean<-baseline[,2]+change_mean_values2[,ncol(change_mean_values2)]
    baseline[,"Radiographic.Progression.Baseline.Mean_edit"]<-actual_mean
    replaced_baseline_mean_with_actual<-baseline
    eee<-cbind(tttt_studyInfo,replaced_baseline_mean_with_actual%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                 select(
                   matches("\\.SD_edit$"),
                   matches("\\.N$"),
                   matches("\\.Mean_edit$")
                 ),timepoint=colnames(replaced_baseline_mean_with_actual)[3])
    colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
    temp_data3[[i]]<- eee
  }
  
}

View(do.call(rbind,empty_entries3))
View(do.call(rbind,temp_data3))
View(do.call(rbind,error_check_naive))


RadiographicP_mtxIR_actual_dataset_analysis<-do.call(rbind,temp_data3)
RadiographicP_mtxIR_actual_dataset_analysis1<-RadiographicP_mtxIR_actual_dataset_analysis
row.names(RadiographicP_mtxIR_actual_dataset_analysis1)<-NULL

##Drop furst study due to high ROB
RadiographicP_mtxIR_actual_dataset_analysis1<- RadiographicP_mtxIR_actual_dataset_analysis1 %>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012"))

RadiographicP_mtxIR_actual_dataset_analysis1$Study_no<-renumber(RadiographicP_mtxIR_actual_dataset_analysis1,"Study_no")                                                                                                                                         
row.names(RadiographicP_mtxIR_actual_dataset_analysis1)<-NULL

#####Formatting the month column
RadiographicP_mtxIR_actual_dataset_analysis1$Timepoint <- gsub(
  "Radiographic.Progression\\.|\\.months\\.{2,3}actual\\.value\\.(N|Mean_edit|Median|SD|Lower\\.CI)|\\.(Mean_edit|Median)",
  "",
  RadiographicP_mtxIR_actual_dataset_analysis1$Timepoint
)



RadiographicP_mtxIR_actual_dataset_analysis1 <- RadiographicP_mtxIR_actual_dataset_analysis1 %>%
  dplyr::rename("treatment" = Intervention)



library(dplyr)
for(i in 1:nrow(RadiographicP_mtxIR_actual_dataset_analysis1)){
  ppp<- master.medtable%>% filter(RadiographicP_mtxIR_actual_dataset_analysis1$treatment[i] ==master.medtable$medname )
  RadiographicP_mtxIR_actual_dataset_analysis1$treatment[i]<-ppp[1,1]
}

##Creating initial legend
RadiographicP.actual.mtxIR.medlegend <- new.medtable(RadiographicP_mtxIR_actual_dataset_analysis1,'treatment')
RadiographicP.actual.mtxIR.medlegend

##NMA plot reason for narrowing of network
nma.plot(RadiographicP_mtxIR_actual_dataset_analysis1[,1],RadiographicP_mtxIR_actual_dataset_analysis1[,3],data = RadiographicP_mtxIR_actual_dataset_analysis1, trtname = as.vector(RadiographicP.actual.mtxIR.medlegend$medname))

################################################################################################################
#trying the narrowing a network function
#############################################################################################
RadiographicP_mtxIR_actual_dataset_analysis1

RadiographicP_mtxIR_actual_dataset_analysis1<- RadiographicP_mtxIR_actual_dataset_analysis1%>%
  dplyr::rename("study" = Study_no,"study_name" = Study.Identifier)
#Narrow 1
mtxIR_RadiographicP_actual_narrow<-narrow.all(RadiographicP_mtxIR_actual_dataset_analysis1)
mtxIR_RadiographicP_actual_narrow

#sub network plot 1
mtxIR.RadiographicP_actual.d<-mtxIR_RadiographicP_actual_narrow[[1]]
length(unique(mtxIR.RadiographicP_actual.d$treatment))
mtxIR.RadiographicP_actual_subl<-new.medtable(mtxIR.RadiographicP_actual.d,'treatment')

write.csv(mtxIR.RadiographicP_actual.d,"Radiographic Progression MTX-IR dataset_AUG2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("Radiographic Progression-Actual MTX-IR Network plot with named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.RadiographicP_actual.d[,1],mtxIR.RadiographicP_actual.d[,3],data = mtxIR.RadiographicP_actual.d, trtname = as.vector(mtxIR.RadiographicP_actual_subl$medname),title = "Radiographic Progression mtxIR Network plot-Actual")

# Close the pdf file
dev.off() 
mtxIR.RadiographicP_actual.d1<-mtxIR.RadiographicP_actual.d
mtxIR.RadiographicP_actual.d1$treatment<-renumber(mtxIR.RadiographicP_actual.d1,"treatment")
# 2. Open jpeg file
jpeg("Radiographic Progression-Actual MTX-IR Network plot without named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.RadiographicP_actual.d1[,1],mtxIR.RadiographicP_actual.d1[,3],data = mtxIR.RadiographicP_actual.d1,title = "Radiographic Progression mtxIR Network plot-Actual")
# Close the pdf file
dev.off() 
##############################################################################
#Performing sub network plot 1
library(coda)
library(gemtc)
library(netmeta)
library(tidyr)
View(mtxIR.RadiographicP_actual.d)
mtxIR.RadiographicP_actual_subl

network_data_RadiographicP_actual_mtxIR_1<-mtxIR.RadiographicP_actual.d%>%
  select(study, study_name, treatment, SD, sampleSize, Mean, Timepoint)


network_data_RadiographicP_actual_mtxIR_1$treatment<-renumber(network_data_RadiographicP_actual_mtxIR_1,"treatment")
network_data_RadiographicP_actual_mtxIR_1$study<-renumber(network_data_RadiographicP_actual_mtxIR_1,"study")
network_data_RadiographicP_actual_mtxIR_1<-network_data_RadiographicP_actual_mtxIR_1[order(network_data_RadiographicP_actual_mtxIR_1$study),]
row.names(network_data_RadiographicP_actual_mtxIR_1)<-NULL
nma.plot(network_data_RadiographicP_actual_mtxIR_1$study,network_data_RadiographicP_actual_mtxIR_1$treatment,data = network_data_RadiographicP_actual_mtxIR_1, trtname = as.vector(mtxIR.RadiographicP_actual_subl$medname),title = "Radiographic Progression mtxIR Network plot-Actual")

network_data_RadiographicP_actual_mtxIR_1<- network_data_RadiographicP_actual_mtxIR_1%>%
  dplyr::rename("std.dev"=SD, "mean"= Mean)


p_RadiographicP_actual_MD<-pairwise(treat=treatment,mean=mean,sd=std.dev,data=network_data_RadiographicP_actual_mtxIR_1,
                                    studlab = study,n=sampleSize,sm="SMD")

#colnames(p_RadiographicP_actual_MD)[c(1,4:5)]<-c("study","diff", "std.err")

p_RadiographicP_actual_MD<- p_RadiographicP_actual_MD%>%
  dplyr::rename("diff" = TE, "std.err" = seTE)

p_RadiographicP_actual_MD %>%
  dplyr::select(c(study_name,diff, std.err,treat1,treat2)) %>%
  pivot_longer(-study_name,
               names_to = c(".value"),
               names_pattern = "(..)") %>% 
  set_colnames(c("study", "diff", 
                 "std.err", "treatment"))  -> data6

data6 <- clean_and_impute_se(data6)
network_data_radiographicP_actual_mtxIR_1<-mtc.network(data.re=data6)
summary(network_data_radiographicP_actual_mtxIR_1)
View(network_data_radiographicP_actual_mtxIR_1$data.re)

##RUN the model
set.seed(5811156)
network_data_radiographicP_actual_mtxIR_1.model <- mtc.model(network_data_radiographicP_actual_mtxIR_1, likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_radiographicP_actual_mtxIR_1.model$inits <- mapply(c, network_data_radiographicP_actual_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_radiographicP_actual<-mtc.run(network_data_radiographicP_actual_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_radiographicP_actual)

View(as.data.frame(relative.effect.table(results_mtxIR_radiographicP_actual)))
View(as.data.frame(relative.effect.table.2DP2(results_mtxIR_radiographicP_actual)))
refME_mtxIR_radiographicP_actual<-as.data.frame(relative.effect.table.2DP2(results_mtxIR_radiographicP_actual, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refME_mtxIR_radiographicP_actual)<-mapvalues(row.names(refME_mtxIR_radiographicP_actual),
                                                       from = c(1:(nrow(mtxIR.RadiographicP_actual_subl))),
                                                       to = as.vector(mtxIR.RadiographicP_actual_subl$medname))
colnames(refME_mtxIR_radiographicP_actual)<-mapvalues(colnames(refME_mtxIR_radiographicP_actual),
                                                      from = c(1:(nrow(mtxIR.RadiographicP_actual_subl))),
                                                      to = as.vector(mtxIR.RadiographicP_actual_subl$medname))
View(refME_mtxIR_radiographicP_actual)

index_diagonal5<-data.frame(index1=1:nrow(refME_mtxIR_radiographicP_actual),index2=1:nrow(refME_mtxIR_radiographicP_actual))
options(warn=-1)
for(i in 1:nrow(refME_mtxIR_radiographicP_actual)){
  refME_mtxIR_radiographicP_actual[index_diagonal5$index1[i],index_diagonal5$index2[i]]<-selected_treatment_name(as.numeric(refME_mtxIR_radiographicP_actual[row(refME_mtxIR_radiographicP_actual)==col(refME_mtxIR_radiographicP_actual)])[i],mtxIR.RadiographicP_actual_subl)
}
options(warn=0)
View(refME_mtxIR_radiographicP_actual)


write.csv(refME_mtxIR_radiographicP_actual,"Radiographic Progression-Actual MTX-IR pairwise comparisons_AUG2026.csv", row.names = T)
#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "Radiographic Progression-Actual")
writeData(wb, "Radiographic Progression-Actual", refME_mtxIR_radiographicP_actual, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)

##GRADE Assessment
nma_estimate_mtxIR_radiographicP_actual_dataframe<-nma.estimate.table.MD(results_mtxIR_radiographicP_actual)#use results
View(nma_estimate_mtxIR_radiographicP_actual_dataframe)

#Direct estimate and Indirect estimate
##Node splitting, not working
network_data_radiographicP_actual_mtxIR_1.ns<-mtc.nodesplit(network_data_radiographicP_actual_mtxIR_1, thin=50,likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 

names(network_data_radiographicP_actual_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_radiographicP_actual_mtxIR_1.ns)

# # Overall summary and plot
# summary.network_data_radiographicP_actual_mtxIR_1.ns <- summary(network_data_radiographicP_actual_mtxIR_1.ns)
# 
# #Direct effects
# direct.data<-summary.network_data_radiographicP_actual_mtxIR_1.ns$dir.effect
# 
# #Indirect effects
# indirect.data<-summary.network_data_radiographicP_actual_mtxIR_1.ns$ind.effect
# 
# #Network estimate
# network.data<-summary.network_data_radiographicP_actual_mtxIR_1.ns$cons.effect
# 
# #P-values
# P_values_nodesplit_radiographicP_actualmtxIR<-summary.network_data_radiographicP_actual_mtxIR_1.ns$p.value
# P_values_nodesplit_radiographicP_actualmtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_radiographicP_actualmtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
# P_values_nodesplit_radiographicP_actualmtxIR2<-P_values_nodesplit_radiographicP_actualmtxIR
# P_values_nodesplit_radiographicP_actualmtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_radiographicP_actualmtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
# P_values_nodesplit_radiographicP_actualmtxIR2$comb<-paste0(P_values_nodesplit_radiographicP_actualmtxIR2$t1,":",P_values_nodesplit_radiographicP_actualmtxIR2$t2)
# 
# 
# direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(direct.data[,c(3,4,5)],2),round(indirect.data[,c(3,4,5)],2))
# colnames(direct_indirect_estimate)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")
# 
# direct_indirect_estimate_list<-list()
# 
# for(i in 1:nrow(direct_indirect_estimate)){
#   if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
#     direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
#   }else{
#     temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
#                 -round(direct_indirect_estimate[i,3:8],2))
#     temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
#     colnames(temp2)<-c("treat1","treat2","MD1","low1","high1","MD2","low2","high2")
#     direct_indirect_estimate_list[[i]]<-temp2
#   }
# }
# 
# direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)

# ##Acquiring both direct and indirect estimate from gemtc
# grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_radiographicP_actual_dataframe,
#                     by=c("treat1","treat2"),all.x=T,all.y=T)
###Adding ROB data
radiographicP_actual_mtxIR<-unique(mtxIR.RadiographicP_actual.d$study_name)
data_qq_mtxIR<-mtxIR.RadiographicP_actual.d
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%radiographicP_actual_mtxIR,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<- current_assessment_for_studies_mtxIR%>%
  select(c(col_names))


list_appraisal_mtxIR3<-list()
for(i in 1:length(radiographicP_actual_mtxIR)){
  temp11<-data_qq_mtxIR[data_qq_mtxIR$study_name==radiographicP_actual_mtxIR[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==radiographicP_actual_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal_mtxIR3[[i]]<-temp13
}
radiographicP_actual_mtxIR_withROB<-do.call(rbind,list_appraisal_mtxIR3)
View(radiographicP_actual_mtxIR_withROB)
radiographicP_actual_mtxIR_withROB2<-radiographicP_actual_mtxIR_withROB%>%
  select(c("study","study_name", "treatment", "SD", "sampleSize",
           "Mean","Timepoint","OVERALL..ALL.OTHER.OUTCOMES"))
View(radiographicP_actual_mtxIR_withROB2)

###Adding Indirectness data
study_unique<-unique(radiographicP_actual_mtxIR_withROB2$study_name)
data_nn<-radiographicP_actual_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp%>%
  dplyr::rename( "Consensus"= `Final Rating`)

View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]

list_appraisal2<-list()

for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
radiographicP_actual_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal2)
radiographicP_actual_mtxIR_withROB2_indrectness<- radiographicP_actual_mtxIR_withROB2_indrectness%>%
  dplyr::rename(
    ROB = OVERALL..ALL.OTHER.OUTCOMES,
    Indirectness = Consensus
  )
radiographicP_actual_mtxIR_withROB2_indrectness$Indirectness<-tolower(radiographicP_actual_mtxIR_withROB2_indrectness$Indirectness)
View(radiographicP_actual_mtxIR_withROB2_indrectness)

radiographicP_actual_mtxIR_data_new_withROB2_indrectness2<-radiographicP_actual_mtxIR_withROB2_indrectness
radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(radiographicP_actual_mtxIR_data_new_withROB2_indrectness2,"treatment")
radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$ROB<-ifelse(radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$ROB=="unclear","some concerns",radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$ROB)
View(radiographicP_actual_mtxIR_data_new_withROB2_indrectness2)

p_radiographicP_actual_SMD2<-p_RadiographicP_actual_MD
p_radiographicP_actual_SMD2<- p_radiographicP_actual_SMD2%>%
  dplyr::rename("TE" = diff,"seTE"= std.err)

net1 <- netmeta(p_radiographicP_actual_SMD2)
dat.net<-netsplit(net1)
View(dat.net$indirect.common)
View(cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)],dat.net$compare.common[,c(2:5)]))
View(dat.net)
dat.net.contrib<-netcontrib(net1, method = "r")

direct_evidence_new<-cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)])


#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")
nma_estimate_mtxIR_radiographicP_actual_dataframe
grade_table1<-cbind(nma_estimate_mtxIR_radiographicP_actual_dataframe[,c(1,2)],
                    as.data.frame(matrix(NA,nrow(nma_estimate_mtxIR_radiographicP_actual_dataframe),6)),
                    nma_estimate_mtxIR_radiographicP_actual_dataframe[,-c(1,2)])
colnames(grade_table1)[3:8]<-c("MD/SMD1","low1","high1","MD/SMD2","low2","high2")


#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]

# ##Getting direct evidence for intervention in only multi-arm trials
# zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
# colnames(zzz)<-c("treat1", "treat2")
# grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
# inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
# grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]





grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR.RadiographicP_actual_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR.RadiographicP_actual_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))

#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))

library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
RadiographicP_actual_output_df2_list<-list()
I_square_list<-list()
radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$study<-renumber(radiographicP_actual_mtxIR_data_new_withROB2_indrectness2,"study")



for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-radiographicP_actual_mtxIR_data_new_withROB2_indrectness2[radiographicP_actual_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select(study, Mean.1, SD.1, Mean.2, SD.2 )
  temp11.1[,c("SD.1","SD.2")]<-temp11.1[,c("SD.1","SD.2")]/(sqrt(temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c("sampleSize.1" , "sampleSize.2")]))
  
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","se.t","x.c","se.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  
  
  
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    ##For one study
    n_i<-temp11[k,c("sampleSize.1" , "sampleSize.2")]
    sds<-temp11[k,c("SD.1","SD.2")]
    
    # Calculate the pooled SD for SMD
    pooled_sd <- sqrt(sum((n_i - 1) * sds^2) / sum(n_i - 1))
    model.file<-"Cont_direct_SMD_onestudy.txt"
    #model.file<-"Cont_direct_onestudy.txt"
    n.chains=3
    set.seed(66451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c,
                                   "pooled.sd.c"=pooled_sd,
                                   "pooled.sd.t"=pooled_sd
                         ),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    # params <- c("diff")
    params <- c("SMD")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(MD=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    
    data4sdPool<-temp11[c("sampleSize.1","SD.1","sampleSize.2","SD.2")]
    # Define a function to calculate the pooled standard deviation
    calculate_pooled_sd <- function(row) {
      sqrt(sum((row[c("sampleSize.1", "sampleSize.2")] - 1) * row[c("SD.1", "SD.2")]^2) / sum(row[c("sampleSize.1", "sampleSize.2")] - 1))
    }
    
    # Apply the function row-wise
    data4sdPool$pooled_sd <- apply(data4sdPool, 1, calculate_pooled_sd)
    
    #model.file<-"Cont_direct_multiplestudy_RE.txt" #Check this again  #"Cont_direct_SMD_multiplestudy_RE.txt"
    model.file<- "Cont_direct_SMD_multiplestudy_RE.txt"
    
    n.chains=3
    set.seed(77451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c, "NS"=length(temp12$NS),
                                   "pooled.sd.c"=data4sdPool[,5],#rep(plllsd,length(temp12$NS)) ,#rep(pooled_sd.c,length(temp12$NS)) ,
                                   "pooled.sd.t"=data4sdPool[,5]#rep(plllsd,length(temp12$NS))#rep(pooled_sd.t,length(temp12$NS))
                         ),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    #params <- c("diff")
    params <- c("SMD")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  
  #Frequentist approach for I^2
  nma_dataRadiographicP_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select(study, sampleSize.1, Mean.1, SD.1, sampleSize.2, Mean.2, SD.2 )
  #nma_dataDAS28_3[,c("SD.1","SD.2")]<-nma_dataDAS28_3[,c("SD.1","SD.2")]/(sqrt(nma_dataDAS28_3[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(2,5)]))
  nma_dataRadiographicP_4<-nma_dataRadiographicP_3
  colnames(nma_dataRadiographicP_4)<-c("NS","n.t","mean.t","sd.t","n.c","mean.c","sd.c")
  row.names(nma_dataRadiographicP_4)<-NULL
  m.cont<-metacont(n.e = n.t, 
                   mean.e = mean.t, 
                   sd.e = sd.t, 
                   n.c = n.c, 
                   mean.c = mean.c, 
                   sd.c =sd.c,
                   data = nma_dataRadiographicP_4,
                   sm = "SMD",
                   method.tau = "PM",
                   random = TRUE,
                   hakn = TRUE,
                   title = "Pain (final)"
                   
  )
  
  summary_m.cont<-summary(m.cont)
  contribution_weight<-as.character(round((summary_m.cont$w.random/sum(summary_m.cont$w.random))*100,2))
  pooled_estimate2<-data.frame(Mean=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.cont$tau2,4),
                          Lower=round(summary_m.cont$lower.tau2,4),
                          Upper=round(summary_m.cont$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.cont$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.cont$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.cont$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.cont$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  Mean_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1,temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(5,12)])
  base_data <- tibble(mean  = Mean_individual_studies$MD, 
                      lower = Mean_individual_studies$Lower,
                      upper = Mean_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      treatment1 = paste0(as.character(temp11.2$x.c),"      ","(",as.character(temp11.2$sampleSize.2),")"),
                      treatment2 = paste0(as.character(temp11.2$x.t),"      ","(",as.character(temp11.2$sampleSize.1),")"),
                      SMD = paste(format(round(Mean_individual_studies$MD,2),nsmall=2),paste("(",format(round(Mean_individual_studies$Lower,2),nsmall=2),", ",format(round(Mean_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$Mean, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   SMD = paste(format(round(pooled_estimate2$Mean,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("Radiographic Progression (Final)", "Study"),
                   treatment1 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.2),mtxIR.RadiographicP_actual_subl)),
                   treatment2 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.1),mtxIR.RadiographicP_actual_subl)),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   SMD = c("", "SMD"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  RadiographicP_actual_output_df <- bind_rows(header,
                                              base_data,
                                              empty_row,
                                              pooled)
  
  
  
  
  RadiographicP_actual_output_df2<-rbind(RadiographicP_actual_output_df,RadiographicP_actual_output_df[c(nrow(RadiographicP_actual_output_df)),])
  RadiographicP_actual_output_df2[nrow(RadiographicP_actual_output_df),1]<-I2
  RadiographicP_actual_output_df2[nrow(RadiographicP_actual_output_df),2]<-NA
  RadiographicP_actual_output_df2[nrow(RadiographicP_actual_output_df),c(6:11)]<-c(NA,NA,NA,NA,NA,NA,NA)
  RadiographicP_actual_output_df2_list[[i]]<-RadiographicP_actual_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-RadiographicP_actual_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, SMD,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(-1, 1), 
               xlog = F, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR.RadiographicP_actual_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR.RadiographicP_actual_subl))
    )
  
  
  
}#END



View(RadiographicP_actual_output_df2_list[[1]])
n_studies<-c()
n_participants<-c()

extract_num<-function(data){
  return(str_extract(data,"(?<=\\()\\d+(?=\\))"))
}
countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}




for(i in 1:length(graph_list)){
  df<-RadiographicP_actual_output_df2_list[[i]][-c(1,2,countdown(length(RadiographicP_actual_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df,1,extract_num))))
}

length(graph_list)
graph_list[[1]]
Grade_assessment_dataset_direct_RadiographicP_actual<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR.RadiographicP_actual_subl),
                                                                 selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR.RadiographicP_actual_subl),
                                                                 Outcome="Radiographic Progression",
                                                                 `Effect measure`="SMD",
                                                                 `Number of studies`=n_studies,
                                                                 `Number of participant`=n_participants,
                                                                 ROB="",Indirectness="",Imprecision="",
                                                                 Inconsistency="",
                                                                 Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_RadiographicP_actual)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                                        "Intevention",
                                                                                        "Effect measure",
                                                                                        "Number of studies",
                                                                                        "Number of participant",
                                                                                        "Risk of Bias",
                                                                                        "Publication Bias",
                                                                                        "Certainty of evidence (ignores imprecision)")

View(Grade_assessment_dataset_direct_RadiographicP_actual)
#Saving grade assessment details for pain mtxIR as an object of lists
grade_details_RadiographicP_actual_mtxIR<-list(graph_list,RadiographicP_actual_output_df2_list,Grade_assessment_dataset_direct_RadiographicP_actual,I_square_list)
#saveRDS(grade_details_HAQ_bio,file="List of HAQ Biologics data components for GRADE_OCT31.rds")
saveRDS(grade_details_RadiographicP_actual_mtxIR,file="List of Radiographic Progression MTX-IR-Actual data components for GRADE_AUG2026.rds")


RadiographicP_actual_output_one_study_list <- list()
RadiographicP_actual_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

n_studies<-c()
for(i in 1:length(graph_list)){
  df<-RadiographicP_actual_output_df2_list[[i]][-c(1,2,countdown(length(RadiographicP_actual_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
}

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
RadiographicP_actual_output_one_study_list  <- RadiographicP_actual_output_df2_list[one_idx]
RadiographicP_actual_output_more_than_2_list <- RadiographicP_actual_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_RadiographicP_actual[
    Grade_assessment_dataset_direct_RadiographicP_actual$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_RadiographicP_actual[
    Grade_assessment_dataset_direct_RadiographicP_actual$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_RadiographicP_actual_mtxIR_one_study <- list(
  graph_list_one_study_list,
  RadiographicP_actual_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_RadiographicP_actual_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  RadiographicP_actual_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_RadiographicP_actual_mtxIR_one_study,
  file = "List of Radiographic Progression-Actual Dmard MTX-IR data components for GRADE_one_study_AUG2026.rds"
)
#No more than 2 studies
saveRDS(
  grade_details_RadiographicP_actual_mtxIR_more_than_2,
  file = "List of Radiographic Progression-Actual Dmard MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)


#####################################
#####################################
# inthreeArmtrial
# inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,naive_HAQ_subl)))
# colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
# Grade_assessment_dataset_direct_HAQ
# get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_HAQ$Treatment, 
#                                                 Grade_assessment_dataset_direct_HAQ$Intevention), 
#                                           paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
# datr_list<-list()
# for(i in 1:length(get_index_of_threeArm)){
#   datr<-as.data.frame(HAQ_output_df2_list[get_index_of_threeArm[i]])
#   datr_list[[i]]<-datr[datr$study=="Summary","MD"][2]
# }

# 
# grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
# grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 

glen_direct_evidence_rating_RadioPA<-read.csv("Glen_2026-09-04 18_38_RADIOActual-MTXIR_complete.csv")
glen_direct_evidence_rating_RadioPA2<-glen_direct_evidence_rating_RadioPA[,-c(13:15)]
colnames(glen_direct_evidence_rating_RadioPA2)[1:2]<-colnames(grade_table4)[1:2]
# #update direct evidence report
# glen_direct_evidence_rating_haq2.1<-cbind(Grade_assessment_dataset_direct_HAQ[,c(1:6)],glen_direct_evidence_rating_haq2[,-c(1,2)])


glen_direct_evidence_rating_RadioPA2
View(merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_RadioPA2,
           by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T))
glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_RadioPA2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"Radiographic Progression"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"SMD"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR.RadiographicP_actual_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR.RadiographicP_actual_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}


##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for Radiographic Progression_MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for Radiographic Progression_MTX-IR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR HAQ OUTCOME
initial_gradeTable_RadioPA<-Temp_GRADE5

estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_RadioPA2<-initial_gradeTable_RadioPA[match(all_nma_estimates,initial_gradeTable_RadioPA$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()


#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_RadioPA2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_RadioPA2.1<-cbind(initial_gradeTable_RadioPA2[,c(1:14)],first_order_loop_contri,
                                     initial_gradeTable_RadioPA2[,-c(1:14)])

initial_gradeTable_RadioPA2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_RadioPA2.1,
                                                                all_nma_estimates,
                                                                estimate_contribution_matrix,
                                                                glen_direct_evidence_rating6_onlyDirect,
                                                                top_n = 2))
# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_RadioPA2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_RadioPA2.1<-result_df
initial_gradeTable_RadioPA2.1$combI<-str_replace(initial_gradeTable_RadioPA2.1$comb, ":", " vs ")

initial_gradeTable_RadioPA2.1<-process_dataset(initial_gradeTable_RadioPA2.1, 2)
initial_gradeTable_RadioPA3<-initial_gradeTable_RadioPA2.1[match(initial_gradeTable_RadioPA$comb,initial_gradeTable_RadioPA2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for Radiographic Progression_Actual_MTX-IR outcome_2026.csv",row.names = T)

###########Assessing Intransitivity
sub_data_intransitivity <- initial_gradeTable_RadioPA3[initial_gradeTable_RadioPA3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_RadioPA3[initial_gradeTable_RadioPA3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_RadioPA3.1<-cbind(initial_gradeTable_RadioPA3[,c(1:23)],
                                     Intransitivity="",
                                     `Certainty of evidence (ignores imprecision)-Indirect`="",
                                     initial_gradeTable_RadioPA3[,-c(1:23)])
initial_gradeTable_RadioPA3.1[initial_gradeTable_RadioPA3.1$Indirect!="-",c("Intransitivity",
                                                                            "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                                  "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_RadioPA4<-initial_gradeTable_RadioPA3.1

# #######Creating percentage contribution of both direct and indirect evidence
# ####
# colnames_vector<-colnames(estimate_contribution_matrix)
# grade_temp<-cbind(initial_gradeTable_RadioPA4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_RadioPA4[,-c(1:25)])
# ###Applying the P-values 
# P_values_nodesplit_radiographicP_actualmtxIR2
# grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])
# 
# drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
# P_values_nodesplit_radiographicP_actualmtxIR3 <- P_values_nodesplit_radiographicP_actualmtxIR2[match(drt$comb, P_values_nodesplit_radiographicP_actualmtxIR2$comb), , drop = FALSE]
# grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_radiographicP_actualmtxIR3$p,2)
# 
# subRADIOPASOF<-grade_temp2[grade_temp2$p_value!="-",]
# subcontriMatRadioPA<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subRADIOPASOF$comb,]
# 
# 
# # Initializing dirPer and IndirPer vectors
# dirPer <- numeric(nrow(subcontriMatRadioPA))
# IndirPer <- numeric(nrow(subcontriMatRadioPA))
# 
# # Loop to calculate dirPer and IndirPer
# for (i in 1:nrow(subcontriMatRadioPA)) {
#   dirPer[i] <- subcontriMatRadioPA[i, colnames(subcontriMatRadioPA) == row.names(subcontriMatRadioPA)[i]]
#   IndirPer[i] <- 1- dirPer[i]
# }
# 
# # Print the results
# cat("dirPer:", dirPer, "\n")
# cat("IndirPer:", IndirPer, "\n")
# 
# grade_temp2[match(row.names(subcontriMatRadioPA), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
# grade_temp2[match(row.names(subcontriMatRadioPA), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)
# 
# 
# ##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
# grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
# #When only direct or only indirect estimate is present
# grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
#                                  grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
#                                  ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
#                                         grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
# #When both direct and indirect are present
# fg<-c()
# sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
# for(i in 1:nrow(sub_grade_temp2.1)){
#   if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
#     fg[i]<-"Flag for manual review"
#   }else{
#     if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
#       if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
#       ){
#         fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
#       }else{
#         fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
#         
#       }
#     }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
#       fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
#     }else{
#       fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
#     }
#   }
#   
#   
# }
# 
# 
# grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg
###Applying the P-values, when no p-values is available because no both direct and indirect evidence available
#P_values_nodesplit_haqnaive2
grade_temp<-cbind(initial_gradeTable_RadioPA4[,c(1:25)],DirectPerContri="-",IndirectPerContri="-",initial_gradeTable_RadioPA4[,-c(1:25)])
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))


####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
# incho<-c()
# sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
# for(i in 1:nrow(sub_grade_temp2.2)){
#   if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
#     if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
#       incho[i]<-"Serious"
#     }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
#       incho[i]<-"Not Serious"
#       
#     }else{
#       incho[i]<-"Some Concerns"
#     }
#   }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
#     if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
#       incho[i]<-"Some Concerns"
#     }else{
#       incho[i]<-"Not Serious"
#     } 
#   }else{
#     if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
#       incho[i]<-"Some Concerns"
#     }else{
#       incho[i]<-"Not Serious"
#     } 
#   }
#   
#   
#   
# }
# 
# 
# grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho


initial_gradeTable_RadioPA5<-grade_temp2.2
initial_gradeTable_RadioPA5<-initial_gradeTable_RadioPA5[,-c(30,33)]
colnames(initial_gradeTable_RadioPA5)[c(29,30)]<-c("Initial NMA evidence rating",
                                                   "P-Value")
###After incoherence evaluation
initial_gradeTable_RadioPA5.1<-cbind(initial_gradeTable_RadioPA5[,c(1:31)],
                                     qoeAfterIncohenrence="",
                                     initial_gradeTable_RadioPA5[,-c(1:31)])
initial_gradeTable_RadioPA5.1[(initial_gradeTable_RadioPA5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_RadioPA5.1[(initial_gradeTable_RadioPA5.1$`P-Value`=="-"),"Initial NMA evidence rating"]


#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_RadioPA5.1[(initial_gradeTable_RadioPA5.1$Direct!="-")&(initial_gradeTable_RadioPA5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_RadioPA5.1[(initial_gradeTable_RadioPA5.1$Direct!="-")&(initial_gradeTable_RadioPA5.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_RadioPA6<-cbind(initial_gradeTable_RadioPA5.1,
                                   p.estimate=apply(as.data.frame(initial_gradeTable_RadioPA5.1$NMA),1,function(x){
                                     as.numeric(sub(" .*", "", x))
                                   }),
                                   t(apply(as.data.frame(initial_gradeTable_RadioPA5.1$NMA),1,function(x){
                                     as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                     
                                   })))

colnames(initial_gradeTable_RadioPA6)[c(36,37)]<-c("l.bound","u.bound")

####Running Continuous NMA model
# load data for MTC
##Getting baseline risk using the placebo intervention
radiographicP_actual_data<-radiographicP_actual_mtxIR_data_new_withROB2_indrectness2
radiographicP_actual_data<-radiographicP_actual_data[order(radiographicP_actual_data$study),c(1:6)]

radiographicP_actual_data_wide2<-reshape(number_of_arms(radiographicP_actual_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(radiographicP_actual_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(radiographicP_actual_data_wide2[,c(3,7,11)],1,xcv)

radiographicP_actual_data_wide3<-cbind(radiographicP_actual_data_wide2,na=numberOfArms)
View(radiographicP_actual_data_wide3)
#finding standard errors
radiographicP_actual_data_wide3$SD.1<-radiographicP_actual_data_wide3$SD.1/sqrt(radiographicP_actual_data_wide3$sampleSize.1)
radiographicP_actual_data_wide3$SD.2<-radiographicP_actual_data_wide3$SD.2/sqrt(radiographicP_actual_data_wide3$sampleSize.2)
radiographicP_actual_data_wide3$SD.3<-radiographicP_actual_data_wide3$SD.3/sqrt(radiographicP_actual_data_wide3$sampleSize.3)
#radiographicP_actual_data_wide3$SD.4<-radiographicP_actual_data_wide3$SD.4/sqrt(radiographicP_actual_data_wide3$sampleSize.4)
colnames(radiographicP_actual_data_wide3)<-c("study", "study_name", "t.1","se.1","n.1", "x.1",
                                             "t.2","se.2","n.2","x.2","t.3","se.3","n.3","x.3", "na")
                                             #"t.4","se.4","n.4","x.4","na" )
MTCData_cont <- as.data.frame(radiographicP_actual_data_wide3)
x <- data.matrix(MTCData_cont[,c(6,10,14)])
n <- data.matrix(MTCData_cont[,c(5,9,13)])
t <- data.matrix(apply(MTCData_cont[,c("t.1","t.2","t.3")],c(1,2),as.numeric))
se <- data.matrix(MTCData_cont[,c(4,8,12)])
na <- data.matrix(MTCData_cont[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- x[t[,'t.1']==1,"x.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
sePBO <- se[t[,'t.1']==1,"se.1"]

#pooled sd
radiographicP_actual_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2")]
# calculate_pooled_sd <- function(row) {
#   sqrt(sum((row[c("sampleSize.1", "sampleSize.2")] - 1) * row[c("SD.1", "SD.2")]^2) / sum(row[c("sampleSize.1", "sampleSize.2")] - 1))
# }
calculate_pooled_sd <- function(row) {
  
  sd_cols <- c("SD.1", "SD.2", "SD.3")
  n_cols  <- c("sampleSize.1", "sampleSize.2",
               "sampleSize.3")
  
  sd <- as.numeric(row[sd_cols])
  n  <- as.numeric(row[n_cols])
  
  keep <- !is.na(sd) & !is.na(n)
  
  sqrt(
    sum((n[keep] - 1) * sd[keep]^2) /
      sum(n[keep] - 1)
  )
}
# Apply the function row-wise
# pooled.sd=data.matrix(cbind(apply(radiographicP_actual_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2")], 1, calculate_pooled_sd),
#                             apply(radiographicP_actual_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2")], 1, calculate_pooled_sd)))
# 
# data_mean.pooled.sd.ref<-radiographicP_actual_data_wide2[radiographicP_actual_data_wide2$treatment.1==1,c("SD.1", "sampleSize.1","SD.2", "sampleSize.2")]
# mean.pooled.sd.ref<- calculate_pooled_sd(data_mean.pooled.sd.ref)

pooled.sd <- data.matrix(
  cbind(
    apply(
      radiographicP_actual_data_wide2[
        , c("SD.1","sampleSize.1",
            "SD.2","sampleSize.2",
            "SD.3","sampleSize.3")
      ],
      1,
      calculate_pooled_sd
    ),
    apply(
      radiographicP_actual_data_wide2[
        , c("SD.1","sampleSize.1",
            "SD.2","sampleSize.2",
            "SD.3","sampleSize.3")
      ],
      1,
      calculate_pooled_sd
    ),
    apply(
      radiographicP_actual_data_wide2[
        , c("SD.1","sampleSize.1",
            "SD.2","sampleSize.2",
            "SD.3","sampleSize.3")
      ],
      1,
      calculate_pooled_sd
    ),
    apply(
      radiographicP_actual_data_wide2[
        , c("SD.1","sampleSize.1",
            "SD.2","sampleSize.2",
            "SD.3","sampleSize.3")
      ],
      1,
      calculate_pooled_sd
    )
  )
)

data_mean.pooled.sd.ref <-
  radiographicP_actual_data_wide2[
    radiographicP_actual_data_wide2$treatment.1 == 1,
    c("SD.1","sampleSize.1",
      "SD.2","sampleSize.2",
      "SD.3","sampleSize.3")
  ]

mean.pooled.sd.ref <- calculate_pooled_sd(
  data_mean.pooled.sd.ref[1, ]
)
# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
#model.file4<-"Cont_multiarm_RE.txt"
model.file4<-"tt.txt"#"Cont_new_multiarm_SMD_RE.txt"
n.chains=3
model.inits <- initv1 

cont.model3<-jags.model(file=model.file4,
                        data= list("x"=x,"t"=t,"se"=se,"na" = na[,1], "NT" = NT, "NS" = NS,
                                   "NTpbo"= NTpbo,"xPBO"=xPBO,"sePBO"=sePBO,
                                   "pooled.sd"=pooled.sd,"mean.pooled.sd.ref"=mean.pooled.sd.ref),
                        n.chains=n.chains,
                        inits = model.inits,
                        n.adapt=5000
)

update(object=cont.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("diff","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( cont.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]


#Converting estimates relative to baseline to absolutes
estimates_relative_to_baseline<-initial_gradeTable_RadioPA6[1:(NT-1),c(1,2,35:37)]
typeof(estimates_relative_to_baseline)

absolutes_from_baseline<-cbind(baseline_risk =round(baseline_risk,2),round(baseline_risk,2)+estimates_relative_to_baseline[,c(3:5)])
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(round(baseline_risk,2),estimates_relative_to_baseline_combine[,7]))
colnames(all_baseline)<-"base"

com.list<-list()
unique_treat2<-unique(initial_gradeTable_RadioPA6$treat2)
d.sub<-initial_gradeTable_RadioPA6[initial_gradeTable_RadioPA6$treat2==unique_treat2[1],c(1,2,35:37)]
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_RadioPA6[initial_gradeTable_RadioPA6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,base+d.sub[,c(3:5)])
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe<-com.dataframe[c(1:nrow(initial_gradeTable_RadioPA6)),]
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9])

higher_lower_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"lower"
  }else{
    text<-"higher"
  }
  return(paste(abs(round(e1-e2,2)),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,13])))


prefix.assumedrisk <- ' Mean '
prefix.corrrisk<-' Mean '
prefix.md<- " SMD "
suffix.assumedrisk <- ''
suffix.corrrisk <- ''
time<-1


result_com.dataframe_withtext2<-cbind(com.dataframe_withtext[,1:2],
                                      paste(prefix.assumedrisk,com.dataframe_withtext[,6], suffix.assumedrisk, sep=''),
                                      paste(prefix.corrrisk, com.dataframe_withtext[,7], suffix.corrrisk, ' (', com.dataframe_withtext[,8], ' to ',com.dataframe_withtext[,9],')', sep=''),
                                      paste(prefix.md, com.dataframe_withtext[,10], ' (', com.dataframe_withtext[,11], ' to ',com.dataframe_withtext[,12],')', sep=''),
                                      # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                      com.dataframe_withtext[,c(3:5,10:12)]
)

colnames(result_com.dataframe_withtext2)[c(3:11)]<-c("Assumed mean (comparator)",	"Corresponding mean (intervention)",
                                                     "Mean difference","p.estimate","l.bound","u.bound","","","")



initial_gradeTable_RadioPA7<-cbind(initial_gradeTable_RadioPA6[,-c(35:37)],
                                   result_com.dataframe_withtext2[,-c(1,2)],pair="")
#colnames(initial_gradeTable_RadioPA7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_RadioPA8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_RadioPA7$treat1,mtxIR.RadiographicP_actual_subl)),
                                   `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_RadioPA7$treat2,mtxIR.RadiographicP_actual_subl)),
                                   initial_gradeTable_RadioPA7)
colnames(initial_gradeTable_RadioPA8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_RadioPA8<-initial_gradeTable_RadioPA8[,-c(46)]
write.csv(initial_gradeTable_RadioPA8,"Radiographic Progression MTX-IR-Actual SoF table without imprecision assessment AUG2026.csv",row.names = F)



effect.size.rating_list<-lapply(split(initial_gradeTable_RadioPA8, seq_len(nrow(initial_gradeTable_RadioPA8))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 5, large.effect = 10, framework = T)
})



as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_RadioPA9<-initial_gradeTable_RadioPA8
initial_gradeTable_RadioPA9<-cbind(initial_gradeTable_RadioPA9[,c(1:34,37:45,35:36)])
initial_gradeTable_RadioPA9.1<-cbind(initial_gradeTable_RadioPA9[,c(1:37)],
                                     Small.effect.threshold= 5,
                                     Large.effect.threshold= 10,
                                     Framework = "NULL",
                                     Imprecision.1 = "",
                                     Certainty.of.evidence = "",
                                     initial_gradeTable_RadioPA9[,c(44,45)])

initial_gradeTable_RadioPA9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_RadioPA_mtxIR<-initial_gradeTable_RadioPA9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_RadioPA_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_RadioPA_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_RadioPA_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_RadioPA_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_RadioPA_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_RadioPA_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_RadioPA_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_RadioPA_mtxIR[is.na(detailed_SOF_RadioPA_mtxIR$Certainty.of.evidence),])
detailed_SOF_RadioPA_mtxIR[is.na(detailed_SOF_RadioPA_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"


write.csv(detailed_SOF_RadioPA_mtxIR,"Radiographic Progression MTX-IR-Actual_Detail_SoF_AUG2026.csv", row.names = F)

###############################################################################
###############################################################################
#########################################################################################################
#########################################################################################################

##Change in baseline studies, with no actuals
radiog_Cbaseline<-final_mtxIRDataset_radiographicProgression2
View(radiog_Cbaseline)
View(radiog_Cbaseline[radiog_Cbaseline$Study_no==140,])
tttt_studyInfo<-radiog_Cbaseline[radiog_Cbaseline$Study_no==124,c(1:3)]
tttt<-radiog_Cbaseline[radiog_Cbaseline$Study_no==124,c(52:141)]
tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
getting_Cbaseline<-tttt[colSums(!is.na(tttt)) > 0]
getting_Cbaseline
colnames(getting_Cbaseline)
getting_Cbaseline2<-getting_Cbaseline[, grepl("change.from.baseline",names(getting_Cbaseline))]
getting_Cbaseline2
View(getting_Cbaseline2)


##Nam_IDEA 2014
months12_cols<-grep(
  "^Radiographic.Progression\\.12\\.months",
  names(radiog_Cbaseline),
  value = TRUE
)
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Nam_IDEA 2014",months12_cols]<-NA
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Nam_IDEA 2014","Radiographic.Progression.6.months...change.from.baseline.N"]<-c(57,55)
#Leirisalo-Repo 2013
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Leirisalo-Repo 2013",c("Radiographic.Progression.6.months...change.from.baseline.SD","Radiographic.Progression.6.months...change.from.baseline.SD_edit")]<-exactp.sd(0.004,46,45,1.4,-0.3)[2]

#Boers 1997
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Boers 1997",months12_cols]<-NA
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Boers 1997",c("Radiographic.Progression.6.months...change.from.baseline.SD","Radiographic.Progression.6.months...change.from.baseline.SD_edit")]<-round(c(iqr2sd(0,44),iqr2sd(0,28)),2)
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Boers 1997","Radiographic.Progression.6.months...change.from.baseline.Mean_edit"]<-radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Boers 1997","Radiographic.Progression.6.months...change.from.baseline.Median"]

#Proudman 2000 
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Proudman 2000","Radiographic.Progression.12.months..change.from.baseline.Mean_edit"]<-radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Proudman 2000","Radiographic.Progression.12.months..change.from.baseline.Median"]

#Drosos 1998
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Drosos 1998",c("Radiographic.Progression.End.of.trial.change.from.baseline.SD","Radiographic.Progression.End.of.trial.change.from.baseline.SD_edit")]<- paste0("8888",round(c(ci2sd(0.9489,6.3669,51),ci2sd(0.7633,6.1521,52)),2))

##Combe 2021
months12_cols<-grep(
  "^Radiographic.Progression\\.12\\.months",
  names(radiog_Cbaseline),
  value = TRUE
)
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Combe 2021",months12_cols]<-NA

##Taylor 2017
months12_cols<-grep(
  "^Radiographic.Progression\\.12\\.months",
  names(radiog_Cbaseline),
  value = TRUE
)
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Taylor 2017",months12_cols]<-NA

#Peterfy 2016
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Peterfy 2016","Radiographic.Progression.12.months..change.from.baseline.SD_edit"] <- exactp.sd(0.002,63,60,1.37,0.29)[2]

#Matsubara 2018, no SD
radiog_Cbaseline[radiog_Cbaseline$Study.Identifier=="Matsubara 2018",months12_cols]<-NA

temp_data_radgCbaseline<-list()
out_entries_radg<-list()
empty_entries_radg<-list()
study_changeBRadio<-unique(radiog_Cbaseline$Study_no)



for(i in 1:length(unique(radiog_Cbaseline$Study_no))){
  tttt<-radiog_Cbaseline[radiog_Cbaseline$Study_no==study_changeBRadio[i],c(52:141)]
  tttt_studyInfo<-radiog_Cbaseline[radiog_Cbaseline$Study_no==study_changeBRadio[i],c(1:3)]
  tttt<-cbind(tttt[,21:60],tttt[,1:20],tttt[,71:90],tttt[61:70])
  
  if(sum(which(colSums(!is.na(tttt)) > 0))==0){
    empty_entries_radg[[i]]<-radiog_Cbaseline[radiog_Cbaseline$Study_no==study_changeBRadio[i],c(1:3,52:141)]
  }else if(!(any(grepl("change.from.baseline",
                       colnames(tttt[colSums(!is.na(tttt)) > 0]))))){
    out_entries_radg[[i]]<-radiog_Cbaseline[radiog_Cbaseline$Study_no==study_changeBRadio[i],c(1:3,52:141)]
  }else{
    ##actual only but need to choose month
    getting_actuals<-tttt[colSums(!is.na(tttt)) > 0]
    getting_actuals2<-getting_actuals[, grepl("change.from.baseline",names(getting_actuals))]
    if(ncol(getting_actuals2)==5){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[3])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==6){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[4])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==7){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[5])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==8){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[6])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==10){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[8])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==11){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[9])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==15){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[13])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==18){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[17])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }else if(ncol(getting_actuals2)==20){
      eee<-cbind(tttt_studyInfo,getting_actuals2%>%select(matches("\\.(SD_edit|N|Mean_edit)$"))%>%
                   select(
                     matches("\\.SD_edit$"),
                     matches("\\.N$"),
                     matches("\\.Mean_edit$")
                   ),timepoint=colnames(getting_actuals2)[18])
      colnames(eee)<-c("Study_no","Study.Identifier","Intervention","SD","sampleSize","Mean","Timepoint")
      temp_data_radgCbaseline[[i]]<- eee
    }
  }
  
  
}

View(do.call(rbind,temp_data_radgCbaseline))
View(do.call(rbind,out_entries_radg))
View(do.call(rbind,empty_entries_radg))


RadiographicPCB_mtxIR_dataset_analysis<-do.call(rbind,temp_data_radgCbaseline)
RadiographicPCB_mtxIR_dataset_analysis1<-RadiographicPCB_mtxIR_dataset_analysis

RadiographicPCB_mtxIR_dataset_analysis1<- RadiographicPCB_mtxIR_dataset_analysis1%>%
  filter(!Study.Identifier=="Fleischmann 2023.a")

row.names(RadiographicPCB_mtxIR_dataset_analysis1)<-NULL

##Removing 8888
RadiographicPCB_mtxIR_dataset_analysis1$SD<-round(as.numeric(gsub(8888,'',RadiographicPCB_mtxIR_dataset_analysis1$SD)),2)

#####Formatting the month column
RadiographicPCB_mtxIR_dataset_analysis1$Timepoint <- gsub(
  "Radiographic.Progression\\.|\\.months\\.{2,3}change.from\\.baseline\\.(N|Mean_edit|Median|SD|SD_edit|Lower\\.CI)|\\.(Mean_edit|Median)",
  "",
  RadiographicPCB_mtxIR_dataset_analysis1$Timepoint
)

RadiographicPCB_mtxIR_dataset_analysis1<- RadiographicPCB_mtxIR_dataset_analysis1%>%
  dplyr::rename("treatment" = Intervention)

##removing choy dataset, its making the model unstable
##RadiographicPCB_mtxIR_dataset_analysis1<-RadiographicPCB_mtxIR_dataset_analysis1[!RadiographicPCB_mtxIR_dataset_analysis1$Study.Identifier=="Choy 2008",]

##Drop the high risk of bias studies
RadiographicPCB_mtxIR_dataset_analysis1<- RadiographicPCB_mtxIR_dataset_analysis1%>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012"))

RadiographicPCB_mtxIR_dataset_analysis1$Study_no<-renumber(RadiographicPCB_mtxIR_dataset_analysis1,"Study_no")
row.names(RadiographicPCB_mtxIR_dataset_analysis1)<-NULL

library(dplyr)
for(i in 1:nrow(RadiographicPCB_mtxIR_dataset_analysis1)){
  ppp<- master.medtable%>% filter(RadiographicPCB_mtxIR_dataset_analysis1$treatment[i] ==master.medtable$medname )
  RadiographicPCB_mtxIR_dataset_analysis1$treatment[i]<-ppp[1,1]
}

##Creating initial legend
RadiographicPCB.mtxIR.medlegend <- new.medtable(RadiographicPCB_mtxIR_dataset_analysis1,'treatment')
RadiographicPCB.mtxIR.medlegend

##NMA plot reason for narrowing of network
nma.plot(RadiographicPCB_mtxIR_dataset_analysis1[,1],RadiographicPCB_mtxIR_dataset_analysis1[,3],data = RadiographicPCB_mtxIR_dataset_analysis1, trtname = as.vector(RadiographicPCB.mtxIR.medlegend$medname))

################################################################################################################
#trying the narrowing a network function
#############################################################################################
RadiographicPCB_mtxIR_dataset_analysis1

RadiographicPCB_mtxIR_dataset_analysis1 <- RadiographicPCB_mtxIR_dataset_analysis1 %>%
  dplyr:: rename("study" = Study_no,"study_name" = Study.Identifier)

#Narrow 1
mtxIR_RadiographicPCB_narrow<-narrow.all(RadiographicPCB_mtxIR_dataset_analysis1)
mtxIR_RadiographicPCB_narrow

#sub network plot 1
mtxIR.RadiographicPCB.d<-mtxIR_RadiographicPCB_narrow[[1]]
length(unique(mtxIR.RadiographicPCB.d$treatment))
mtxIR_RadiographicPCB_subl<-new.medtable(mtxIR.RadiographicPCB.d,'treatment')

##To download acr50 data
mtxIR.RadiographicPCB.d_data_download<-mtxIR.RadiographicPCB.d
mtxIR.RadiographicPCB.d_data_download$treatment<-renumber(mtxIR.RadiographicPCB.d_data_download, 'treatment')
mtxIR.RadiographicPCB.d_data_download$treatment<-selected_treatment_name(mtxIR.RadiographicPCB.d_data_download$treatment,mtxIR_RadiographicPCB_subl)
write.csv(mtxIR.RadiographicPCB.d_data_download,"Radiographic progression Change from baseline DMARD MTX-IR dataset_AUG2026.csv",row.names = F)
# 1. Open jpeg file
jpeg("Radiographic Progression CB_mtxIR Network plot with named interventions_JUNE2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.RadiographicPCB.d[,1],mtxIR.RadiographicPCB.d[,3],data = mtxIR.RadiographicPCB.d, trtname = as.vector(mtxIR_RadiographicPCB_subl$medname),title = "Radiographic Progression CB mtxIR Network plot")
# Close the jpeg file
dev.off() 
mtxIR.RadiographicPCB.d1<-mtxIR.RadiographicPCB.d
mtxIR.RadiographicPCB.d1$treatment<-renumber(mtxIR.RadiographicPCB.d1,"treatment")
# 1. Open jpeg file
jpeg("Radiographic Progression CB_MTX-IR Network plot without named interventions_AUG2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.RadiographicPCB.d1[,1],mtxIR.RadiographicPCB.d1[,3],data = mtxIR.RadiographicPCB.d1,title = "Radiographic Progression CB mtxIR Network plot1")
# Close the jpeg file
dev.off() 
##############################################################################
##############################################################################
#Performing sub network plot 1
library(coda)
library(gemtc)
library(netmeta)
library(tidyr)
library(netmeta)
library(magrittr)
library(tidyr)
View(mtxIR.RadiographicPCB.d)
mtxIR_RadiographicPCB_subl

network_data_RadiographicPCB_mtxIR_1<-mtxIR.RadiographicPCB.d%>%
  select("study", "study_name","treatment","SD", "sampleSize","Mean","Timepoint" )
network_data_RadiographicPCB_mtxIR_1$treatment<-renumber(network_data_RadiographicPCB_mtxIR_1,"treatment")
network_data_RadiographicPCB_mtxIR_1$study<-renumber(network_data_RadiographicPCB_mtxIR_1,"study")
network_data_RadiographicPCB_mtxIR_1<-network_data_RadiographicPCB_mtxIR_1[order(network_data_RadiographicPCB_mtxIR_1$study),]
row.names(network_data_RadiographicPCB_mtxIR_1)<-NULL
nma.plot(network_data_RadiographicPCB_mtxIR_1$study,network_data_RadiographicPCB_mtxIR_1$treatment,data = network_data_RadiographicPCB_mtxIR_1, trtname = as.vector(mtxIR_RadiographicPCB_subl$medname),title = "Radiographic Progression mtxIR Network plot1")

network_data_RadiographicPCB_mtxIR_1<-network_data_RadiographicPCB_mtxIR_1%>%
  dplyr::rename("std.dev" = SD, "mean"=Mean)

p_RadiographicPCB_SMD<-pairwise(treat=treatment,mean=mean,sd=std.dev,data=network_data_RadiographicPCB_mtxIR_1,
                                studlab = study,n=sampleSize,sm="SMD")
#colnames(p_RadiographicPCB_SMD)[c(1,4:5)]<-c("study","diff", "std.err")
p_RadiographicPCB_SMD<-p_RadiographicPCB_SMD %>%
  dplyr::rename("diff"= TE, "std.err"=seTE)

p_RadiographicPCB_SMD %>%
  dplyr::select(c(study_name,diff, std.err,treat1,treat2)) %>%
  pivot_longer(-study_name,
               names_to = c(".value"),
               names_pattern = "(..)") %>% 
  set_colnames(c("study", "diff", 
                 "std.err", "treatment"))-> data5


data5 <- clean_and_impute_se(data5)

network_data_radiographicPCB_mtxIR_1<-mtc.network(data.re=data5)
summary(network_data_radiographicPCB_mtxIR_1)
View(network_data_radiographicPCB_mtxIR_1$data.re)


##RUN the model
set.seed(5811122)
network_data_radiographicPCB_mtxIR_1.model <- mtc.model(network_data_radiographicPCB_mtxIR_1, likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

# By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_radiographicPCB_mtxIR_1.model$inits <- mapply(c, network_data_radiographicPCB_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_radiographicPCB<-mtc.run(network_data_radiographicPCB_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_radiographicPCB)

View(as.data.frame(relative.effect.table(results_mtxIR_radiographicPCB)))
View(as.data.frame(relative.effect.table.2DP(results_mtxIR_radiographicPCB)))
refME_mtxIR_radiographicPCB<-as.data.frame(relative.effect.table.2DP2(results_mtxIR_radiographicPCB, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refME_mtxIR_radiographicPCB)<-mapvalues(row.names(refME_mtxIR_radiographicPCB),
                                                  from = c(1:(nrow(mtxIR_RadiographicPCB_subl))),
                                                  to = as.vector(mtxIR_RadiographicPCB_subl$medname))
colnames(refME_mtxIR_radiographicPCB)<-mapvalues(colnames(refME_mtxIR_radiographicPCB),
                                                 from = c(1:(nrow(mtxIR_RadiographicPCB_subl))),
                                                 to = as.vector(mtxIR_RadiographicPCB_subl$medname))
View(refME_mtxIR_radiographicPCB)

index_diagonal5<-data.frame(index1=1:nrow(refME_mtxIR_radiographicPCB),index2=1:nrow(refME_mtxIR_radiographicPCB))
options(warn=-1)
for(i in 1:nrow(refME_mtxIR_radiographicPCB)){
  refME_mtxIR_radiographicPCB[index_diagonal5$index1[i],index_diagonal5$index2[i]]<-selected_treatment_name(as.numeric(refME_mtxIR_radiographicPCB[row(refME_mtxIR_radiographicPCB)==col(refME_mtxIR_radiographicPCB)])[i],mtxIR_RadiographicPCB_subl)
}
options(warn=0)
View(refME_mtxIR_radiographicPCB)

##########
write.csv(refME_mtxIR_radiographicPCB, "Radiographic Progression-Change in Baseline pairwise comparisons_MTX-IR_AUG2026.csv", row.names = T)
#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_AUG2026.xlsx")

addWorksheet(wb, "Radiographic Progression-CB")
writeData(wb, "Radiographic Progression-CB", refME_mtxIR_radiographicPCB, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_AUG2026.xlsx", overwrite = TRUE)

##GRADE Assessment
nma_estimate_mtxIR_radiographicPCB_dataframe<-nma.estimate.table.MD(results_mtxIR_radiographicPCB)#use results
View(nma_estimate_mtxIR_radiographicPCB_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
# network_data_radiographicPCB_mtxIR_1.ns<-mtc.nodesplit(network_data_radiographicPCB_mtxIR_1, thin=50,likelihood='normal', link='identity',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 
# 
# names(network_data_radiographicPCB_mtxIR_1.ns)
# mtc.nodesplit.comparisons(network_data_radiographicPCB_mtxIR_1.ns)
# 
# # Overall summary and plot
# summary.network_data_radiographicPCB_mtxIR_1.ns <- summary(network_data_radiographicPCB_mtxIR_1.ns)
# 
# #Direct effects
# direct.data<-summary.network_data_radiographicPCB_mtxIR_1.ns$dir.effect
# 
# #Indirect effects
# indirect.data<-summary.network_data_radiographicPCB_mtxIR_1.ns$ind.effect
# 
# #Network estimate
# network.data<-summary.network_data_radiographicPCB_mtxIR_1.ns$cons.effect
# 
##Acquiring both direct and indirect estimate from gemtc
# grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_radiographicPCB_dataframe,
#                     by=c("treat1","treat2"),all.x=T,all.y=T)




#Adding ROB Data
radioCB_mtxIR<-unique(network_data_RadiographicPCB_mtxIR_1$study_name)
data_radioCB_mtxIR<-network_data_RadiographicPCB_mtxIR_1

current_assessment_for_studies<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%radioCB_mtxIR,]
#current_assessment_for_studies<-current_assessment_for_studies[,c(1,15,2,6,5,4,3,7:9,14,10,12,13,11)]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies<- current_assessment_for_studies %>%
  select(c(col_names))
list_appraisal<-list()
for(i in 1:length(radioCB_mtxIR)){
  temp11<-data_radioCB_mtxIR[data_radioCB_mtxIR$study_name==radioCB_mtxIR[i],]
  temp12<-current_assessment_for_studies[current_assessment_for_studies$Study.Id==radioCB_mtxIR[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal[[i]]<-temp13
}
radioCB_mtxIR_withROB<-do.call(rbind,list_appraisal)
View(radioCB_mtxIR_withROB)
radioCB_mtxIR_withROB<- radioCB_mtxIR_withROB%>%
  dplyr:: rename("SD" = std.dev, "Mean" = mean)
radioCB_mtxIR_withROB2<-radioCB_mtxIR_withROB%>%
  select(c("study","study_name", "treatment", "SD", "sampleSize",
           "Mean","Timepoint","OVERALL..ALL.OTHER.OUTCOMES"))
View(radioCB_mtxIR_withROB2)

###Adding Indirectness data
study_unique<-unique(radioCB_mtxIR_withROB2$study_name)
data_nn<-radioCB_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp%>%
  dplyr::rename( "Consensus"= `Final Rating`)

View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]



list_appraisal2<-list()

for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$study_name==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(length(temp12)==0){
    temp13<-cbind(temp11,temp12="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
radioCB_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal2)
radioCB_mtxIR_withROB2_indrectness<- radioCB_mtxIR_withROB2_indrectness%>%
  dplyr::rename(
    ROB = OVERALL..ALL.OTHER.OUTCOMES,
    Indirectness = Consensus
  )

radioCB_mtxIR_withROB2_indrectness$Indirectness<-tolower(radioCB_mtxIR_withROB2_indrectness$Indirectness)
View(radioCB_mtxIR_withROB2_indrectness)

radioCB_mtxIR_data_new_withROB2_indrectness2<-radioCB_mtxIR_withROB2_indrectness
radioCB_mtxIR_data_new_withROB2_indrectness2$treatment<-renumber(radioCB_mtxIR_data_new_withROB2_indrectness2,"treatment")
radioCB_mtxIR_data_new_withROB2_indrectness2$ROB<-ifelse(radioCB_mtxIR_data_new_withROB2_indrectness2$ROB=="unclear","some concerns",radioCB_mtxIR_data_new_withROB2_indrectness2$ROB)

p_RadiographicPCB_SMD2<-p_RadiographicPCB_SMD
p_RadiographicPCB_SMD2<-p_RadiographicPCB_SMD2%>%
  dplyr:: rename("TE"=diff,"seTE"=std.err)


net1 <- netmeta(p_RadiographicPCB_SMD2)
dat.net<-netsplit(net1)
View(dat.net$indirect.common)
View(cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)],dat.net$compare.common[,c(2:5)]))
View(dat.net)
dat.net.contrib<-netcontrib(net1, method = "r")

direct_evidence_new<-cbind(dat.net$direct.common[,c(1:5)],dat.net$indirect.common[,c(2:5)])
View(direct_evidence_new)
#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")

grade_table1<-cbind(nma_estimate_mtxIR_radiographicPCB_dataframe[,c(1,2)],
                    as.data.frame(matrix(NA,nrow(nma_estimate_mtxIR_radiographicPCB_dataframe),6)),
                    nma_estimate_mtxIR_radiographicPCB_dataframe[,-c(1,2)])
colnames(grade_table1)[3:8]<-c("MD/SMD1","low1","high1","MD/SMD2","low2","high2")

#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]


##Getting direct evidence for intervention in only multi-arm trials
zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
colnames(zzz)<-c("treat1", "treat2")
grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]


grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_RadiographicPCB_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_RadiographicPCB_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))


#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
RadiographicPCB_output_df2_list<-list()
I_square_list<-list()

for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-radioCB_mtxIR_data_new_withROB2_indrectness2[radioCB_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select(study, Mean.1, SD.1, Mean.2, SD.2 )
  temp11.1[,c("SD.1","SD.2")]<-temp11.1[,c("SD.1","SD.2")]/(sqrt(temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c("sampleSize.1" , "sampleSize.2")]))
  
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","x.t","se.t","x.c","se.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  
  
  
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    ##For one study
    n_i<-temp11[k,c("sampleSize.1" , "sampleSize.2")]
    sds<-temp11[k,c("SD.1","SD.2")]
    
    # Calculate the pooled SD for SMD
    pooled_sd <- sqrt(sum((n_i - 1) * sds^2) / sum(n_i - 1))
    model.file<-"Cont_direct_SMD_onestudy.txt"
    #model.file<-"Cont_direct_onestudy.txt"
    n.chains=3
    set.seed(66451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c,
                                   "pooled.sd.c"=pooled_sd,
                                   "pooled.sd.t"=pooled_sd
                         ),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    # params <- c("diff")
    params <- c("SMD")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    fgh[[k]]<-data.frame(MD=summary_entry1$statistics[1],
                         Lower=summary_entry1$quantiles[c(1)],
                         Upper=summary_entry1$quantiles[c(5)])
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
  }else{
    temp12<-temp11.1
    
    data4sdPool<-temp11[c("sampleSize.1","SD.1","sampleSize.2","SD.2")]
    # Define a function to calculate the pooled standard deviation
    calculate_pooled_sd <- function(row) {
      sqrt(sum((row[c("sampleSize.1", "sampleSize.2")] - 1) * row[c("SD.1", "SD.2")]^2) / sum(row[c("sampleSize.1", "sampleSize.2")] - 1))
    }
    
    # Apply the function row-wise
    data4sdPool$pooled_sd <- apply(data4sdPool, 1, calculate_pooled_sd)
    
    #model.file<-"Cont_direct_multiplestudy_RE.txt" #Check this again  #"Cont_direct_SMD_multiplestudy_RE.txt"
    model.file<- "Cont_direct_SMD_multiplestudy_RE.txt"
    
    n.chains=3
    set.seed(77451)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("x.t"=temp12$x.t,
                                   "se.t"=temp12$se.t,"x.c"=temp12$x.c,"se.c"=temp12$se.c, "NS"=length(temp12$NS),
                                   "pooled.sd.c"=data4sdPool[,5],#rep(plllsd,length(temp12$NS)) ,#rep(pooled_sd.c,length(temp12$NS)) ,
                                   "pooled.sd.t"=data4sdPool[,5]#rep(plllsd,length(temp12$NS))#rep(pooled_sd.t,length(temp12$NS))
                         ),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    #params <- c("diff")
    params <- c("SMD")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry<-summary(samps)
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_dataRadiographicP_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),c(1,5,6,4,12,13,11)]
  #nma_dataDAS28_3[,c("SD.1","SD.2")]<-nma_dataDAS28_3[,c("SD.1","SD.2")]/(sqrt(nma_dataDAS28_3[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(2,5)]))
  nma_dataRadiographicP_4<-nma_dataRadiographicP_3
  colnames(nma_dataRadiographicP_4)<-c("NS","n.t","mean.t","sd.t","n.c","mean.c","sd.c")
  row.names(nma_dataRadiographicP_4)<-NULL
  m.cont<-metacont(n.e = n.t, 
                   mean.e = mean.t, 
                   sd.e = sd.t, 
                   n.c = n.c, 
                   mean.c = mean.c, 
                   sd.c =sd.c,
                   data = nma_dataRadiographicP_4,
                   sm = "SMD",
                   method.tau = "PM",
                   random = TRUE,
                   hakn = TRUE,
                   title = "Radiographic Progression (CB)"
                   
  )
  
  summary_m.cont<-summary(m.cont)
  contribution_weight<-as.character(round((summary_m.cont$w.random/sum(summary_m.cont$w.random))*100,2))
  pooled_estimate2<-data.frame(Mean=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.cont$tau2,4),
                          Lower=round(summary_m.cont$lower.tau2,4),
                          Upper=round(summary_m.cont$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.cont$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.cont$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.cont$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.cont$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  Mean_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1,temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),c(5,12)])
  base_data <- tibble(mean  = Mean_individual_studies$MD, 
                      lower = Mean_individual_studies$Lower,
                      upper = Mean_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      treatment1 = paste0(as.character(temp11.2$x.c),"      ","(",as.character(temp11.2$sampleSize.2),")"),
                      treatment2 = paste0(as.character(temp11.2$x.t),"      ","(",as.character(temp11.2$sampleSize.1),")"),
                      SMD = paste(format(round(Mean_individual_studies$MD,2),nsmall=2),paste("(",format(round(Mean_individual_studies$Lower,2),nsmall=2),", ",format(round(Mean_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$Mean, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   SMD = paste(format(round(pooled_estimate2$Mean,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("Radiographic Progression (CB)", "Study"),
                   treatment1 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.2),mtxIR_RadiographicPCB_subl)),
                   treatment2 = c("mean  (n)",selected_treatment_name(unique(temp11$treatment.1),mtxIR_RadiographicPCB_subl)),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   SMD = c("", "SMD"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  RadiographicPCB_output_df <- bind_rows(header,
                                         base_data,
                                         empty_row,
                                         pooled)
  
  
  
  
  RadiographicPCB_output_df2<-rbind(RadiographicPCB_output_df,RadiographicPCB_output_df[c(nrow(RadiographicPCB_output_df)),])
  RadiographicPCB_output_df2[nrow(RadiographicPCB_output_df),1]<-I2
  RadiographicPCB_output_df2[nrow(RadiographicPCB_output_df),2]<-NA
  RadiographicPCB_output_df2[nrow(RadiographicPCB_output_df),c(6:11)]<-c(NA,NA,NA,NA,NA,NA,NA)
  RadiographicPCB_output_df2_list[[i]]<-RadiographicPCB_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-RadiographicPCB_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, SMD,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(-1, 1), 
               xlog = F, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_RadiographicPCB_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_RadiographicPCB_subl))
    )
  
  
  
}#END



View(RadiographicPCB_output_df2_list[[3]])
n_studies<-c()
n_participants<-c()

extract_num<-function(data){
  return(str_extract(data,"(?<=\\()\\d+(?=\\))"))
}
countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}




for(i in 1:length(graph_list)){
  df<-RadiographicPCB_output_df2_list[[i]][-c(1,2,countdown(length(RadiographicPCB_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df,1,extract_num))))
}

length(graph_list)
graph_list[[3]]
Grade_assessment_dataset_direct_RadiographicPCB<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_RadiographicPCB_subl),
                                                            selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_RadiographicPCB_subl),
                                                            Outcome="Radiographic Progression_CB",
                                                            `Effect measure`="SMD",
                                                            `Number of studies`=n_studies,
                                                            `Number of participant`=n_participants,
                                                            ROB="",Indirectness="",Imprecision="",
                                                            Inconsistency="",
                                                            Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_RadiographicPCB)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                                   "Intevention",
                                                                                   "Effect measure",
                                                                                   "Number of studies",
                                                                                   "Number of participant",
                                                                                   "Risk of Bias",
                                                                                   "Publication Bias",
                                                                                   "Certainty of evidence (ignores imprecision)")

View(Grade_assessment_dataset_direct_RadiographicPCB)
#Saving grade assessment details for pain mtxIR as an object of lists
grade_details_RadiographicPCB_mtxIR<-list(graph_list,RadiographicPCB_output_df2_list,Grade_assessment_dataset_direct_RadiographicPCB,I_square_list)
#saveRDS(grade_details_HAQ_bio,file="List of HAQ Biologics data components for GRADE_OCT31.rds")
saveRDS(grade_details_RadiographicPCB_mtxIR,file="List of Radiographic Progression Change in Baseline mtxIR data components for GRADE_AUG2026.rds")

RadiographicP_CB_output_one_study_list <- list()
RadiographicP_CB_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

n_studies<-c()
for(i in 1:length(graph_list)){
  df<-RadiographicPCB_output_df2_list[[i]][-c(1,2,countdown(length(RadiographicPCB_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
}

# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
RadiographicP_CB_output_one_study_list  <- RadiographicPCB_output_df2_list[one_idx]
RadiographicP_CB_output_more_than_2_list <- RadiographicPCB_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_RadiographicPCB[
    Grade_assessment_dataset_direct_RadiographicPCB$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_RadiographicPCB[
    Grade_assessment_dataset_direct_RadiographicPCB$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_RadiographicP_CB_mtxIR_one_study <- list(
  graph_list_one_study_list,
  RadiographicP_CB_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_RadiographicP_CB_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  RadiographicP_CB_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_RadiographicP_CB_mtxIR_one_study,
  file = "List of Radiographic Progression Change in Baseline Dmard MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_RadiographicP_CB_mtxIR_more_than_2,
  file = "List of Radiographic Progression Change in Baseline Dmard MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)



#####################################
#####################################
#####################################
#####################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_RadiographicPCB_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct_RadiographicPCB
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_RadiographicPCB$Treatment,
                                                Grade_assessment_dataset_direct_RadiographicPCB$Intevention),
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(RadiographicPCB_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","SMD"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 

glen_direct_evidence_rating_RadioPCB<-read.csv("Glen_2026-09-05 01_20_RADIOCB-MTXIR_complete.csv")
glen_direct_evidence_rating_RadioPCB2<-glen_direct_evidence_rating_RadioPCB[,-c(13:15)]
colnames(glen_direct_evidence_rating_RadioPCB2)[1:2]<-colnames(grade_table4)[1:2]
# #update direct evidence report
# glen_direct_evidence_rating_haq2.1<-cbind(Grade_assessment_dataset_direct_HAQ[,c(1:6)],glen_direct_evidence_rating_haq2[,-c(1,2)])


glen_direct_evidence_rating_RadioPCB2
View(merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_RadioPCB2,
           by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T))
glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_RadioPCB2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"Radiographic Progression"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"SMD"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR_RadiographicPCB_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR_RadiographicPCB_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}

##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for Radiographic Progression MTX-IR CB_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for Radiographic Progression MTX-IR CB_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR RadioPCB OUTCOME
initial_gradeTable_RadioPCB<-Temp_GRADE5

estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_RadioPCB2<-initial_gradeTable_RadioPCB[match(all_nma_estimates,initial_gradeTable_RadioPCB$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()

#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_RadioPCB2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_RadioPCB2.1<-cbind(initial_gradeTable_RadioPCB2[,c(1:14)],first_order_loop_contri,
                                      initial_gradeTable_RadioPCB2[,-c(1:14)])

initial_gradeTable_RadioPCB2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_RadioPCB2.1,
                                                                 all_nma_estimates,
                                                                 estimate_contribution_matrix,
                                                                 glen_direct_evidence_rating6_onlyDirect,
                                                                 top_n = 2))

# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_RadioPCB2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_RadioPCB2.1<-result_df
initial_gradeTable_RadioPCB2.1$combI<-str_replace(initial_gradeTable_RadioPCB2.1$comb, ":", " vs ")

initial_gradeTable_RadioPCB2.1<-process_dataset(initial_gradeTable_RadioPCB2.1, 2)
initial_gradeTable_RadioPCB3<-initial_gradeTable_RadioPCB2.1[match(initial_gradeTable_RadioPCB$comb,initial_gradeTable_RadioPCB2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for Radiographic Progression CB_MTX-IR outcome_2026.csv",row.names = T)
#
###########Assessing Intransitivity
sub_data_intransitivity <- initial_gradeTable_RadioPCB3[initial_gradeTable_RadioPCB3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_RadioPCB3[initial_gradeTable_RadioPCB3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_RadioPCB3.1<-cbind(initial_gradeTable_RadioPCB3[,c(1:23)],
                                      Intransitivity="",
                                      `Certainty of evidence (ignores imprecision)-Indirect`="",
                                      initial_gradeTable_RadioPCB3[,-c(1:23)])
initial_gradeTable_RadioPCB3.1[initial_gradeTable_RadioPCB3.1$Indirect!="-",c("Intransitivity",
                                                                              "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                                    "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_RadioPCB4<-initial_gradeTable_RadioPCB3.1


###Applying the P-values, when no p-values is available because no both direct and indirect evidence available
#P_values_nodesplit_haqmtxIR2
#P_values_nodesplit_radiographicPCBmtxIR2
grade_temp<-cbind(initial_gradeTable_RadioPCB4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_RadioPCB4[,-c(1:25)])
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])


# #######Creating percentage contribution of both direct and indirect evidence
# ####
# colnames_vector<-colnames(estimate_contribution_matrix)
# 
# subradioCBSOF<-grade_temp2[grade_temp2$p_value!="-",]
# subcontriMatradioCB<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subradioCBSOF$comb,]
# 
# 
# # Initializing dirPer and IndirPer vectors
# dirPer <- numeric(nrow(subcontriMatradioCB))
# IndirPer <- numeric(nrow(subcontriMatradioCB))
# 
# # Loop to calculate dirPer and IndirPer
# for (i in 1:nrow(subcontriMatradioCB)) {
#   dirPer[i] <- subcontriMatradioCB[i, colnames(subcontriMatradioCB) == row.names(subcontriMatradioCB)[i]]
#   IndirPer[i] <- 1-dirPer[i]
#   #IndirPer[i] <- sum(subcontriMatHAQ[i, colnames(subcontriMatHAQ) != row.names(subcontriMatHAQ)[i]])
# }
# 
# # Print the results
# cat("dirPer:", dirPer, "\n")
# cat("IndirPer:", IndirPer, "\n")
# 
# 
# grade_temp2[match(row.names(subcontriMatradioCB), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
# grade_temp2[match(row.names(subcontriMatradioCB), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)
# 
# ##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
# grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
# #When only direct or only indirect estimate is present
# grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
#                                  grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
#                                  ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
#                                         grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
# #When both direct and indirect are present
# 
# fg<-c()
# sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
# for(i in 1:nrow(sub_grade_temp2.1)){
#   if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
#     fg[i]<-"Flag for manual review"
#   }else{
#     if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
#       if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
#       ){
#         fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
#       }else{
#         fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
#         
#       }
#     }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
#       fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
#     }else{
#       fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
#     }
#   }
#   
#   
# }
# 
# 
# grade_temp2.1[(grade_temp2.1$initialqoe==" "),"initialqoe"]<-fg


##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))


####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="-",grade_temp2.1[,-c(1:31)])

initial_gradeTable_RadioPCB5<-grade_temp2.2
initial_gradeTable_RadioPCB5<-initial_gradeTable_RadioPCB5[,-c(30,33)]
colnames(initial_gradeTable_RadioPCB5)[c(29,30)]<-c("Initial NMA evidence rating",
                                                    "P-Value")
###After incoherence evaluation
initial_gradeTable_RadioPCB5.1<-cbind(initial_gradeTable_RadioPCB5[,c(1:31)],
                                      qoeAfterIncohenrence="",
                                      initial_gradeTable_RadioPCB5[,-c(1:31)])
initial_gradeTable_RadioPCB5.1[(initial_gradeTable_RadioPCB5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_RadioPCB5.1[(initial_gradeTable_RadioPCB5.1$`P-Value`=="-"),"Initial NMA evidence rating"]


#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_RadioPCB5.1[(initial_gradeTable_RadioPCB5.1$Direct!="-")&(initial_gradeTable_RadioPCB5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_RadioPCB5.1[(initial_gradeTable_RadioPCB5.1$Direct!="-")&(initial_gradeTable_RadioPCB5.1$Indirect=="-"),"Direct"]

###Rating NMA after imprecision
initial_gradeTable_RadioPCB6<-cbind(initial_gradeTable_RadioPCB5.1,
                                    p.estimate=apply(as.data.frame(initial_gradeTable_RadioPCB5.1$NMA),1,function(x){
                                      as.numeric(sub(" .*", "", x))
                                    }),
                                    t(apply(as.data.frame(initial_gradeTable_RadioPCB5.1$NMA),1,function(x){
                                      as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                      
                                    })))

colnames(initial_gradeTable_RadioPCB6)[c(36,37)]<-c("l.bound","u.bound")


###HERE
to_higher_lower_text <- function(e1) {
  sign_val <- sign(e1)
  abs_val <- abs(round(e1, 2))
  text <- ifelse(sign_val < 0, "lower", "higher")
  return(paste(abs_val, text))
}


####Running Continuous NMA model
# load data for MTC
##Getting baseline risk using the placebo intervention
# HAQ_bio_withROB2_indrectness$treatment<-renumber(HAQ_bio_withROB2_indrectness,"treatment")
# HAQ_bio_withROB2_indrectness$study<-renumber(HAQ_bio_withROB2_indrectness,"study")
radiographicP_CB_data<-radioCB_mtxIR_data_new_withROB2_indrectness2
radiographicP_CB_data<-radiographicP_CB_data[order(radiographicP_CB_data$study),c(1:6)]


radiographicP_CB_data_wide2<-reshape(number_of_arms(radiographicP_CB_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(radiographicP_CB_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(radiographicP_CB_data_wide2[,c(3,7,11)],1,xcv)

radiographicP_CB_data_wide3<-cbind(radiographicP_CB_data_wide2,na=numberOfArms)
View(radiographicP_CB_data_wide3)
#finding standard errors
radiographicP_CB_data_wide3$SD.1<-radiographicP_CB_data_wide3$SD.1/sqrt(radiographicP_CB_data_wide3$sampleSize.1)
radiographicP_CB_data_wide3$SD.2<-radiographicP_CB_data_wide3$SD.2/sqrt(radiographicP_CB_data_wide3$sampleSize.2)
radiographicP_CB_data_wide3$SD.3<-radiographicP_CB_data_wide3$SD.3/sqrt(radiographicP_CB_data_wide3$sampleSize.3)
#radiographicP_CB_data_wide3$SD.4<-radiographicP_CB_data_wide3$SD.4/sqrt(radiographicP_CB_data_wide3$sampleSize.4)
colnames(radiographicP_CB_data_wide3)<-c("study", "study_name", "t.1","se.1","n.1", "x.1",
                                         "t.2","se.2","n.2","x.2","t.3","se.3","n.3", "x.3", "na")
                                         #"t.4","se.4","n.4", "x.4","na" )
MTCData_cont <- as.data.frame(radiographicP_CB_data_wide3)
x <- data.matrix(MTCData_cont[,c(6,10,14)])
n <- data.matrix(MTCData_cont[,c(5,9,13)])
t <- data.matrix(apply(MTCData_cont[,c("t.1","t.2","t.3")],c(1,2),as.numeric))
se <- data.matrix(MTCData_cont[,c(4,8,12)])
na <- data.matrix(MTCData_cont[,"na"])
NT <- max(t, na.rm=TRUE)
NS <- nrow(x)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- x[t[,'t.1']==1,"x.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
sePBO <- se[t[,'t.1']==1,"se.1"]

#pooled sd
radiographicP_CB_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")]
calculate_pooled_sd2 <- function(row) {
  sqrt(sum((row[c("sampleSize.1", "sampleSize.2","sampleSize.3")] - 1) * row[c("SD.1", "SD.2","SD.3")]^2,na.rm = TRUE) / sum(row[c("sampleSize.1", "sampleSize.2","sampleSize.3")] - 1,na.rm = TRUE))
}

# Apply the function row-wise
pooled.sd=data.matrix(cbind(apply(radiographicP_CB_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")], 1, calculate_pooled_sd2),
                            apply(radiographicP_CB_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")], 1, calculate_pooled_sd2),
                            apply(radiographicP_CB_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")], 1, calculate_pooled_sd2),
                            apply(radiographicP_CB_data_wide2[,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")], 1, calculate_pooled_sd2)))

data_mean.pooled.sd.ref<-radiographicP_CB_data_wide2[radiographicP_CB_data_wide2$treatment.1==1,c("SD.1","sampleSize.1","SD.2","sampleSize.2","SD.3","sampleSize.3")]
mean.pooled.sd.ref<- calculate_pooled_sd2(data_mean.pooled.sd.ref)


# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,NT-1)), mu=rep(0,NS),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
#model.file4<-"Cont_multiarm_RE.txt"
model.file4<-"tt.txt"#"Cont_new_multiarm_SMD_RE.txt"
n.chains=3
model.inits <- initv1 

cont.model3<-jags.model(file=model.file4,
                        data= list("x"=x,"t"=t,"se"=se,"na" = na[,1], "NT" = NT, "NS" = NS,
                                   "NTpbo"= NTpbo,"xPBO"=xPBO,"sePBO"=sePBO,
                                   "pooled.sd"=pooled.sd,"mean.pooled.sd.ref"=mean.pooled.sd.ref),
                        n.chains=n.chains,
                        inits = model.inits,
                        n.adapt=5000
)

update(object=cont.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("diff","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( cont.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]


#Converting estimates relative to baseline to absolutes
estimates_relative_to_baseline<-initial_gradeTable_RadioPCB6[1:(NT-1),c(1,2,35:37)]
typeof(estimates_relative_to_baseline)

absolutes_from_baseline<-cbind(baseline_risk =round(baseline_risk,2),round(baseline_risk,2)+estimates_relative_to_baseline[,c(3:5)])
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(round(baseline_risk,2),estimates_relative_to_baseline_combine[,7]))
colnames(all_baseline)<-"base"

com.list<-list()
unique_treat2<-unique(initial_gradeTable_RadioPCB6$treat2)
d.sub<-initial_gradeTable_RadioPCB6[initial_gradeTable_RadioPCB6$treat2==unique_treat2[1],c(1,2,35:37)]
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_RadioPCB6[initial_gradeTable_RadioPCB6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,base+d.sub[,c(3:5)])
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe<-com.dataframe[c(1:nrow(initial_gradeTable_RadioPCB6)),]
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9])

higher_lower_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"lower"
  }else{
    text<-"higher"
  }
  return(paste(abs(round(e1-e2,2)),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(higher_lower_text,com.dataframe2[,10],com.dataframe2[,13])))


prefix.assumedrisk <- ' Mean '
prefix.corrrisk<-' Mean '
prefix.md<- " SMD "
suffix.assumedrisk <- ''
suffix.corrrisk <- ''
time<-1


result_com.dataframe_withtext2<-cbind(com.dataframe_withtext[,1:2],
                                      paste(prefix.assumedrisk,com.dataframe_withtext[,6], suffix.assumedrisk, sep=''),
                                      paste(prefix.corrrisk, com.dataframe_withtext[,7], suffix.corrrisk, ' (', com.dataframe_withtext[,8], ' to ',com.dataframe_withtext[,9],')', sep=''),
                                      paste(prefix.md, com.dataframe_withtext[,10], ' (', com.dataframe_withtext[,11], ' to ',com.dataframe_withtext[,12],')', sep=''),
                                      # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                      com.dataframe_withtext[,c(3:5,10:12)]
)

colnames(result_com.dataframe_withtext2)[c(3:11)]<-c("Assumed mean (comparator)",	"Corresponding mean (intervention)",
                                                     "Mean difference","p.estimate","l.bound","u.bound","","","")



initial_gradeTable_RadioPCB7<-cbind(initial_gradeTable_RadioPCB6[,-c(35:37)],
                                    result_com.dataframe_withtext2[,-c(1,2)],pair="")
#colnames(initial_gradeTable_RadioPCB7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_RadioPCB8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_RadioPCB7$treat1,mtxIR_RadiographicPCB_subl)),
                                    `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_RadioPCB7$treat2,mtxIR_RadiographicPCB_subl)),
                                    initial_gradeTable_RadioPCB7)
colnames(initial_gradeTable_RadioPCB8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_RadioPCB8<-initial_gradeTable_RadioPCB8[,-c(46)]
write.csv(initial_gradeTable_RadioPCB8,"Radiographic Progression MTXIR-CB SoF table without imprecision assessment AUG_2026.csv",row.names = F)



effect.size.rating_list<-lapply(split(initial_gradeTable_RadioPCB8, seq_len(nrow(initial_gradeTable_RadioPCB8))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 5, large.effect = 10, framework = T)
})



as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_RadioPCB9<-initial_gradeTable_RadioPCB8
initial_gradeTable_RadioPCB9<-cbind(initial_gradeTable_RadioPCB9[,c(1:34,37:45,35:36)])
initial_gradeTable_RadioPCB9.1<-cbind(initial_gradeTable_RadioPCB9[,c(1:37)],
                                      Small.effect.threshold= 5,
                                      Large.effect.threshold= 10,
                                      Framework = "NULL",
                                      Imprecision.1 = "",
                                      Certainty.of.evidence = "",
                                      initial_gradeTable_RadioPCB9[,c(44,45)])


initial_gradeTable_RadioPCB9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_RadioPCB_mtxIR<-initial_gradeTable_RadioPCB9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_RadioPCB_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_RadioPCB_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_RadioPCB_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_RadioPCB_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_RadioPCB_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_RadioPCB_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_RadioPCB_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_RadioPCB_mtxIR[is.na(detailed_SOF_RadioPCB_mtxIR$Certainty.of.evidence),])
detailed_SOF_RadioPCB_mtxIR[is.na(detailed_SOF_RadioPCB_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

detailed_SOF_RadioPCB_mtxIR$comb<-initial_gradeTable_RadioPCB7$comb

write.csv(detailed_SOF_RadioPCB_mtxIR,"Radiographic Progression MTX-IR-CB_Detail_SoF_AUG_2026.csv", row.names = F)

###############################################################################
################################################################################################
#############Withdrawal (WDAE) Analysis#####################
withdrawal_dataset_mtxIR<-cbind(final_mtx_irDataset_actual[,c(1:51)],  final_mtx_irDataset_actual%>% 
                            select(matches("^Withdrawals.Due.to.Adverse.Events..WDAE\\..*")) %>%
                            select(-contains("Reported"))%>%
                              select(-contains("..4.months"))
                          
)
colnames(withdrawal_dataset_mtxIR)
View(withdrawal_dataset_mtxIR)

library(tibble)
clean_withdrawal <- function(dataset, month) {
  
  if(month == 12){
    
    temp <- dataset %>%
      select(matches("^Withdrawals.Due.to.Adverse.Events..WDAE\\..*")) %>%
      select(contains("12.months"))
    
  }
  
  if(month == 6){
    
    temp <- dataset %>%
      select(matches("^Withdrawals.Due.to.Adverse.Events..WDAE\\..*")) %>%
      select(contains("6.months"))
    
  }
  
  if(month == 3){
    
    temp <- dataset %>%
      select(matches("^Withdrawals.Due.to.Adverse.Events..WDAE\\..*")) %>%
      select(contains("3.months"))
    
  }
  
  if(month == "EndOfTrial"){
    
    temp <- dataset %>%
      select(matches("^Withdrawals.Due.to.Adverse.Events..WDAE\\..*")) %>%
      select(contains("End.of.trial"))
  }
  
  return(temp)
}

View(clean_withdrawal(withdrawal_dataset_mtxIR,3))

clean_withdrawal_all<-function(dataset){
  #3months
  clean_withdrawal_3months<-clean_withdrawal(dataset,3)
  clean_withdrawal_3months<-add_column(clean_withdrawal_3months, "edit" = clean_withdrawal_3months[,colnames(clean_withdrawal_3months)[4]], .after = colnames(clean_withdrawal_3months)[4])
  colnames(clean_withdrawal_3months)[5]<-paste0(colnames(clean_withdrawal_3months)[4],"_edit")
  index3months<-which(is.na(clean_withdrawal_3months[,colnames(clean_withdrawal_3months)[4]]))
  sub_clean_withdrawal_3months<-clean_withdrawal_3months[index3months,]
  clean_withdrawal_3months[index3months,colnames(clean_withdrawal_3months)[5]]<-sub_clean_withdrawal_3months[,3]
  #6months
  clean_withdrawal_6months<-clean_withdrawal(dataset,6)
  clean_withdrawal_6months<-add_column(clean_withdrawal_6months, "edit" = clean_withdrawal_6months[,colnames(clean_withdrawal_6months)[4]], .after = colnames(clean_withdrawal_6months)[4])
  colnames(clean_withdrawal_6months)[5]<-paste0(colnames(clean_withdrawal_6months)[4],"_edit")
  index6months<-which(is.na(clean_withdrawal_6months[,colnames(clean_withdrawal_6months)[4]]))
  sub_clean_withdrawal_6months<-clean_withdrawal_6months[index6months,]
  clean_withdrawal_6months[index6months,colnames(clean_withdrawal_6months)[5]]<-sub_clean_withdrawal_6months[,3]
  #12months
  clean_withdrawal_12months<-clean_withdrawal(dataset,12)
  clean_withdrawal_12months<-add_column(clean_withdrawal_12months, "edit" = clean_withdrawal_12months[,colnames(clean_withdrawal_12months)[4]], .after = colnames(clean_withdrawal_12months)[4])
  colnames(clean_withdrawal_12months)[5]<-paste0(colnames(clean_withdrawal_12months)[4],"_edit")
  index12months<-which(is.na(clean_withdrawal_12months[,colnames(clean_withdrawal_12months)[4]]))
  sub_clean_withdrawal_12months<-clean_withdrawal_12months[index12months,]
  clean_withdrawal_12months[index12months,colnames(clean_withdrawal_12months)[5]]<-sub_clean_withdrawal_12months[,3]
  #EndOfTrial
  clean_withdrawal_EndOfTrial<-clean_withdrawal(dataset,"EndOfTrial")
  clean_withdrawal_EndOfTrial<-add_column(clean_withdrawal_EndOfTrial, "edit" = clean_withdrawal_EndOfTrial[,colnames(clean_withdrawal_EndOfTrial)[4]], .after = colnames(clean_withdrawal_EndOfTrial)[4])
  colnames(clean_withdrawal_EndOfTrial)[5]<-paste0(colnames(clean_withdrawal_EndOfTrial)[4],"_edit")
  indexEndOfTrial<-which(is.na(clean_withdrawal_EndOfTrial[,colnames(clean_withdrawal_EndOfTrial)[4]]))
  sub_clean_withdrawal_EndOfTrial<-clean_withdrawal_EndOfTrial[indexEndOfTrial,]
  clean_withdrawal_EndOfTrial[indexEndOfTrial,colnames(clean_withdrawal_EndOfTrial)[5]]<-sub_clean_withdrawal_EndOfTrial[,3]
  
  
  return(cbind(clean_withdrawal_3months,clean_withdrawal_6months,clean_withdrawal_12months,clean_withdrawal_EndOfTrial))
}


View(clean_withdrawal_all(withdrawal_dataset_mtxIR))

###Get new withdrawal dataset
withdrawal_dataset_mtxIR1<-cbind(final_mtx_irDataset_actual[,c(1:51)],clean_withdrawal_all(withdrawal_dataset_mtxIR))
View(withdrawal_dataset_mtxIR1)
withdrawal_dataset_mtxIR1[withdrawal_dataset_mtxIR1=="NA"]<-NA

withdrawal_dataset_mtxIR2<- withdrawal_dataset_mtxIR1

# Columns by time point
wdae_3m_cols <- grep(
  "^Withdrawals\\.Due\\.to\\.Adverse\\.Events\\.\\.WDAE\\.\\.3\\.months",
  names(withdrawal_dataset_mtxIR2),
  value = TRUE
)

wdae_6m_cols <- grep(
  "^Withdrawals\\.Due\\.to\\.Adverse\\.Events\\.\\.WDAE\\.\\.6\\.months",
  names(withdrawal_dataset_mtxIR2),
  value = TRUE
)

wdae_12m_cols <- grep(
  "^Withdrawals\\.Due\\.to\\.Adverse\\.Events\\.\\.WDAE\\.\\.12\\.months",
  names(withdrawal_dataset_mtxIR2),
  value = TRUE
)

wdae_eot_cols <- grep(
  "^Withdrawals\\.Due\\.to\\.Adverse\\.Events\\.\\.WDAE\\.\\.End\\.of\\.trial",
  names(withdrawal_dataset_mtxIR2),
  value = TRUE
)

# End of trial
withdrawal_dataset_mtxIR2[
  withdrawal_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,6,9,11,12,13,18),
  wdae_eot_cols
] <- NA

# 12 months
withdrawal_dataset_mtxIR2[
  withdrawal_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,23,24,63),
  wdae_12m_cols
] <- NA

# 6 months
withdrawal_dataset_mtxIR2[
  withdrawal_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,9,11,12,13,18,23,24,63),
  wdae_6m_cols
] <- NA

# 3 months
withdrawal_dataset_mtxIR2[
  withdrawal_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(6,9,11,12,13,18,23,24,63),
  wdae_3m_cols
] <- NA


write.csv(withdrawal_dataset_mtxIR2,"Withdrawal dataset for MTX-IR_AUG2026.csv",row.names = F)

#Temporary reading of data
withdrawal_dataset_mtxIR2<-as.data.frame(read_excel('/Users/mujaab/Downloads/Withdrawal dataset for MTX-IR_AUG2026_withExposure.xlsx'))

##Kremer 2006, no exposure data
withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Kremer 2006",wdae_12m_cols]<-NA
##Kremer 2005, no exposure data
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Kremer 2005",wdae_6m_cols]<-NA

##Genovese 2011
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Genovese 2011","Withdrawals.Due.to.Adverse.Events..WDAE..6.months.Exposure..patient.years."]<-c(((736*24)+(17*(24/2)))/52,((721*24)+(25*(24/2)))/52)

##Choy 2012
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Choy 2012","Withdrawals.Due.to.Adverse.Events..WDAE..6.months.Exposure..patient.years."]<-c(((119*24)+(6*(24/2)))/52,((124*24)+(7*(24/2)))/52)

#Weinblatt 2003, no exposure data
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Weinblatt 2003",wdae_6m_cols]<-NA

#Chen 2009, no exposure data
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Chen 2009",wdae_3m_cols]<-NA

#Keystone 2008, no exposure data
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Keystone 2008",wdae_12m_cols]<-NA

##Kameda 2011
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Kameda 2011","Withdrawals.Due.to.Adverse.Events..WDAE..12.months.Exposure..patient.years."]<-c(((69*52)+(4*(52/2)))/52,((73*52)+(1*(52/2)))/52)

#Weinblatt 1999
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Weinblatt 1999",wdae_6m_cols]<-NA

#Lipsky 2000
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Lipsky 2000","Withdrawals.Due.to.Adverse.Events..WDAE..12.months.Exposure..patient.years."]<-c(((88*52)+(7*(52/2)))/52,((86*52)+(5*(52/2)))/52)

#Tanaka 2012
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Tanaka 2012","Withdrawals.Due.to.Adverse.Events..WDAE..6.months.Exposure..patient.years."]<-c(((88*24)+(1*(24/2)))/52,((86*24)+(4*(24/2)))/52)

#Kay 2008
withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Kay 2008",wdae_12m_cols]<-NA
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Kay 2008","Withdrawals.Due.to.Adverse.Events..WDAE..6.months.Exposure..patient.years."]<-c(((34*24)+(3*(24/2)))/52,((37*24)+(2*(24/2)))/52)

#Zhang 2006
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Zhang 2006",wdae_3m_cols]<-NA

#Westhovens 2006
#withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Westhovens 2006","Withdrawals.Due.to.Adverse.Events..WDAE..6.months.Exposure..patient.years."]<-c(((363*24)+(8*(22/2)))/52,((360*24)+(18*(22/2)))/52)
#Emery 2010
withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Emery 2010",wdae_12m_cols]<-NA
#Schiff 2008
withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Schiff 2008",wdae_12m_cols]<-NA
#Maini 2006
withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study.Identifier=="Maini 2006",wdae_eot_cols]<-NA


colindex<-as.vector(which(colSums(!is.na(withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==16,c(1:3,52:75)])) > 0))
temp_withdrawal<-withdrawal_dataset_mtxIR2[,c(1:3,52:75)]
temp_withdrawal2<-temp_withdrawal[temp_withdrawal$Study_no==16,colindex]
temp_withdrawal2
temp_withdrawal2[,c(1:3,6,4,7)]
colnames(temp_withdrawal2)
temp_withdrawal2[,c(1:3,12,13,15)]
temp_withdrawal2[,c(1:3,14,10,11)]

withdrawal_dataset_mtxIR2[52:75]<-lapply(withdrawal_dataset_mtxIR2[52:75], as.numeric)


empty_withdrawalStudy<-list()  
keep_withdrawalStudy_for_errorCheck<-list()
withdrawalStudy_for_analysis<-list()

for(i in 1:length(unique(withdrawal_dataset_mtxIR2$Study_no))){
  if(ncol(withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,c(52:75)][colSums(!is.na(withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,c(52:75)])) > 0])==0){
    empty_withdrawalStudy[[i]]<-withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,]
  }else if(!(any(c(1,7,13,19) %in% as.vector(which(colSums(!is.na(withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,c(52:75)])) > 0))))){
    keep_withdrawalStudy_for_errorCheck[[i]]<-withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,]
  }else{
    colindex<-as.vector(which(colSums(!is.na(withdrawal_dataset_mtxIR2[withdrawal_dataset_mtxIR2$Study_no==i,c(1:3,52:75)])) > 0))
    temp_withdrawal<-withdrawal_dataset_mtxIR2[,c(1:3,52:75)]
    temp_withdrawal2<-temp_withdrawal[temp_withdrawal$Study_no==i,colindex]
    df_wdae_lastTP <- temp_withdrawal2 %>%
      select(matches("WDAE.*(Exposure|N|events_edit)")) %>%
      {
        cols <- names(.)
        
        suffix <- str_extract(cols, "(Exposure|N|events_edit)")
        timept <- str_extract(cols, "(\\d+\\.months|End\\.of\\.trial)")
        
        time_order   <- c("3.months","6.months","12.months","End.of.trial")
        suffix_order <- c("Exposure","N","events_edit")
        
        # order columns correctly first
        ord <- order(match(timept, time_order),
                     match(suffix, suffix_order))
        
        df <- .[, ord]
        
        # ---- NEW PART: keep only LAST available timepoint ----
        ordered_time_present <- time_order[time_order %in% timept]
        last_time <- tail(ordered_time_present, 1)
        
        df[, str_detect(names(df), last_time)]
      } #%>%dplyr::select(-dplyr::matches("Total\\.number\\.of\\.events$"))
    
    temp_withdrawal3 <- cbind(temp_withdrawal2[,c(1:3)],df_wdae_lastTP[c(1:3)], colnames(temp_withdrawal2)[4])
    colnames(temp_withdrawal3)<-c("Study_no","Study.Identifier","Intervention","exposure","sampleSize","responders","month")
    withdrawalStudy_for_analysis[[i]]<-temp_withdrawal3
  }
  
}

View(do.call(rbind,withdrawalStudy_for_analysis))
View(do.call(rbind,empty_withdrawalStudy))
current_withdrawal_dataset<-do.call(rbind,withdrawalStudy_for_analysis)
##Original needed to compare after checking for errors
# current_withdrawal_dataset_original<-distinct(current_withdrawal_dataset)
# current_withdrawal_dataset<-distinct(current_withdrawal_dataset)
#current_withdrawal_dataset<-current_withdrawal_dataset[current_withdrawal_dataset$Study.Identifier %in%actual_included_studies,]

current_withdrawal_dataset$Study_no<-renumber(current_withdrawal_dataset,"Study_no")
row.names(current_withdrawal_dataset)<-NULL

current_withdrawal_dataset$exposure<-round(as.numeric(current_withdrawal_dataset$exposure),2)
current_withdrawal_dataset$sampleSize<-as.numeric(current_withdrawal_dataset$sampleSize)
current_withdrawal_dataset$responders<-as.numeric(current_withdrawal_dataset$responders)

#Fleischmann 2019
current_withdrawal_dataset[current_withdrawal_dataset$Study.Identifier=="Fleischmann 2019","sampleSize"][3]<-327

#Kim 2000
current_withdrawal_dataset[current_withdrawal_dataset$Study.Identifier=="Kim 2000","responders"][2]<-0


#Drop rows with NA
current_withdrawal_dataset<-current_withdrawal_dataset%>%
  filter(!is.na(exposure))

#####Formatting the month column
current_withdrawal_dataset$month<- stringr::str_extract(current_withdrawal_dataset$month, "\\d+(?=\\.months)") %>%
  as.numeric()
current_withdrawal_dataset$month<- ifelse(is.na(current_withdrawal_dataset$month), "EndOfTrial",current_withdrawal_dataset$month)

##Mean and Median of time frame
mean(current_withdrawal_dataset$month)
median(current_withdrawal_dataset$month)

View(current_withdrawal_dataset)
current_withdrawal_dataset$Study_no<-renumber(current_withdrawal_dataset,"Study_no")
current_withdrawal_dataset_mtxIR<-current_withdrawal_dataset
current_withdrawal_dataset_mtxIR <- current_withdrawal_dataset_mtxIR %>%
  dplyr::rename("treatment" = Intervention)

##Drop the high risk of bias studies
current_withdrawal_dataset_mtxIR <- current_withdrawal_dataset_mtxIR%>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012","Moreland 2002","Nishimoto 2007"))

current_withdrawal_dataset_mtxIR$Study_no<-renumber(current_withdrawal_dataset_mtxIR,"Study_no")

library(dplyr)
for(i in 1:nrow(current_withdrawal_dataset_mtxIR)){
  ppp<- master.medtable %>% filter(current_withdrawal_dataset_mtxIR$treatment[i] ==master.medtable$medname )
  current_withdrawal_dataset_mtxIR$treatment[i]<-ppp[1,1]
}



View(do.call(rbind,empty_withdrawalStudy))
View(do.call(rbind,keep_withdrawalStudy_for_errorCheck))

##Creating initial legend
mtxIR.withdrawal.medlegend <- new.medtable(current_withdrawal_dataset_mtxIR,'treatment')
mtxIR.withdrawal.medlegend


###Adding ROB data
ffff<-unique(current_withdrawal_dataset_mtxIR$Study.Identifier)
data_qq<-current_withdrawal_dataset_mtxIR
View(ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff,])
current_assessment_for_studies<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies<- current_assessment_for_studies %>%
  select(c(col_names))

list_appraisal<-list()
for(i in 1:length(ffff)){
  temp11<-data_qq[data_qq$Study.Identifier==ffff[i],]
  temp12<-current_assessment_for_studies[current_assessment_for_studies$Study.Id==ffff[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal[[i]]<-temp13
}
current_WDAE_dataset_mtxIR_withROB<-do.call(rbind,list_appraisal)
View(current_WDAE_dataset_mtxIR_withROB)
current_WDAE_dataset_mtxIR_withROB2<-current_WDAE_dataset_mtxIR_withROB%>%
  select(c(Study_no,Study.Identifier, "treatment", "exposure", "sampleSize",
           "responders","month","OVERALL..WITHDRAWALS"))
View(current_WDAE_dataset_mtxIR_withROB2)

#Drop further high risk bias
current_WDAE_dataset_mtxIR_withROB2 <- current_WDAE_dataset_mtxIR_withROB2%>%
  filter(!OVERALL..WITHDRAWALS =="high")

##check Cohen 2020 for bias data
current_WDAE_dataset_mtxIR_withROB2[current_WDAE_dataset_mtxIR_withROB2$Study.Identifier=="Cohen 2020","OVERALL..WITHDRAWALS"]<-"low"


###Adding Indirectness data
study_unique<-unique(current_WDAE_dataset_mtxIR_withROB2$Study.Identifier)
data_nn<-current_WDAE_dataset_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")



resolved_indirectness_consensus<-resolved_indirectness_consensus_temp%>%
  dplyr::rename( "Consensus"= `Final Rating`)

View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%study_unique,]



list_appraisal2<-list()

for(i in 1:length(study_unique)){
  temp11<-data_nn[data_nn$Study.Identifier==study_unique[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==study_unique[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(nrow(temp12)==0){
    temp13<-cbind(temp11,Consensus="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
current_WDAE_dataset_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal2)
#colnames(current_WDAE_dataset_mtxIR_withROB2_indrectness)[9]<-"Indirectness"
current_WDAE_dataset_mtxIR_withROB2_indrectness<- current_WDAE_dataset_mtxIR_withROB2_indrectness%>%
  dplyr::rename("Indirectness"= Consensus, "study"= Study_no, "study_name" = Study.Identifier)
#ACR50_try_data_new_withROB2_indrectness<-cbind(ACR50_try_data_new_withROB2,Indirectness="low")
current_WDAE_dataset_mtxIR_withROB2_indrectness$Indirectness<-tolower(current_WDAE_dataset_mtxIR_withROB2_indrectness$Indirectness)
View(current_WDAE_dataset_mtxIR_withROB2_indrectness)


current_WDAE_dataset_mtxIR_withROB2_indrectness<- current_WDAE_dataset_mtxIR_withROB2_indrectness[order(current_WDAE_dataset_mtxIR_withROB2_indrectness$study),]
row.names(current_WDAE_dataset_mtxIR_withROB2_indrectness)<-NULL



###removing the zero cells
#current_WDAE_dataset_mtxIR_withROB2_indrectness<-current_WDAE_dataset_mtxIR_withROB2_indrectness[!current_WDAE_dataset_mtxIR_withROB2_indrectness$study_name%in%c("Ma 2024", "Verschueren 2015","Marcora 2006"),]
mtxIR.withdrawal.medlegend2<-new.medtable(current_WDAE_dataset_mtxIR_withROB2_indrectness,'treatment')
##NMA plot reason for narrowing of network
nma.plot(current_WDAE_dataset_mtxIR_withROB2_indrectness[,1],current_WDAE_dataset_mtxIR_withROB2_indrectness[,3],data = current_WDAE_dataset_mtxIR_withROB2_indrectness, trtname = as.vector(mtxIR.withdrawal.medlegend2$medname))


################################################################################################################
#trying the narrowing a network function
#############################################################################################
current_WDAE_dataset_mtxIR_withROB2_indrectness
current_WDAE_dataset_mtxIR_withROB2_indrectness<-current_WDAE_dataset_mtxIR_withROB2_indrectness[!current_WDAE_dataset_mtxIR_withROB2_indrectness$study_name%in%c("Conaghan 2013","Tanaka 2016","Chen 2011","MacIsaac 2014"),]

# current_WDAE_dataset_mtxIR_withROB2_indrectness3<-current_WDAE_dataset_mtxIR_withROB2_indrectness3%>%
#   dplyr::rename("study"=Study_no,"study_name"=Study.Identifier)
#Narrow 1
mtxIR_withdrawal_narrow<-narrow.all(current_WDAE_dataset_mtxIR_withROB2_indrectness)
mtxIR_withdrawal_narrow


#sub network plot 1
mtxIR.withdrawal.d<-mtxIR_withdrawal_narrow[[1]]
length(unique(mtxIR.withdrawal.d$treatment))
mtxIR_withdrawal_subl<-new.medtable(mtxIR.withdrawal.d,'treatment')
mtxIR.withdrawal.d$treatment<-renumber(mtxIR.withdrawal.d,"treatment")
#Data to download
mtxIR.withdrawal.d_todownload<- mtxIR.withdrawal.d
write.csv(mtxIR.withdrawal.d_todownload,"WDAE MTX-IR dataset_AUG2026.csv",row.names = F)

# 1. Open jpeg file
jpeg("Withdrawal Due to Adverse Events (WDAE) MTX-IR Network plot with named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.withdrawal.d[,1],mtxIR.withdrawal.d[,3],data = mtxIR.withdrawal.d, trtname = as.vector(mtxIR_withdrawal_subl$medname), title = "Withdrawal due to adverse events (WDAE) Network plot-MTX-IR")

# Close the pdf file
dev.off() 
# 2. Open jpeg file
jpeg("Withdrawal Due to Adverse Events (WDAE) MTX-IR Network plot without named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.withdrawal.d[,1],mtxIR.withdrawal.d[,3],data = mtxIR.withdrawal.d, title = "Withdrawal due to adverse events (WDAE) Network plot-MTX-IR")

# Close the pdf file
dev.off() 
#############################################################################
#Performing sub network plot 3
library(coda)
library(gemtc)
View(mtxIR.withdrawal.d)
mtxIR_withdrawal_subl
#Correcting for zero cells
mtxIR.withdrawal.d3<-cbind(mtxIR.withdrawal.d,
                             responders_corr_zero_cells=ifelse(mtxIR.withdrawal.d$responders==0,0.5,
                                                               mtxIR.withdrawal.d$responders))
##Sorting to start with 1 as reference treatment
##Function
reorder_reference<-function(data,Study_no){
  ss<- data[data$study==Study_no,]
  return(ss[order(ss$treatment), ])
}


mtxIR.withdrawal.d3$study<-renumber(mtxIR.withdrawal.d3,"study")
for(i in 1:length(unique(mtxIR.withdrawal.d3$study))){
  mtxIR.withdrawal.d3[mtxIR.withdrawal.d3$study==i,]<-reorder_reference(mtxIR.withdrawal.d3,i)
}
row.names(mtxIR.withdrawal.d3)<-NULL

data_withdrawal_mtxIR_1<-mtxIR.withdrawal.d3
#data_withdrawal_mtxIR_1[data_withdrawal_mtxIR_1$study_name=="Clark 1993","exposure"][2]<-((65*6)+(0*(6/2)))/12
#colnames(data_withdrawal_mtxIR_1)[c(6,8)]<-c("responders_corr_zero_cells","responders" )
data_withdrawal_mtxIR_1$treatment<-renumber(data_withdrawal_mtxIR_1,"treatment")
data_withdrawal_mtxIR_1$study<-renumber(data_withdrawal_mtxIR_1,"study")
nma.plot(data_withdrawal_mtxIR_1$study,data_withdrawal_mtxIR_1$treatment,data = data_withdrawal_mtxIR_1, trtname = as.vector(mtxIR_withdrawal_subl$medname), title = "Withdrawal due to adverse events (WDAE) Network plot-MTX-IR")
network_data_withdrawal_mtxIR_1<-mtc.network(data_withdrawal_mtxIR_1)

##RUN the model
set.seed(51221)
#mtc.hy.prior("prec", "dgamma", 0.01, 0.01)
#mtc.hy.prior("std.dev", "dunif", 0, 2)
network_data_withdrawal_mtxIR_1.model <- mtc.model(network_data_withdrawal_mtxIR_1, likelihood='poisson', link='log',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
#network_data_withdrawal_mtxIR_1.model <- mtc.model(network_data_withdrawal_mtxIR_1, likelihood='poisson', link='log',hy.prior=mtc.hy.prior("prec", "dgamma", 0.1, 0.1))
#network_data_withdrawal_mtxIR_1.model <- mtc.model(network_data_withdrawal_mtxIR_1, likelihood='poisson', link='log',hy.prior=mtc.hy.prior("std.dev", "dunif", 0, 2))


#By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_withdrawal_mtxIR_1.model$inits <- mapply(c, network_data_withdrawal_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_withdrawal<-mtc.run(network_data_withdrawal_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_withdrawal)
View(as.data.frame(relative.effect.table.OR(results_mtxIR_withdrawal)))
refHR_mtxIR_withdrawal<-as.data.frame(relative.effect.table.OR2(results_mtxIR_withdrawal,comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refHR_mtxIR_withdrawal)<-mapvalues(row.names(refHR_mtxIR_withdrawal),
                                               from = c(1:(nrow(mtxIR_withdrawal_subl))),
                                               to = as.vector(mtxIR_withdrawal_subl$medname))
colnames(refHR_mtxIR_withdrawal)<-mapvalues(colnames(refHR_mtxIR_withdrawal),
                                              from = c(1:(nrow(mtxIR_withdrawal_subl))),
                                              to = as.vector(mtxIR_withdrawal_subl$medname))
View(refHR_mtxIR_withdrawal)

index_diagonal3<-data.frame(index1=1:nrow(refHR_mtxIR_withdrawal),index2=1:nrow(refHR_mtxIR_withdrawal))
options(warn=-1)
for(i in 1:nrow(refHR_mtxIR_withdrawal)){
  refHR_mtxIR_withdrawal[index_diagonal3$index1[i],index_diagonal3$index2[i]]<-selected_treatment_name(as.numeric(refHR_mtxIR_withdrawal[row(refHR_mtxIR_withdrawal)==col(refHR_mtxIR_withdrawal)])[i],mtxIR_withdrawal_subl)
}
options(warn=0)
View(refHR_mtxIR_withdrawal)
#write.xlsx(refHR_mtxIR_withdrawal,"Results for mtxIR_2026.xlsx",sheetName = "WDAE",row.names = T, append = T)
#write.xlsx(refHR_bio_withdrawal,"Results for Biologics Informative Prior.xlsx",sheetName = "WDAE",row.names = T, append = T)


##########
write.csv(refHR_mtxIR_withdrawal,"WDAE MTX-IR pairwise comparisons_2026.csv", row.names = T)
#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_OCT2026.xlsx")

addWorksheet(wb, "WDAE")
writeData(wb, "WDAE", refHR_mtxIR_withdrawal, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_OCT2026.xlsx", overwrite = TRUE)

##GRADE Assessment
nma_estimate_mtxIR_WDAE_dataframe<-nma.estimate.table.OR(results_mtxIR_withdrawal)
View(nma_estimate_mtxIR_WDAE_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_WDAE_mtxIR_1.ns<-mtc.nodesplit(network_data_withdrawal_mtxIR_1, thin=50,likelihood='poisson', link='log',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 
names(network_data_WDAE_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_WDAE_mtxIR_1.ns)

# Overall summary and plot
summary.network_data_WDAE_mtxIR_1.ns <- summary(network_data_WDAE_mtxIR_1.ns)


#Direct effects
direct.data<-summary.network_data_WDAE_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_WDAE_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_WDAE_mtxIR_1.ns$cons.effect

#P-values
P_values_nodesplit_wdaemtxIR<-summary.network_data_WDAE_mtxIR_1.ns$p.value
P_values_nodesplit_wdaemtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_wdaemtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_wdaemtxIR2<-P_values_nodesplit_wdaemtxIR
P_values_nodesplit_wdaemtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_wdaemtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_wdaemtxIR2$comb<-paste0(P_values_nodesplit_wdaemtxIR2$t1,":",P_values_nodesplit_wdaemtxIR2$t2)



direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(exp(direct.data[,c(3:5)]),2),round(exp(indirect.data[,c(3:5)]),2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","RR1","low1","high1","RR2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                round(1/direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","RR1","low1","high1","RR2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)

grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_WDAE_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)

current_WDAE_dataset_mtxIR_withROB2_indrectness3<-mtxIR.withdrawal.d
current_WDAE_dataset_mtxIR_withROB2_indrectness3$treatment<-renumber(current_WDAE_dataset_mtxIR_withROB2_indrectness3,"treatment")
current_WDAE_dataset_mtxIR_withROB2_indrectness3$OVERALL..WITHDRAWALS<-ifelse(current_WDAE_dataset_mtxIR_withROB2_indrectness3$OVERALL..WITHDRAWALS=="unclear","some concerns",current_WDAE_dataset_mtxIR_withROB2_indrectness3$OVERALL..WITHDRAWALS)


pairwise_WDAE_RR<-pairwise(treat=treatment,event = responders_corr_zero_cells,data=network_data_withdrawal_mtxIR_1$data.ab,
                           studlab = study,n=sampleSize,sm="RR")


netPairwiseWDAE<-netmeta(pairwise_WDAE_RR, sm="RR")
netPairwiseWDAE_split<-netsplit(netPairwiseWDAE)
View(netPairwiseWDAE_split$direct.common)
dat.net.contrib<-netcontrib(netPairwiseWDAE, method = "r")
View(dat.net.contrib$random)

direct_evidence_new<-cbind(netPairwiseWDAE_split$direct.common[,c(1:5)],netPairwiseWDAE_split$indirect.common[,c(2:5)])
View(direct_evidence_new)
#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")
#direct_only4%in%paste(nma_estimate_bio_HAQ_dataframe$treat1,nma_estimate_bio_HAQ_dataframe$treat2,sep = ":")

#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]

##Getting direct evidence for intervention in only multi-arm trials
zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
colnames(zzz)<-c("treat1", "treat2")
grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]





grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_withdrawal_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_withdrawal_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))


#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)

graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
WDAE_output_df2_list<-list()
I_square_list<-list()
all_direct_evidence2[,c(1,2)]

WDAE_mtxIR_data_new_withROB2_indrectness2<-current_WDAE_dataset_mtxIR_withROB2_indrectness3

WDAE_mtxIR_data_new_withROB2_indrectness2<- WDAE_mtxIR_data_new_withROB2_indrectness2%>%
  dplyr::rename("ROB" = OVERALL..WITHDRAWALS)

for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-WDAE_mtxIR_data_new_withROB2_indrectness2[WDAE_mtxIR_data_new_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select("study", "exposure.1", "responders.1","exposure.2", "responders.2" )
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","E.t","x.t","E.c","x.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    # model.file<-"Rate_direct_onestudy.txt"
    # n.chains=3
    # set.seed(2453)
    # seeds <- sample.int(n.chains, n = .Machine$integer.max)
    # model.inits <- list(
    #   list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
    #   list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
    #   list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    # jj.model<-jags.model(file=model.file,
    #                      data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
    #                                "x.c"=temp12$x.c,
    #                                "E.c"=temp12$E.c),
    #                      n.chains=n.chains,
    #                      inits = model.inits,
    #                      n.adapt=500
    # )
    # update(object=jj.model,n.iter=10000,progress.bar='text')
    # ## specify parameters to be monitored
    # params <- c("RR")
    # ## run JAGS and save posterior samples
    # samps <- coda.samples( jj.model, params, n.iter=10000 )
    model.file<-"Rate_direct_multiplestudy_FE.txt"
    n.chains=3
    set.seed(2553)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
                                   "x.c"=temp12$x.c,
                                   "E.c"=temp12$E.c,
                                   "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=1000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("RR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=1000)
    summary_entry1<-summary(samps)
    if(round(summary_entry1$statistics[1],1)==round(summary_entry1$quantiles[3],1)){
      fgh[[k]]<-data.frame(RR=summary_entry1$statistics[1],
                           Lower=summary_entry1$quantiles[c(1)],
                           Upper=summary_entry1$quantiles[c(5)])
    }else{
      fgh[[k]]<-data.frame(RR=summary_entry1$quantiles[3],
                           Lower=summary_entry1$quantiles[c(1)],
                           Upper=summary_entry1$quantiles[c(5)])
    }
    
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    if(round(summary_entry1$statistics[1],1)==round(summary_entry1$quantiles[3],1)){
      pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
      
    }else{
      pooled_estimate<-c(summary_entry1$quantiles[3],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
      
    }  }else{
      temp12<-temp11.1
      model.file<-"Rate_direct_multiplestudy_RE.txt"
      n.chains=3
      set.seed(2553)
      seeds <- sample.int(n.chains, n = .Machine$integer.max)
      model.inits <- list(
        list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
        list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
        list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
      jj.model<-jags.model(file=model.file,
                           data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
                                     "x.c"=temp12$x.c,
                                     "E.c"=temp12$E.c,
                                     "NS"=length(temp12$NS)),
                           n.chains=n.chains,
                           inits = model.inits,
                           n.adapt=5000
      )
      update(object=jj.model,n.iter=10000,progress.bar='text')
      ## specify parameters to be monitored
      params <- c("RR")
      ## run JAGS and save posterior samples
      samps <- coda.samples( jj.model, params, n.iter=10000 )
      summary_entry<-summary(samps)
      pooled_estimate<-c(summary_entry$quantiles[3],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
    }
  
  
  #Frequentist approach for I^2
  nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select("study", "exposure.1", "responders.1","exposure.2", "responders.2" )
  nma_data50_4<-nma_data50_3
  colnames(nma_data50_4)<-c("NS","py.t","x.t","py.c","x.c")
  row.names(nma_data50_4)<-NULL
  m.inc<- metainc(event.e =x.t,
                  event.c = x.c,
                  time.e = py.t,
                  time.c = py.c,
                  data = nma_data50_4,
                  sm = "IRR",
                  method = "Inverse",
                  MH.exact = TRUE,
                  fixed = FALSE,
                  random = TRUE,
                  method.tau = "PM",
                  hakn = TRUE,
                  title = "WDAE"
  )
  
  
  summary_m.inc<-summary(m.inc)
  contribution_weight<-as.character(round((summary_m.inc$w.random/sum(summary_m.inc$w.random))*100,2))
  pooled_estimate2<-data.frame(RR=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.inc$tau2,4),
                          Lower=round(summary_m.inc$lower.tau2,4),
                          Upper=round(summary_m.inc$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.inc$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.inc$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.inc$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.inc$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  
  
  
  
  RR_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1)
  base_data <- tibble(mean  = RR_individual_studies$RR, 
                      lower = RR_individual_studies$Lower,
                      upper = RR_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      # treatment1 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      # treatment2 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      # treatment1 = paste0(as.character(temp11.2$E.c),"          ","(",as.character(temp11.2$x.c),")"),
                      # treatment2 = paste0(as.character(temp11.2$E.t),"          ","(",as.character(temp11.2$x.t),")"),
                      treatment1 = paste0(as.character(temp11.2$x.t),"          ","(",as.character(temp11.2$E.t),")"),
                      treatment2 = paste0(as.character(temp11.2$x.c),"          ","(",as.character(temp11.2$E.c),")"),
                      RR = paste(format(round(RR_individual_studies$RR,2),nsmall=2),paste("(",format(round(RR_individual_studies$Lower,2),nsmall=2),", ",format(round(RR_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$RR, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   #                   RR = paste(round(as.numeric(pooled_estimate2$RR),2),paste("(",round(as.numeric(pooled_estimate2$Lower),2),", ",round(as.numeric(pooled_estimate2$Upper),2),")",sep="")),
                   RR = paste(format(round(pooled_estimate2$RR,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("WDAE", "Study"),
                   # treatment1 = c("",selected_treatment_name(unique(temp11$treatment.1),subl)),
                   # treatment2 = c("",selected_treatment_name(unique(temp11$treatment.2),subl)),
                   treatment1 = c("Events (exposure patient-years)",selected_treatment_name(unique(temp11$treatment.1),mtxIR_withdrawal_subl)),
                   treatment2 = c("Events (exposure patient-years)",selected_treatment_name(unique(temp11$treatment.2),mtxIR_withdrawal_subl)),
                   #ROB=c(paste("Favours",selected_treatment_name(unique(temp11$treatment.1),subl)),"ROB"),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   # OR = c(paste("Favours",selected_treatment_name(unique(temp11$treatment.2),subl)), "OR"),
                   RR = c("", "RR"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  WDAE_output_df <- bind_rows(header,
                              base_data,
                              empty_row,
                              pooled)
  
  
  
  
  WDAE_output_df2<-rbind(WDAE_output_df,WDAE_output_df[c(nrow(WDAE_output_df)),])
  WDAE_output_df2[nrow(WDAE_output_df),1]<-I2
  WDAE_output_df2[nrow(WDAE_output_df),2]<-NA
  WDAE_output_df2[nrow(WDAE_output_df),c(6:10)]<-c(NA,NA,NA,NA,NA)
  
  WDAE_output_df2_list[[i]]<-WDAE_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-WDAE_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, RR,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(0.1, 2.4), 
               xlog = T, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_withdrawal_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_withdrawal_subl))
    )
  
}

length(graph_list)
graph_list[[1]]
n_studies<-c()
n_participants<-c()

countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}


View(WDAE_output_df2_list[[19]])

for(i in 1:length(graph_list)){
  df<-WDAE_output_df2_list[[i]][-c(1,2,countdown(length(WDAE_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df, 1, function(x) gsub("\\s*\\(.*", "", x)))))
}


Grade_assessment_dataset_direct_wdae<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2[,1])),mtxIR_withdrawal_subl),
                                                 selected_treatment_name(as.numeric(as.vector(all_direct_evidence2[,2])),mtxIR_withdrawal_subl),
                                                 Outcome="WDAE",
                                                 `Effect measure`="RR",
                                                 `Number of studies`=n_studies,
                                                 `Number of participant`=n_participants,
                                                 ROB="",Indirectness="",Imprecision="",
                                                 Inconsistency="",
                                                 Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_wdae)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                        "Intervention",
                                                                        "Effect measure",
                                                                        "Number of studies",
                                                                        "Number of participant",
                                                                        "Risk of Bias",
                                                                        "Publication Bias",
                                                                        "Certainty of evidence (ignores imprecision)")
View(Grade_assessment_dataset_direct_wdae)
#Saving grade assessment details for WDAE mtxIR as an object of lists
grade_details_WDAE_mtxIR<-list(graph_list,WDAE_output_df2_list,Grade_assessment_dataset_direct_wdae,I_square_list)
saveRDS(grade_details_WDAE_mtxIR,file="List of WDAE MTX-IR data components for GRADE_AUG2026.rds")

wdae_output_one_study_list <- list()
wdae_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

# Count studies for each comparison
n_studies <- c()
for(i in 1:length(graph_list)){
  df<-WDAE_output_df2_list[[i]][-c(1,2,countdown(length(WDAE_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
}



# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
wdae_output_one_study_list  <- WDAE_output_df2_list[one_idx]
wdae_output_more_than_2_list <- WDAE_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_wdae[
    Grade_assessment_dataset_direct_wdae$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_wdae[
    Grade_assessment_dataset_direct_wdae$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_wdae_mtxIR_one_study <- list(
  graph_list_one_study_list,
  wdae_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_wdae_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  wdae_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_wdae_mtxIR_one_study,
  file = "List of WDAE MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_wdae_mtxIR_more_than_2,
  file = "List of WDAE MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)


###############################################################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_withdrawal_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct_wdae
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_wdae$Treatment, 
                                                Grade_assessment_dataset_direct_wdae$Intervention), 
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
#Getting the actual direct estimate for missing three arm trial
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(WDAE_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","RR"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 
glen_direct_evidence_rating_wdae<-read.csv("GLen_2026-09-03 22_16_WDAE-MTXIR_complete.csv")
glen_direct_evidence_rating_wdae2<-glen_direct_evidence_rating_wdae[,-c(13:15)]
colnames(glen_direct_evidence_rating_wdae2)[1:2]<-colnames(grade_table4)[1:2]

glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_wdae2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"WDAE"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"RR"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR_withdrawal_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR_withdrawal_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}


##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb
##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for WDAE_MTXIR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for WDAE_MTXIR_AUG2026 modified.csv",row.names = F)

##GETTING FINAL GRADE TABLE FOR WDAE OUTCOME
initial_gradeTable_wdae<-Temp_GRADE5

estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_wdae2<-initial_gradeTable_wdae[match(all_nma_estimates,initial_gradeTable_wdae$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()

#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_wdae2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_wdae2.1<-cbind(initial_gradeTable_wdae2[,c(1:14)],first_order_loop_contri,
                                  initial_gradeTable_wdae2[,-c(1:14)])

initial_gradeTable_wdae2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_wdae2.1,
                                                             all_nma_estimates,
                                                             estimate_contribution_matrix,
                                                             glen_direct_evidence_rating6_onlyDirect,
                                                             top_n = 2))

# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_wdae2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_wdae2.1<-result_df
initial_gradeTable_wdae2.1$combI<-str_replace(initial_gradeTable_wdae2.1$comb, ":", " vs ")

initial_gradeTable_wdae2.1<-process_dataset(initial_gradeTable_wdae2.1, 2)
initial_gradeTable_wdae3<-initial_gradeTable_wdae2.1[match(initial_gradeTable_wdae$comb,initial_gradeTable_wdae2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for WDAE_MTX-IR outcome_2026.csv",row.names = T)
#
###########Assessing Intransitivity
sub_data_intransitivity <- initial_gradeTable_wdae3[initial_gradeTable_wdae3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_wdae3[initial_gradeTable_wdae3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_wdae3.1<-cbind(initial_gradeTable_wdae3[,c(1:23)],
                                  Intransitivity="",
                                  `Certainty of evidence (ignores imprecision)-Indirect`="",
                                  initial_gradeTable_wdae3[,-c(1:23)])
initial_gradeTable_wdae3.1[initial_gradeTable_wdae3.1$Indirect!="-",c("Intransitivity",
                                                                      "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                            "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_wdae4<-initial_gradeTable_wdae3.1


#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)
grade_temp<-cbind(initial_gradeTable_wdae4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_wdae4[,-c(1:25)])
###Applying the P-values
P_values_nodesplit_wdaemtxIR2
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_wdaemtxIR3 <- P_values_nodesplit_wdaemtxIR2[match(drt$comb, P_values_nodesplit_wdaemtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_wdaemtxIR3$p,2)

subWDAESOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatWDAE<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subWDAESOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(nrow(subcontriMatWDAE))
IndirPer <- numeric(nrow(subcontriMatWDAE))

# Loop to calculate dirPer and IndirPer
for (i in 1:nrow(subcontriMatWDAE)) {
  dirPer[i] <- subcontriMatWDAE[i, colnames(subcontriMatWDAE) == row.names(subcontriMatWDAE)[i]]
  IndirPer[i] <- 1- dirPer[i]
}

# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")

grade_temp2[match(row.names(subcontriMatWDAE), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(row.names(subcontriMatWDAE), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)



##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
#When both direct and indirect are present
fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg

####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho


initial_gradeTable_wdae5<-grade_temp2.2
initial_gradeTable_wdae5<-initial_gradeTable_wdae5[,-c(30,33)]
colnames(initial_gradeTable_wdae5)[c(29,30)]<-c("Initial NMA evidence rating",
                                                "P-Value")

###After incoherence evaluation
initial_gradeTable_wdae5.1<-cbind(initial_gradeTable_wdae5[,c(1:31)],
                                  qoeAfterIncohenrence="",
                                  initial_gradeTable_wdae5[,-c(1:31)])
initial_gradeTable_wdae5.1[(initial_gradeTable_wdae5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_wdae5.1[(initial_gradeTable_wdae5.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_wdae5.1[initial_gradeTable_wdae5.1$Incoherence!="",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_wdae5.1[(initial_gradeTable_wdae5.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence

#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_wdae5.1[(initial_gradeTable_wdae5.1$Direct!="-")&(initial_gradeTable_wdae5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_wdae5.1[(initial_gradeTable_wdae5.1$Direct!="-")&(initial_gradeTable_wdae5.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_wdae6<-cbind(initial_gradeTable_wdae5.1,
                                p.estimate=apply(as.data.frame(initial_gradeTable_wdae5.1$NMA),1,function(x){
                                  as.numeric(sub(" .*", "", x))
                                }),
                                t(apply(as.data.frame(initial_gradeTable_wdae5.1$NMA),1,function(x){
                                  as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                  
                                })))

colnames(initial_gradeTable_wdae6)[c(36,37)]<-c("l.bound","u.bound")




row.names(initial_gradeTable_wdae6)<-NULL
##DATA
#WDAE_data<-current_WDAE_dataset_mtxIR_withROB2_indrectness2
WDAE_data<-WDAE_mtxIR_data_new_withROB2_indrectness2
WDAE_data<-WDAE_data[order(WDAE_data$study),c(1:6)]


WDAE_data_wide2<-reshape(number_of_arms(WDAE_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(WDAE_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(WDAE_data_wide2[,c(3,7,11)],1,xcv)

WDAE_data_wide3<-cbind(WDAE_data_wide2,na=numberOfArms)
View(WDAE_data_wide3)
colnames(WDAE_data_wide3)<-c("study", "study_name", "t.1","e.1","n.1", "r.1",
                             "t.2","e.2","n.2","r.2","t.3","e.3","n.3","r.3", "na")
                             #"t.4","e.4","n.4","r.4","na" )
#######
####Running Rates NMA model
# load data for MTC
MTCData_rates <- as.data.frame(WDAE_data_wide3)
r <- data.matrix(MTCData_rates[,c(6,10,14)])
n <- data.matrix(MTCData_rates[,c(5,9,13)])
t <- data.matrix(MTCData_rates[,c(3,7,11)])
E <- data.matrix(MTCData_rates[,c(4,8,12)])
na <- data.matrix(MTCData_rates[,"na"])
nt <- max(t, na.rm=TRUE)
ns <- nrow(r)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- r[t[,'t.1']==1,"r.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
E.PBO <- E[t[,'t.1']==1,"e.1"]

# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,nt-1)), mu=rep(0,ns),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
model.file4<-"Rate_nma_analysis.txt"
n.chains=3
model.inits <- initv1 

rates.model3<-jags.model(file=model.file4,
                         data= list("r"=r,"t"=t,"E"=E,"na" = na[,1], "nt" = nt, "ns" = ns,
                                    "NTpbo"= NTpbo,"xPBO"=xPBO,"E.PBO"=E.PBO),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=5000
)

update(object=rates.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("RR","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( rates.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]
typeof(estimates_relative_to_baseline)

#Converting estimates relative to baseline to absolutes
#nma_estimate_bio_dataframe3<-as.data.frame(apply(t(do.call(rbind,nma_estimate_bio_dataframe2)),2,as.numeric))
estimates_relative_to_baseline<-initial_gradeTable_wdae6[1:(nt-1),c(1,2,35:37)]
#estimates_relative_to_baseline<-initial_gradeTable_wdae6[initial_gradeTable_wdae6$treat2=="MTX",c(1,2,35:37)]

#absolutes_from_baseline<-cbind(baseline_risk,(estimates_relative_to_baseline[,3:5]*baseline_risk)/((1-baseline_risk)+(estimates_relative_to_baseline[,3:5]*baseline_risk)))
absolutes_from_baseline<-cbind(baseline_risk,1-exp(-estimates_relative_to_baseline[,3:5]*baseline_risk))
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(baseline_risk,estimates_relative_to_baseline_combine[,7]))
#1-exp(-RR*ACR)
colnames(all_baseline)<-"base"
com.list<-list()
unique_treat2<-unique(initial_gradeTable_wdae6$treat2)
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_wdae6[initial_gradeTable_wdae6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,1-exp(-d.sub[,3:5]*base))
  #absolutes_conversion<-cbind(base,(d.sub[,3:5]*base)/((1-base)+(d.sub[,3:5]*base)))
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)
com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9]*1000)
more_fewer_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"fewer"
  }else{
    text<-"more"
  }
  return(paste(round(abs(e1-e2),0),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,13]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,13])))

prefix.assumedrisk <- ''
prefix.corrrisk<-''
suffix.assumedrisk <- ' per 1000 pt-years'
suffix.corrrisk <- ' per 1000 pt-years'
result_com.dataframe_withtext<-cbind(com.dataframe_withtext[,1:2],
                                     paste(prefix.assumedrisk,round(com.dataframe_withtext[,6]*1000,0), suffix.assumedrisk, sep=''),
                                     paste(prefix.corrrisk, round(com.dataframe_withtext[,7]*1000,0), suffix.corrrisk, ' (', round(com.dataframe_withtext[,8]*1000,0), ' to ',round(com.dataframe_withtext[,9]*1000,0),')', sep=''),
                                     paste(prefix.corrrisk, com.dataframe_withtext[,13], ' (', com.dataframe_withtext[,14], ' to ',com.dataframe_withtext[,15],')', sep=''),
                                     # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                     com.dataframe_withtext[,10:15]
)

colnames(result_com.dataframe_withtext)[c(3:11)]<-c("Assumed risk (comparator)",	"Corresponding risk (intervention)",
                                                    "Rate difference","p.estimate","l.bound","u.bound","","","")



initial_gradeTable_wdae7<-cbind(initial_gradeTable_wdae6[,-c(35:37)],
                                result_com.dataframe_withtext[,-c(1,2)],pair="")
#colnames(initial_gradeTable_wdae7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_wdae8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_wdae7$treat1,mtxIR_withdrawal_subl)),
                                `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_wdae7$treat2,mtxIR_withdrawal_subl)),
                                initial_gradeTable_wdae7)
colnames(initial_gradeTable_wdae8)[3:4]<-c("Intervention",	"Comparator")

#write.csv(initial_gradeTable_ACR508,"ACR50 SoF table without imprecision assessment May2023.csv",row.names = F)

initial_gradeTable_wdae8<-initial_gradeTable_wdae8[,-c(46)]
write.csv(initial_gradeTable_wdae8,"WDAE MTXIR SoF table without imprecision assessment AUG2026.csv",row.names = F)
effect.size.rating_list<-lapply(split(initial_gradeTable_wdae7, seq_len(nrow(initial_gradeTable_wdae7))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 50, large.effect = 200, framework = T)
})

as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_wdae9<-initial_gradeTable_wdae8
initial_gradeTable_wdae9<-cbind(initial_gradeTable_wdae9[,c(1:34,37:45,35:36)])
initial_gradeTable_wdae9.1<-cbind(initial_gradeTable_wdae9[,c(1:37)],
                                  Small.effect.threshold= 50,
                                  Large.effect.threshold= 200,
                                  Framework = "NULL",
                                  Imprecision.1 = "",
                                  Certainty.of.evidence = "",
                                  initial_gradeTable_wdae9[,c(44,45)])


initial_gradeTable_wdae9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)
###Final certainty rating of the NMA
detailed_SOF_WDAE_mtxIR<-initial_gradeTable_wdae9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_WDAE_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_WDAE_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_WDAE_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_WDAE_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_WDAE_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_WDAE_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_WDAE_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_WDAE_mtxIR[is.na(detailed_SOF_WDAE_mtxIR$Certainty.of.evidence),])
detailed_SOF_WDAE_mtxIR[is.na(detailed_SOF_WDAE_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

#detailed_SOF_WDAE_mtxIR$comb<-initial_gradeTable_wdae7$comb

write.csv(detailed_SOF_WDAE_mtxIR,"WDAE_MTX-IR_Detail_SoF_OCT2026.csv", row.names = F)

###############################################################################
###############################################################################
#####SAE

#SAE_dataset_naive<-final_dmardDataset[,c(1:51,942:961)]
SAE_dataset_mtxIR<-final_mtx_irDataset_actual%>% 
  select(matches("^Serious.Adverse.Events..SAE\\..*")) %>%
  select(-contains("Reported"))%>%
  select(-contains(c("..4.months","22.Ws..")))
colnames(SAE_dataset_mtxIR)
View(SAE_dataset_mtxIR)
####
SAE_dataset_mtxIR[SAE_dataset_mtxIR=="NA"]<-NA
# ###############################################################################
library(tibble)
clean_withdrawal <- function(dataset, month) {
  
  if(month == 12){
    
    temp <- dataset %>%
      select(matches("^Serious.Adverse.Events..SAE\\..*")) %>%
      select(contains("12.months"))
    
  }
  
  if(month == 6){
    
    temp <- dataset %>%
      select(matches("^Serious.Adverse.Events..SAE\\..*")) %>%
      select(contains("6.months"))
    
  }
  
  if(month == 3){
    
    temp <- dataset %>%
      select(matches("^Serious.Adverse.Events..SAE\\..*")) %>%
      select(contains("3.months"))
    
  }
  
  if(month == "EndOfTrial"){
    
    temp <- dataset %>%
      select(matches("^Serious.Adverse.Events..SAE\\..*")) %>%
      select(contains("End.of.trial"))
  }
  
  return(temp)
}

View(clean_withdrawal(SAE_dataset_mtxIR,3))

clean_withdrawal_all<-function(dataset){
  #3months
  clean_withdrawal_3months<-clean_withdrawal(dataset,3)
  clean_withdrawal_3months<-add_column(clean_withdrawal_3months, "edit" = clean_withdrawal_3months[,colnames(clean_withdrawal_3months)[4]], .after = colnames(clean_withdrawal_3months)[4])
  colnames(clean_withdrawal_3months)[5]<-paste0(colnames(clean_withdrawal_3months)[4],"_edit")
  index3months<-which(is.na(clean_withdrawal_3months[,colnames(clean_withdrawal_3months)[4]]))
  sub_clean_withdrawal_3months<-clean_withdrawal_3months[index3months,]
  clean_withdrawal_3months[index3months,colnames(clean_withdrawal_3months)[5]]<-sub_clean_withdrawal_3months[,3]
  #6months
  clean_withdrawal_6months<-clean_withdrawal(dataset,6)
  clean_withdrawal_6months<-add_column(clean_withdrawal_6months, "edit" = clean_withdrawal_6months[,colnames(clean_withdrawal_6months)[4]], .after = colnames(clean_withdrawal_6months)[4])
  colnames(clean_withdrawal_6months)[5]<-paste0(colnames(clean_withdrawal_6months)[4],"_edit")
  index6months<-which(is.na(clean_withdrawal_6months[,colnames(clean_withdrawal_6months)[4]]))
  sub_clean_withdrawal_6months<-clean_withdrawal_6months[index6months,]
  clean_withdrawal_6months[index6months,colnames(clean_withdrawal_6months)[5]]<-sub_clean_withdrawal_6months[,3]
  #12months
  clean_withdrawal_12months<-clean_withdrawal(dataset,12)
  clean_withdrawal_12months<-add_column(clean_withdrawal_12months, "edit" = clean_withdrawal_12months[,colnames(clean_withdrawal_12months)[4]], .after = colnames(clean_withdrawal_12months)[4])
  colnames(clean_withdrawal_12months)[5]<-paste0(colnames(clean_withdrawal_12months)[4],"_edit")
  index12months<-which(is.na(clean_withdrawal_12months[,colnames(clean_withdrawal_12months)[4]]))
  sub_clean_withdrawal_12months<-clean_withdrawal_12months[index12months,]
  clean_withdrawal_12months[index12months,colnames(clean_withdrawal_12months)[5]]<-sub_clean_withdrawal_12months[,3]
  #EndOfTrial
  clean_withdrawal_EndOfTrial<-clean_withdrawal(dataset,"EndOfTrial")
  clean_withdrawal_EndOfTrial<-add_column(clean_withdrawal_EndOfTrial, "edit" = clean_withdrawal_EndOfTrial[,colnames(clean_withdrawal_EndOfTrial)[4]], .after = colnames(clean_withdrawal_EndOfTrial)[4])
  colnames(clean_withdrawal_EndOfTrial)[5]<-paste0(colnames(clean_withdrawal_EndOfTrial)[4],"_edit")
  indexEndOfTrial<-which(is.na(clean_withdrawal_EndOfTrial[,colnames(clean_withdrawal_EndOfTrial)[4]]))
  sub_clean_withdrawal_EndOfTrial<-clean_withdrawal_EndOfTrial[indexEndOfTrial,]
  clean_withdrawal_EndOfTrial[indexEndOfTrial,colnames(clean_withdrawal_EndOfTrial)[5]]<-sub_clean_withdrawal_EndOfTrial[,3]
  
  
  return(cbind(clean_withdrawal_3months,clean_withdrawal_6months,clean_withdrawal_12months,clean_withdrawal_EndOfTrial))
}
View(clean_withdrawal_all(SAE_dataset_mtxIR))

###Get new SAE dataset
SAE_dataset_mtxIR1<-cbind(final_mtx_irDataset_actual[,c(1:51)],clean_withdrawal_all(SAE_dataset_mtxIR))
View(SAE_dataset_mtxIR1)

SAE_dataset_mtxIR1[SAE_dataset_mtxIR1=="NA"]<-NA
write.csv(SAE_dataset_mtxIR1,"SAE dataset_mtxIR_AUG2026.csv",row.names = F)

SAE_dataset_mtxIR2<- SAE_dataset_mtxIR1

# Columns by time point
sae_3m_cols <- grep(
  "^Serious\\.Adverse\\.Events\\.\\.SAE\\.\\.3\\.months",
  names(SAE_dataset_mtxIR2),
  value = TRUE
)

sae_6m_cols <- grep(
  "^Serious\\.Adverse\\.Events\\.\\.SAE\\.\\.6\\.months",
  names(SAE_dataset_mtxIR2),
  value = TRUE
)

sae_12m_cols <- grep(
  "^Serious\\.Adverse\\.Events\\.\\.SAE\\.\\.12\\.months",
  names(SAE_dataset_mtxIR2),
  value = TRUE
)

sae_eot_cols <- grep(
  "^Serious\\.Adverse\\.Events\\.\\.SAE\\.\\.End\\.of\\.trial",
  names(SAE_dataset_mtxIR2),
  value = TRUE
)

# End of trial
SAE_dataset_mtxIR2[
  SAE_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,6,9,11,12,13,18),
  sae_eot_cols
] <- NA

# 12 months
SAE_dataset_mtxIR2[
  SAE_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,23,24),
  sae_12m_cols
] <- NA

# 6 months
SAE_dataset_mtxIR2[
  SAE_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(3,4,9,11,12,13,18,23,24),
  sae_6m_cols
] <- NA

# 3 months
SAE_dataset_mtxIR2[
  SAE_dataset_mtxIR2$`Trial.duration..in.months` %in%
    c(6,9,11,12,13,18,23,24),
  sae_3m_cols
] <- NA



empty_SAEstudy_mtxIR<-list()  
keep_SAEStudy_mtxIR_for_errorCheck<-list()
SAEStudy_for_analysis_mtxIR<-list()


i=37
colindex<-as.vector(which(colSums(!is.na(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,c(1:3,52:75)])) > 0))
temp_SAE<-SAE_dataset_mtxIR2[,c(1:3,52:75)]
temp_SAE2<-temp_SAE[temp_SAE$Study_no==i,colindex]
temp_SAE2

# ##Fleischmann 2017
# SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Fleischmann 2017","Serious.Adverse.Events..SAE..12.months.Exposure..patient.years."]<-c(((161*52)+((10)*(52/2)))/52,((136*52)+((9)*(52/2)))/52,((173*52)+((23)*(52/2)))/52)
# 
# ##Breedveld 2006 
# SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Breedveld 2006",c("Serious.Adverse.Events..SAE..End.of.trial.Rate.of.events..per.100.patient.years.","Serious.Adverse.Events..SAE..End.of.trial.Total.number.of.events_edit")]<- round(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Breedveld 2006",c("Serious.Adverse.Events..SAE..End.of.trial.Rate.of.events..per.100.patient.years.","Serious.Adverse.Events..SAE..End.of.trial.Total.number.of.events_edit")],0)
#Kremer 2006
SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Kremer 2006",sae_12m_cols]<-NA

#Kay 2008 
SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Kay 2008",sae_12m_cols]<-NA

#Emery 2010 
SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Emery 2010",sae_12m_cols]<-NA

#Schiff 2008
SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="Schiff 2008",sae_12m_cols]<-NA

#vanderHeijde 2013
SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study.Identifier=="vanderHeijde 2013",sae_12m_cols]<-NA

for(i in 1:length(unique(SAE_dataset_mtxIR2$Study_no))){
  if(ncol(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,c(52:75)][colSums(!is.na(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,c(52:75)])) > 0])==0){
    empty_SAEstudy_mtxIR[[i]]<-SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,]
  }else if(!(any(c(1,7,13,19) %in% as.vector(which(colSums(!is.na(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,c(52:75)])) > 0))))){
    keep_SAEStudy_mtxIR_for_errorCheck[[i]]<-SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,]
  }else{
    colindex<-as.vector(which(colSums(!is.na(SAE_dataset_mtxIR2[SAE_dataset_mtxIR2$Study_no==i,c(1:3,52:75)])) > 0))
    temp_SAE<-SAE_dataset_mtxIR2[,c(1:3,52:75)]
    temp_SAE2<-temp_SAE[temp_SAE$Study_no==i,colindex]
    df_sae_lastTP <- temp_SAE2 %>%
      select(matches("SAE.*(Exposure|N|events_edit)")) %>%
      {
        cols <- names(.)
        
        suffix <- str_extract(cols, "(Exposure|N|events_edit)")
        timept <- str_extract(cols, "(\\d+\\.months|End\\.of\\.trial)")
        
        time_order   <- c("3.months","6.months","12.months","End.of.trial")
        suffix_order <- c("Exposure","N","events_edit")
        
        # order columns correctly first
        ord <- order(match(timept, time_order),
                     match(suffix, suffix_order))
        
        df <- .[, ord]
        
        # ---- NEW PART: keep only LAST available timepoint ----
        ordered_time_present <- time_order[time_order %in% timept]
        last_time <- tail(ordered_time_present, 1)
        
        df[, str_detect(names(df), last_time)]
      } #%>%dplyr::select(-dplyr::matches("Total\\.number\\.of\\.events$"))
    
    temp_SAE3 <- cbind(temp_SAE2[,c(1:3)],df_sae_lastTP[c(1:3)], colnames(temp_SAE2)[4])
    colnames(temp_SAE3)<-c("Study_no","Study.Identifier","Intervention","exposure","sampleSize","responders","month")
    SAEStudy_for_analysis_mtxIR[[i]]<-temp_SAE3
  }
  
  
}

View(do.call(rbind,SAEStudy_for_analysis_mtxIR))
View(do.call(rbind,empty_SAEstudy_mtxIR))


current_SAE_dataset<-do.call(rbind,SAEStudy_for_analysis_mtxIR)

#Corrections
current_SAE_dataset[current_SAE_dataset$Study.Identifier=="Emery 2010",c("exposure","responders")][2,]<-c(77.54,6)


##Original needed to compare after checking for errors
current_SAE_dataset$Study_no<-renumber(current_SAE_dataset,"Study_no")
row.names(current_SAE_dataset)<-NULL




# #exclude studies Emery 2013(cross over design)
# current_SAE_dataset<-current_SAE_dataset[!current_SAE_dataset$Study_no%in%c(8),]
# current_SAE_dataset$Study_no<-renumber(current_SAE_dataset,"Study_no")
# row.names(current_SAE_dataset)<-NULL

current_SAE_dataset$exposure<-round(as.numeric(current_SAE_dataset$exposure),2)
current_SAE_dataset$sampleSize<-as.numeric(current_SAE_dataset$sampleSize)
current_SAE_dataset$responders<-round(as.numeric(current_SAE_dataset$responders),0)


current_SAE_dataset3<- current_SAE_dataset[order(current_SAE_dataset$Study_no),]
#####Formatting the month column
current_SAE_dataset3$month<- stringr::str_extract(current_SAE_dataset3$month, "\\d+(?=\\.months)") %>%
  as.numeric()
current_SAE_dataset3$month<- ifelse(is.na(current_SAE_dataset3$month),"EndOfTrial",current_SAE_dataset3$month)

##Mean and Median of time frame
mean(current_SAE_dataset3$month)
median(current_SAE_dataset3$month)

View(current_SAE_dataset3)
current_SAE_dataset_mtxIR<-current_SAE_dataset3
#colnames(current_SAE_dataset_mtxIR)[3]<-"treatment"
current_SAE_dataset_mtxIR <- current_SAE_dataset_mtxIR%>%
  dplyr::rename("treatment" = Intervention)

##Drop the high risk of bias studies
current_SAE_dataset_mtxIR<- current_SAE_dataset_mtxIR %>%
  filter(!Study.Identifier %in%c("Keystone 2004", "Smolen 2009", "Maini 2006", "Gashi 2014", "Yamanaka 2016", "Kim 2012", "Maini 1998", "Dougados 2005", "Fleischmann 2009",
                                 "Porter 2016", "Kim 2000", "Karanikolas 2006", "DeFilippis 2006", "Wong 2009", "Huang 2013", "Takeuchi 2013.b", "Heimans 2014", "Jobanputra 2012",
                                 "Dougados 2017", "Jeurissen 1991", "Kang 2018", "Matsubara 2018", "Behrens 2021", "Peterfy 2013", "Solomon 2023", "KowoonJoo 2013", "Hartman 2023.a",
                                 "Kameda 2011","Li 2016","Tam 2012"))



library(dplyr)
for(i in 1:nrow(current_SAE_dataset_mtxIR)){
  ppp<- master.medtable %>% filter(current_SAE_dataset_mtxIR$treatment[i] ==master.medtable$medname )
  current_SAE_dataset_mtxIR$treatment[i]<-ppp[1,1]
}

View(do.call(rbind,empty_SAEstudy_mtxIR))
View(do.call(rbind,keep_SAEStudy_mtxIR_for_errorCheck))

##Creating initial legend
mtxIR.SAE.medlegend <- new.medtable(current_SAE_dataset_mtxIR,'treatment')
mtxIR.SAE.medlegend


#Getting the quality assessment for each study
###Adding ROB data
ffff<-unique(current_SAE_dataset_mtxIR$Study.Identifier)
data_qq<-current_SAE_dataset_mtxIR
View(ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff,])
current_assessment_for_studies_mtxIR<-ROB_quality_assess_data_consensus_mtxIR[ROB_quality_assess_data_consensus_mtxIR$Study.Id%in%ffff,]
col_names<-c("Study.Id", "SEQUENCE.GENERATION", "ALLOCATION.CONCEALMENT","BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..RADIOGRAPHIC.OUTCOMES",    
             "BLINDING.OF.PARTICIPANTS.AND.PERSONNEL..NON.RADIOGRAPHIC.OUTCOMES", "BLINDING.OF.OUTCOME.ASSESSORS..RADIOGRAPHIC.OUTCOMES",             
             "BLINDING.OF.OUTCOME.ASSESSORS..NON.RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..ALL.OTHER.OUTCOMES",                      
             "INCOMPLETE.OUTCOME.DATA..RADIOGRAPHIC.OUTCOMES", "INCOMPLETE.OUTCOME.DATA..WITHDRAWALS",                             
             "SELECTIVE.OUTCOME.REPORTING", "OTHER.BIAS", "OVERALL..RADIOGRAPHIC.OUTCOMES", "OVERALL..WITHDRAWALS",                                             
             "OVERALL..ALL.OTHER.OUTCOMES"  )
current_assessment_for_studies_mtxIR<- current_assessment_for_studies_mtxIR %>%
  select(c(col_names))
list_appraisal<-list()
for(i in 1:length(ffff)){
  temp11<-data_qq[data_qq$Study.Identifier==ffff[i],]
  temp12<-current_assessment_for_studies_mtxIR[current_assessment_for_studies_mtxIR$Study.Id==ffff[i],-1]
  temp13<-cbind(temp11,temp12)
  list_appraisal[[i]]<-temp13
}
current_SAE_dataset_mtxIR_withROB<-do.call(rbind,list_appraisal)
View(current_SAE_dataset_mtxIR_withROB)
current_SAE_dataset_mtxIR_withROB2<-current_SAE_dataset_mtxIR_withROB%>%
  select(Study_no,Study.Identifier,treatment,exposure,sampleSize,responders,month,OVERALL..ALL.OTHER.OUTCOMES)
View(current_SAE_dataset_mtxIR_withROB2)

###Adding Indirectness data
nnnn<-unique(current_SAE_dataset_mtxIR_withROB2$Study.Identifier)
data_nn<-current_SAE_dataset_mtxIR_withROB2
resolved_indirectness_consensus_temp<-read_excel("Indirectness Assessment DMARD-IR final.xlsx")
resolved_indirectness_consensus<-resolved_indirectness_consensus_temp%>%
  dplyr::rename( "Consensus"= `Final Rating`)

View(resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%nnnn,])
current_indirect_assessment_for_studies<-resolved_indirectness_consensus[resolved_indirectness_consensus$`Study Identifier`%in%nnnn,]
#current_indirect_assessment_for_studies<-current_indirect_assessment_for_studies[-19,] #remove dejong one arm
list_appraisal2<-list()

for(i in 1:length(nnnn)){
  temp11<-data_nn[data_nn$Study.Identifier==nnnn[i],]
  temp12<-current_indirect_assessment_for_studies[current_indirect_assessment_for_studies$`Study Identifier` ==nnnn[i],-1][,22]
  # temp13<-cbind(temp11,temp12)
  if(nrow(temp12)==0){
    temp13<-cbind(temp11,Consensus="")
  }else{
    temp13<-cbind(temp11,temp12)
  }
  list_appraisal2[[i]]<-temp13
}
current_SAE_dataset_mtxIR_withROB2_indrectness<-do.call(rbind,list_appraisal2)
#colnames(current_SAE_dataset_mtxIR_withROB2_indrectness)[9]<-"Indirectness"
current_SAE_dataset_mtxIR_withROB2_indrectness <- current_SAE_dataset_mtxIR_withROB2_indrectness%>%
  dplyr:: rename("Indirectness" = Consensus)
#ACR50_try_data_new_withROB2_indrectness<-cbind(ACR50_try_data_new_withROB2,Indirectness="low")
current_SAE_dataset_mtxIR_withROB2_indrectness$Indirectness<-tolower(current_SAE_dataset_mtxIR_withROB2_indrectness$Indirectness)
View(current_SAE_dataset_mtxIR_withROB2_indrectness)



##NMA plot reason for narrowing of network
nma.plot(current_SAE_dataset_mtxIR_withROB2_indrectness[,1],current_SAE_dataset_mtxIR_withROB2_indrectness[,3],data = current_SAE_dataset_mtxIR_withROB2_indrectness, trtname = as.vector(mtxIR.SAE.medlegend$medname))


################################################################################################################
#trying the narrowing a network function
#############################################################################################
current_SAE_dataset_mtxIR_withROB2_indrectness
###Removed these studies with 0 responder arms because the made the model to fail to converge
current_SAE_dataset_mtxIR_withROB2_indrectness<-current_SAE_dataset_mtxIR_withROB2_indrectness[!current_SAE_dataset_mtxIR_withROB2_indrectness$Study.Identifier%in%c("MacIsaac 2014"),]
#colnames(current_SAE_dataset_mtxIR_withROB2_indrectness)[c(1,2)]<-c("study","study_name")
current_SAE_dataset_mtxIR_withROB2_indrectness <- current_SAE_dataset_mtxIR_withROB2_indrectness%>%
  dplyr::rename("study" = Study_no,"study_name"= Study.Identifier)

#Narrow 1
mtxIR_SAE_narrow<-narrow.all(current_SAE_dataset_mtxIR_withROB2_indrectness)
mtxIR_SAE_narrow


#sub network plot 1
mtxIR.SAE.d<-mtxIR_SAE_narrow[[1]]
length(unique(mtxIR.SAE.d$treatment))
mtxIR_SAE_subl<-new.medtable(mtxIR.SAE.d,'treatment')

##To download dataset
mtxIR.SAE.d_todownload <- mtxIR.SAE.d
write.csv(mtxIR.SAE.d_todownload, "SAE mtxIR dataset AUG2026.csv")
# 1. Open jpeg file
jpeg("Serious Adverse Events (SAE) Network plot with named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.SAE.d[,1],mtxIR.SAE.d[,3],data = mtxIR.SAE.d, trtname = as.vector(mtxIR_SAE_subl$medname), title = "Serious adverse events (SAE) Network plot-DMARD-mtxIR")

# Close the pdf file
dev.off() 
# 2. Open jpeg file
jpeg("Serious Adverse Events (SAE) Network plot without named interventions_2026.jpg", width = 3050, height = 2050, res=300)
nma.plot(mtxIR.SAE.d[,1],mtxIR.SAE.d[,3],data = mtxIR.SAE.d, title = "Serious adverse events (SAE) Network plot-DMARD-mtxIR")

# Close the pdf file
dev.off() 

##############################################################################
#Performing sub network plot 3
library(coda)
library(gemtc)
View(mtxIR.SAE.d)
mtxIR_SAE_subl
#Correcting for zero cells
mtxIR.SAE.d3<-cbind(mtxIR.SAE.d,
                    responders_corr_zero_cells=ifelse(mtxIR.SAE.d$responders==0,0.5,
                                                      mtxIR.SAE.d$responders))

#Sorting to start with 1 as reference treatment
##Function
reorder_reference<-function(data,Study_no){
  ss<- data[data$study==Study_no,]
  return(ss[order(ss$treatment), ])
}


mtxIR.SAE.d3$study<-renumber(mtxIR.SAE.d3,"study")
for(i in 1:length(unique(mtxIR.SAE.d3$study))){
  mtxIR.SAE.d3[mtxIR.SAE.d3$study==i,]<-reorder_reference(mtxIR.SAE.d3,i)
}
row.names(mtxIR.SAE.d3)<-NULL

data_SAE_mtxIR_1<-mtxIR.SAE.d3
data_SAE_mtxIR_1$treatment<-renumber(data_SAE_mtxIR_1,"treatment")
data_SAE_mtxIR_1$study<-renumber(data_SAE_mtxIR_1,"study")
#colnames(data_SAE_mtxIR_1)[c(6,10)]<-c("responders_corr_zero_cells","responders")
nma.plot(data_SAE_mtxIR_1$study,data_SAE_mtxIR_1$treatment,data = data_SAE_mtxIR_1, trtname = as.vector(mtxIR_SAE_subl$medname), title = "Serious adverse events (SAE) Network plot-DMARD-mtxIR")
network_data_SAE_mtxIR_1<-mtc.network(data_SAE_mtxIR_1)

##RUN the model
set.seed(522)
network_data_SAE_mtxIR_1.model <- mtc.model(network_data_SAE_mtxIR_1, likelihood='poisson', link='log',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))

#By default, the model will have 4 chains - generate a seed for each
seeds <- sample.int(4, n = .Machine$integer.max)
# Apply JAGS RNG settings to each chain
network_data_SAE_mtxIR_1.model$inits <- mapply(c, network_data_SAE_mtxIR_1.model$inits, list(
  list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
  list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
  list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
  list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
results_mtxIR_SAE<-mtc.run(network_data_SAE_mtxIR_1.model,n.iter = 50000, n.adapt = 30000)
summary(results_mtxIR_SAE)


View(as.data.frame(relative.effect.table.OR(results_mtxIR_SAE)))
refHR_mtxIR_SAE<-as.data.frame(relative.effect.table.OR2(results_mtxIR_SAE, comparator = "1"))

#library for mapvalues
library(plyr)
row.names(refHR_mtxIR_SAE)<-mapvalues(row.names(refHR_mtxIR_SAE),
                                      from = c(1:(nrow(mtxIR_SAE_subl))),
                                      to = as.vector(mtxIR_SAE_subl$medname))
colnames(refHR_mtxIR_SAE)<-mapvalues(colnames(refHR_mtxIR_SAE),
                                     from = c(1:(nrow(mtxIR_SAE_subl))),
                                     to = as.vector(mtxIR_SAE_subl$medname))
View(refHR_mtxIR_SAE)

index_diagonal3<-data.frame(index1=1:nrow(refHR_mtxIR_SAE),index2=1:nrow(refHR_mtxIR_SAE))
options(warn=-1)
for(i in 1:nrow(refHR_mtxIR_SAE)){
  refHR_mtxIR_SAE[index_diagonal3$index1[i],index_diagonal3$index2[i]]<-selected_treatment_name(as.numeric(refHR_mtxIR_SAE[row(refHR_mtxIR_SAE)==col(refHR_mtxIR_SAE)])[i],mtxIR_SAE_subl)
}
options(warn=0)
View(refHR_mtxIR_SAE)
#write.xlsx(refHR_mtxIR_SAE,"Results for Biologics.xlsx",sheetName = "SAE",row.names = T, append = T)
#write.xlsx(refHR_bio_withdrawal,"Results for Biologics Informative Prior.xlsx",sheetName = "WDAE",row.names = T, append = T)
write.csv(refHR_mtxIR_SAE,"SAE MTX-IR pairwise comparisons_OCT2026.csv", row.names = T)
#Add to existing workbook sheet of results
wb <- loadWorkbook("Results for MTX-IR_OCT2026.xlsx")

addWorksheet(wb, "SAE")
writeData(wb, "SAE", refHR_mtxIR_SAE, rowNames = TRUE)

saveWorkbook(wb, "Results for MTX-IR_OCT2026.xlsx", overwrite = TRUE)

##GRADE Assessment
nma_estimate_mtxIR_SAE_dataframe<-nma.estimate.table.OR(results_mtxIR_SAE)
View(nma_estimate_mtxIR_SAE_dataframe)

#Direct estimate and Indirect estimate
##Node splitting
network_data_SAE_mtxIR_1.ns<-mtc.nodesplit(network_data_SAE_mtxIR_1, thin=50,likelihood='poisson', link='log',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma")) 
names(network_data_SAE_mtxIR_1.ns)
mtc.nodesplit.comparisons(network_data_SAE_mtxIR_1.ns)


# Overall summary and plot
summary.network_data_SAE_mtxIR_1.ns <- summary(network_data_SAE_mtxIR_1.ns)


#Direct effects
direct.data<-summary.network_data_SAE_mtxIR_1.ns$dir.effect

#Indirect effects
indirect.data<-summary.network_data_SAE_mtxIR_1.ns$ind.effect

#Network estimate
network.data<-summary.network_data_SAE_mtxIR_1.ns$cons.effect

#P-values
P_values_nodesplit_saemtxIR<-summary.network_data_SAE_mtxIR_1.ns$p.value
P_values_nodesplit_saemtxIR[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_saemtxIR[,c(1,2)],1,switch_func)),1,as.numeric))
P_values_nodesplit_saemtxIR2<-P_values_nodesplit_saemtxIR
P_values_nodesplit_saemtxIR2[,c(1,2)]<-t(apply(t(apply(P_values_nodesplit_saemtxIR2[,c(1,2)],1,switch_func2)),1,as.numeric))
P_values_nodesplit_saemtxIR2$comb<-paste0(P_values_nodesplit_saemtxIR2$t1,":",P_values_nodesplit_saemtxIR2$t2)


direct_indirect_estimate<-cbind(direct.data[,c(2,1)],round(exp(direct.data[,c(3:5)]),2),round(exp(indirect.data[,c(3:5)]),2))
colnames(direct_indirect_estimate)<-c("treat1","treat2","RR1","low1","high1","RR2","low2","high2")

direct_indirect_estimate_list<-list()

for(i in 1:nrow(direct_indirect_estimate)){
  if(as.numeric(as.character(direct_indirect_estimate$treat1[i]))>as.numeric(as.character(direct_indirect_estimate$treat2[i]))){
    direct_indirect_estimate_list[[i]]<-direct_indirect_estimate[i,]
  }else{
    temp<-cbind(direct_indirect_estimate[i,2],direct_indirect_estimate[i,1],
                round(1/direct_indirect_estimate[i,3:8],2))
    temp2<-cbind(temp[,1:3],temp[,c(5,4)],temp[,c(6)],temp[,c(8,7)])
    colnames(temp2)<-c("treat1","treat2","RR1","low1","high1","RR2","low2","high2")
    direct_indirect_estimate_list[[i]]<-temp2
  }
}

direct_indirect_estimate<-do.call(rbind,direct_indirect_estimate_list)

grade_table1<-merge(direct_indirect_estimate,nma_estimate_mtxIR_SAE_dataframe,
                    by=c("treat1","treat2"),all.x=T,all.y=T)

current_SAE_dataset_mtxIR_withROB2_indrectness2<-mtxIR.SAE.d3
current_SAE_dataset_mtxIR_withROB2_indrectness2$treatment<-renumber(current_SAE_dataset_mtxIR_withROB2_indrectness2,"treatment")
current_SAE_dataset_mtxIR_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES<-ifelse(current_SAE_dataset_mtxIR_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES=="unclear","some concerns",current_SAE_dataset_mtxIR_withROB2_indrectness2$OVERALL..ALL.OTHER.OUTCOMES)


pairwise_SAE_RR<-pairwise(treat=treatment,event = responders_corr_zero_cells,data=network_data_SAE_mtxIR_1$data.ab,
                          studlab = study,n=sampleSize,sm="RR")



netPairwiseSAE<-netmeta(pairwise_SAE_RR)
netPairwiseSAE_split<-netsplit(netPairwiseSAE)
View(netPairwiseSAE_split$direct.common)
dat.net.contrib<-netcontrib(netPairwiseSAE, method = "r")
View(dat.net.contrib$random)

direct_evidence_new<-cbind(netPairwiseSAE_split$direct.common[,c(1:5)],netPairwiseSAE_split$indirect.common[,c(2:5)])
View(direct_evidence_new)
#getting only direct estimates
#direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,4:6]), 1,sum)==3,]
direct_only<-direct_evidence_new[apply(is.na(direct_evidence_new[,6:9]), 1,sum)==4,]

library(stringr)
direct_only1<-str_split_fixed(direct_only$comparison, ":", 2)
direct_only2<-cbind(as.numeric(direct_only1[,1]),as.numeric(direct_only1[,2]))
switch_func<-function(x){
  if(which.min(x)==1){
    return(sort(x,decreasing = T))
  }else{
    return(x)
  }
}
direct_only3<-t(apply(direct_only2,1,switch_func))
direct_only4<-paste(direct_only3[,1],direct_only3[,2],sep = ":")
#direct_only4%in%paste(nma_estimate_bio_HAQ_dataframe$treat1,nma_estimate_bio_HAQ_dataframe$treat2,sep = ":")

#Setting up initial grade table
##Assigning nma estimate to direct only
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%direct_only4,c(9:11)]

##Getting direct evidence for intervention in only multi-arm trials
zzz<-t(apply(apply(str_split_fixed(colnames(dat.net.contrib$random), ":", 2),2,as.numeric),1,switch_func))
colnames(zzz)<-c("treat1", "treat2")
grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)]
inthreeArmtrial<-anti_join(as.data.frame(zzz), as.data.frame(apply(grade_table1[apply(!is.na(grade_table1[,3:5]),1,sum)==3,c(1:2)],2,as.numeric)), by = c("treat1","treat2"))
grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(3:5)]<-grade_table1[paste(grade_table1$treat1,grade_table1$treat2,sep = ":")%in%paste(inthreeArmtrial$treat1,inthreeArmtrial$treat2,sep = ":"),c(9:11)]

##Assigning nma estimate to indirect only
grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(6:8)]<-grade_table1[apply(is.na(grade_table1[,5:7]),1,sum)==3,c(9:11)]





grade_table2<-cbind(grade_table1[,c(1,2)],Direct=paste(grade_table1[,3]," (",grade_table1[,4],", ",grade_table1[,5],")",sep = ""),
                    Indirect=paste(grade_table1[,6]," (",grade_table1[,7],", ",grade_table1[,8],")",sep = ""),
                    NMA=paste(grade_table1[,9]," (",grade_table1[,10],", ",grade_table1[,11],")",sep = "")
)
grade_table3<-grade_table2
grade_table3$treat1<-selected_treatment_name(as.numeric(as.character(grade_table3$treat1)),mtxIR_SAE_subl)
grade_table3$treat2<-selected_treatment_name(as.numeric(as.character(grade_table3$treat2)),mtxIR_SAE_subl)
replace_fun<-function(x){
  ifelse(x=="NA (NA, NA)","-",x)
}
grade_table4<-as.data.frame(t(apply(grade_table3, 1,replace_fun )))

#####################################
#####################################
###Getting only the direct evidence and forest plots for the grade assessment
all_direct_evidence<-grade_table1[apply(is.na(grade_table1[,3:5]),1,sum)!=3,]
all_direct_evidence2<-cbind(all_direct_evidence,comb=paste0(all_direct_evidence$treat1,":",all_direct_evidence$treat2))
library(Matrix)
library(metafor)
library(coda)
library(rjags)
library(tidyr)
graph_list<-list()
pooled_estimate2<-data.frame()
tau_square2<-data.frame()
SAE_output_df2_list<-list()
I_square_list<-list()
samps2<-list()




for(i in 1:nrow(all_direct_evidence2)){
  fgh<-list()
  tau_square<-list()
  Studies_contributing_to_direct_evidence<-current_SAE_dataset_mtxIR_withROB2_indrectness2[current_SAE_dataset_mtxIR_withROB2_indrectness2$treatment%in%as.numeric(c(t(all_direct_evidence2[i,c(1,2)]))),-10]
  
  Studies_contributing_to_direct_evidence_wide<-to.wide(as.data.frame(Studies_contributing_to_direct_evidence%>%
                                                                        group_by(study) %>%
                                                                        filter(n()>1)),
                                                        study="study_name", grp="treatment",grpvars = 4:6)
  
  temp11<-Studies_contributing_to_direct_evidence_wide
  study_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"study_name"]
  temp11.1<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),]%>%
    select("study", "exposure.1", "responders.1","exposure.2", "responders.2" )
  ROB_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"ROB"]
  Indirectness_ordered<-temp11[order(as.numeric(gsub(".*? ", "", temp11$study_name))),"Indirectness"]
  colnames(temp11.1)<-c("NS","E.t","x.t","E.c","x.c")
  #temp11.1<-temp11.1[,-1]
  row.names(temp11.1)<-NULL
  for(k in 1:nrow(temp11.1)){
    temp12<-temp11.1[k,]
    # 
    # if(temp12[,"x.t"]==0 & temp12[,"x.c"]==0){
    #   temp12[,c(3,5)]<-0.5
    # }
    
    model.file<-"Rate_direct_onestudy.txt"
    n.chains=3
    set.seed(2453)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
                                   "x.c"=temp12$x.c,
                                   "E.c"=temp12$E.c),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("RR")
    ## run JAGS and save posterior samples
    samps <- coda.samples( jj.model, params, n.iter=10000 )
    # model.file<-"Rate_direct_multiplestudy_FE.txt"
    # n.chains=3
    # set.seed(2553)
    # seeds <- sample.int(n.chains, n = .Machine$integer.max)
    # model.inits <- list(
    #   list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
    #   list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
    #   list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    # jj.model<-jags.model(file=model.file,
    #                      data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
    #                                "x.c"=temp12$x.c,
    #                                "E.c"=temp12$E.c,
    #                                "NS"=length(temp12$NS)),
    #                      n.chains=n.chains,
    #                      inits = model.inits,
    #                      n.adapt=500
    # )
    # update(object=jj.model,n.iter=10000,progress.bar='text')
    # ## specify parameters to be monitored
    # params <- c("RR")
    # ## run JAGS and save posterior samples
    # samps <- coda.samples( jj.model, params, n.iter=10000 )
    summary_entry1<-summary(samps)
    if(round(summary_entry1$statistics[1],1)==round(summary_entry1$quantiles[3],1)){
      fgh[[k]]<-data.frame(RR=summary_entry1$statistics[1],
                           Lower=summary_entry1$quantiles[c(1)],
                           Upper=summary_entry1$quantiles[c(5)])
    }else{
      fgh[[k]]<-data.frame(RR=summary_entry1$quantiles[3],
                           Lower=summary_entry1$quantiles[c(1)],
                           Upper=summary_entry1$quantiles[c(5)])
    }
    tau_square[[k]]<-data.frame(tau2=paste("0.00"),
                                Lower=paste("0.00"),
                                Upper=paste("0.00")
    )
    
  }
  
  ###To calculate direct estimate pooled
  #pooled_estimate<-as.numeric(c(t(all_direct_evidence2[i,c(3:5)])))
  if(nrow(temp11.1)==1){
    if(round(summary_entry1$statistics[1],1)==round(summary_entry1$quantiles[3],1)){
      pooled_estimate<-c(summary_entry1$statistics[1],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
      
    }else{
      pooled_estimate<-c(summary_entry1$quantiles[3],summary_entry1$quantiles[c(1)],summary_entry1$quantiles[c(5)])
      
    }
  }else{
    temp12<-temp11.1
    model.file<-"Rate_direct_multiplestudy_RE.txt"
    n.chains=3
    set.seed(4553)
    seeds <- sample.int(n.chains, n = .Machine$integer.max)
    model.inits <- list(
      list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
      list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
      list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]))
    jj.model<-jags.model(file=model.file,
                         data=list("E.t"=temp12$E.t,"x.t"=temp12$x.t,
                                   "x.c"=temp12$x.c,
                                   "E.c"=temp12$E.c,
                                   "NS"=length(temp12$NS)),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=500
    )
    update(object=jj.model,n.iter=10000,progress.bar='text')
    ## specify parameters to be monitored
    params <- c("RR")
    ## run JAGS and save posterior samples
    samps2[[i]] <- coda.samples( jj.model, params, n.iter=5000 )
    # par(mfrow = c(2,2))
    # gelman.plot(samps2[[i]])
    # dev.off()
    summary_entry<-summary(samps2[[i]])
    pooled_estimate<-c(summary_entry$statistics[1],summary_entry$quantiles[c(1)],summary_entry$quantiles[c(5)])
  }
  
  
  #Frequentist approach for I^2
  nma_data50_3<-Studies_contributing_to_direct_evidence_wide[order(as.numeric(gsub(".*? ", "", Studies_contributing_to_direct_evidence_wide$study_name))),]%>%
    select("study", "exposure.1", "responders.1","exposure.2", "responders.2" )
  nma_data50_3[,c(2,4)]<-nma_data50_3[,c(2,4)]/100
  nma_data50_4<-nma_data50_3
  colnames(nma_data50_4)<-c("NS","py.t","x.t","py.c","x.c")
  row.names(nma_data50_4)<-NULL
  m.inc<- metainc(event.e =x.t,
                  event.c = x.c,
                  time.e = py.t,
                  time.c = py.c,
                  data = nma_data50_4,
                  sm = "IRR",
                  method = "Inverse",
                  MH.exact = TRUE,
                  fixed = FALSE,
                  random = TRUE,
                  method.tau = "PM",
                  hakn = TRUE,
                  title = "SAE (final)"
  )
  
  
  summary_m.inc<-summary(m.inc)
  contribution_weight<-as.character(round((summary_m.inc$w.random/sum(summary_m.inc$w.random))*100,2))
  pooled_estimate2<-data.frame(RR=pooled_estimate[1],
                               Lower=pooled_estimate[2],
                               Upper=pooled_estimate[3])
  tau_square2<-data.frame(tau2=round( summary_m.inc$tau2,4)
                          ,
                          Lower=round(summary_m.inc$lower.tau2,4),
                          Upper=round(summary_m.inc$upper.tau2,4)
  )
  
  
  
  I_square_list[[i]]<-round((summary_m.inc$I2*100),2)
  I2<-paste("I^2 = ",paste(format(round((summary_m.inc$I2*100),2),nsmall=2),"% "),paste("(",format(round((summary_m.inc$lower.I2*100),2),nsmall=2),"%, ",format(round((summary_m.inc$upper.I2*100),2),nsmall=2),"%)",sep=""))
  
  
  
  
  
  RR_individual_studies<-do.call(rbind,fgh)
  temp11.2<-cbind(study_ordered,temp11.1)
  base_data <- tibble(mean  = RR_individual_studies$RR, 
                      lower = RR_individual_studies$Lower,
                      upper = RR_individual_studies$Upper,
                      study = temp11.2$study_ordered,
                      # treatment1 = paste0(as.character(temp11.2$x.t),"/",as.character(temp11.2$n.t)),
                      # treatment2 = paste0(as.character(temp11.2$x.c),"/",as.character(temp11.2$n.c)),
                      # treatment1 = paste0(as.character(temp11.2$E.c),"          ","(",as.character(temp11.2$x.c),")"),
                      # treatment2 = paste0(as.character(temp11.2$E.t),"          ","(",as.character(temp11.2$x.t),")"),
                      treatment1 = paste0(as.character(temp11.2$x.t),"          ","(",as.character(temp11.2$E.t),")"),
                      treatment2 = paste0(as.character(temp11.2$x.c),"          ","(",as.character(temp11.2$E.c),")"),
                      
                      RR = paste(format(round(RR_individual_studies$RR,2),nsmall=2),paste("(",format(round(RR_individual_studies$Lower,2),nsmall=2),", ",format(round(RR_individual_studies$Upper,2),nsmall=2),")",sep = "")),
                      ROB=ROB_ordered,
                      Indirectness=Indirectness_ordered,
                      Weight= contribution_weight
  )
  
  
  pooled <- tibble(mean  = pooled_estimate2$RR, 
                   lower = pooled_estimate2$Lower,
                   upper = pooled_estimate2$Upper,
                   study = "Summary",
                   treatment1=paste("Tau squared =",tau_square2$tau2,paste("(",tau_square2$Lower,", ",tau_square2$Upper,")",sep = "")),
                   #                   RR = paste(round(as.numeric(pooled_estimate2$RR),2),paste("(",round(as.numeric(pooled_estimate2$Lower),2),", ",round(as.numeric(pooled_estimate2$Upper),2),")",sep="")),
                   RR = paste(format(round(pooled_estimate2$RR,2),nsmall=2),paste("(",format(round(pooled_estimate2$Lower,2),nsmall=2),", ",format(round(pooled_estimate2$Upper,2),nsmall=2),")",sep="")),
                   summary = TRUE)
  
  header <- tibble(study = c("SAE (final)", "Study"),
                   # treatment1 = c("",selected_treatment_name(unique(temp11$treatment.1),subl)),
                   # treatment2 = c("",selected_treatment_name(unique(temp11$treatment.2),subl)),
                   treatment1 = c("Events (exposure patient-years)",selected_treatment_name(unique(temp11$treatment.1),mtxIR_SAE_subl)),
                   treatment2 = c("Events (exposure patient-years)",selected_treatment_name(unique(temp11$treatment.2),mtxIR_SAE_subl)),
                   #ROB=c(paste("Favours",selected_treatment_name(unique(temp11$treatment.1),subl)),"ROB"),
                   ROB=c("","ROB"),
                   Indirectness=c("","Indirectness"),
                   Weight=c("","Weight (%)"),
                   # OR = c(paste("Favours",selected_treatment_name(unique(temp11$treatment.2),subl)), "OR"),
                   RR = c("", "RR"),
                   summary = TRUE)
  empty_row <- tibble(mean = NA_real_)
  
  SAE_output_df <- bind_rows(header,
                             base_data,
                             empty_row,
                             pooled)
  
  
  
  
  SAE_output_df2<-rbind(SAE_output_df,SAE_output_df[c(nrow(SAE_output_df)),])
  SAE_output_df2[nrow(SAE_output_df),1]<-I2
  SAE_output_df2[nrow(SAE_output_df),2]<-NA
  SAE_output_df2[nrow(SAE_output_df),c(6:10)]<-c(NA,NA,NA,NA,NA)
  # acr50_output_df2[nrow(acr50_output_df),2]<-c("Prediction Interval")
  # acr50_output_df2[nrow(acr50_output_df),6]<-paste(format(round(prediction_interval$Predictive,2),nsmall=2),paste("(",format(round(prediction_interval$Lower,2),nsmall=2),", ",format(round(prediction_interval$Upper,2),nsmall=2),")",sep=""))
  # acr50_output_df2[nrow(acr50_output_df),8]<-c(round(prediction_interval$Predictive,2))
  # acr50_output_df2[nrow(acr50_output_df),9]<-c(round(prediction_interval$Lower,2))
  # acr50_output_df2[nrow(acr50_output_df),10]<-c(round(prediction_interval$Upper,2))
  
  SAE_output_df2_list[[i]]<-SAE_output_df2
  
  
  library(metafor)
  library(forestplot)
  graph_list[[i]]<-SAE_output_df2 %>% 
    forestplot(labeltext = c(study, treatment1, treatment2, RR,ROB,Indirectness,Weight), 
               is.summary = summary,
               graph.pos=5,
               #clip = c(0.1, 2.4), 
               xlog = T, 
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue"),
               # grid = structure(c(1.2), 
               #                  gp = gpar(lty = 2, col = "#F73831")),
               #txt_gp = fpTxtGp(cex=0.85),
               txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks = gpar(cex=0.6),cex=0.85),
               # xlab = paste("Outcome: ACR50-Direct Comparison",selected_treatment_name(unique(temp11$treatment.1),subl),"and",selected_treatment_name(unique(temp11$treatment.2),subl))
               xlab = paste("Favours",selected_treatment_name(unique(temp11$treatment.2),mtxIR_SAE_subl),"                  ","                    Favours",selected_treatment_name(unique(temp11$treatment.1),mtxIR_SAE_subl))
    )
  
}

length(graph_list)
graph_list[[3]]
n_studies<-c()
n_participants<-c()

countdown <- function(length) {
  count <- rev(pmax(length - 2:0, 0))
  return(count)
}


View(SAE_output_df2_list[[1]])

for(i in 1:length(graph_list)){
  df<-SAE_output_df2_list[[i]][-c(1,2,countdown(length(SAE_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  #Get the number of participants
  n_participants<-c(n_participants,sum(as.numeric(apply(df, 1, function(x) gsub("\\s*\\(.*", "", x)))))
}



Grade_assessment_dataset_direct_SAE<-data.frame(selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat1)),mtxIR_SAE_subl),
                                                selected_treatment_name(as.numeric(as.vector(all_direct_evidence2$treat2)),mtxIR_SAE_subl),
                                                Outcome="SAE",
                                                `Effect measure`="RR",
                                                `Number of studies`=n_studies,
                                                `Number of participant`=n_participants,
                                                ROB="",Indirectness="",Imprecision="",
                                                Inconsistency="",
                                                Publication_Bias="",Overall_Confidence="")

colnames(Grade_assessment_dataset_direct_SAE)[c(1,2,4,5,6,7,11,12)]<-c("Treatment",
                                                                       "Intevention",
                                                                       "Effect measure",
                                                                       "Number of studies",
                                                                       "Number of participant",
                                                                       "Risk of Bias",
                                                                       "Publication Bias",
                                                                       "Certainty of evidence (ignores imprecision)")

View(Grade_assessment_dataset_direct_SAE)
#Saving grade assessment details for pain mtxIR as an object of lists
grade_details_SAE_mtxIR<-list(graph_list,SAE_output_df2_list,Grade_assessment_dataset_direct_SAE,I_square_list)
#saveRDS(grade_details_HAQ_bio,file="List of HAQ Biologics data components for GRADE_OCT31.rds")
saveRDS(grade_details_SAE_mtxIR,file="List of SAE MTX-IR data components for GRADE_AUG2026.rds")

###
SAE_output_one_study_list <- list()
SAE_output_more_than_2_list <- list()
graph_list_one_study_list <- list()
graph_list_more_than_2_list <- list()
I_square_list_one_study_list <- list()
I_square_list_more_than_2_list <- list()

# Count studies for each comparison
n_studies<-c()
for(i in 1:length(graph_list)){
  df<-SAE_output_df2_list[[i]][-c(1,2,countdown(length(SAE_output_df2_list[[i]]$treatment1))),c(2,3)]
  #Get number of studies for each direct evidence
  n_studies<-c(n_studies,nrow(df))
  
}


# Indices
one_idx <- which(n_studies == 1)
more_idx <- which(n_studies >= 2)

# Split lists
SAE_output_one_study_list  <- SAE_output_df2_list[one_idx]
SAE_output_more_than_2_list <- SAE_output_df2_list[more_idx]

graph_list_one_study_list <- graph_list[one_idx]
graph_list_more_than_2_list <- graph_list[more_idx]

I_square_list_one_study_list <- I_square_list[one_idx]
I_square_list_more_than_2_list <- I_square_list[more_idx]

# Split GRADE datasets
Grade_assessment_dataset_direct_one_study <-
  Grade_assessment_dataset_direct_SAE[
    Grade_assessment_dataset_direct_SAE$`Number of studies` == 1,
  ]

Grade_assessment_dataset_direct_more_than_2 <-
  Grade_assessment_dataset_direct_SAE[
    Grade_assessment_dataset_direct_SAE$`Number of studies` != 1,
  ]

# Assemble outputs
grade_details_SAE_mtxIR_one_study <- list(
  graph_list_one_study_list,
  SAE_output_one_study_list,
  Grade_assessment_dataset_direct_one_study,
  I_square_list_one_study_list
)

grade_details_SAE_mtxIR_more_than_2 <- list(
  graph_list_more_than_2_list,
  SAE_output_more_than_2_list,
  Grade_assessment_dataset_direct_more_than_2,
  I_square_list_more_than_2_list
)

# Save
saveRDS(
  grade_details_SAE_mtxIR_one_study,
  file = "List of SAE MTX-IR data components for GRADE_one_study_AUG2026.rds"
)

saveRDS(
  grade_details_SAE_mtxIR_more_than_2,
  file = "List of SAE MTX-IR data components for GRADE_more_than_2_studies_AUG2026.rds"
)

###############################################################################################
###############################################################################################
inthreeArmtrial
inthreeArmtrial2names<-as.data.frame(apply(inthreeArmtrial,2,function(x)selected_treatment_name(x,mtxIR_SAE_subl)))
colnames(inthreeArmtrial2names)<-c("Treatment","Intevention")
Grade_assessment_dataset_direct_SAE
get_index_of_threeArm<-which(!is.na(match(paste(Grade_assessment_dataset_direct_SAE$Treatment,
                                                Grade_assessment_dataset_direct_SAE$Intevention),
                                          paste(inthreeArmtrial2names$Treatment, inthreeArmtrial2names$Intevention))))
datr_list<-list()
for(i in 1:length(get_index_of_threeArm)){
  datr<-as.data.frame(SAE_output_df2_list[get_index_of_threeArm[i]])
  datr_list[[i]]<-datr[datr$study=="Summary","RR"][2]
}


grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(3)]<-as.vector(do.call(rbind,datr_list))
grade_table4[paste(grade_table4$treat1,grade_table4$treat2,sep = ":")%in%paste(inthreeArmtrial2names$Treatment,inthreeArmtrial2names$Intevention,sep = ":"),c(4)]<-"NE"


####Get direct evidence report 
glen_direct_evidence_rating_sae<-read.csv("GLen_2026-09-03 23_16_SAE-MTXIR_complete.csv")
glen_direct_evidence_rating_sae2<-glen_direct_evidence_rating_sae[,-c(13:15)]
colnames(glen_direct_evidence_rating_sae2)[1:2]<-colnames(grade_table4)[1:2]

glen_direct_evidence_rating3<-merge(grade_table4[,c(1:3)],glen_direct_evidence_rating_sae2,
                                    by=c(colnames(grade_table4)[1:2]),all.x=T,all.y=T)
glen_direct_evidence_rating4<-cbind(glen_direct_evidence_rating3[,c(1,2,4:7)],
                                    glen_direct_evidence_rating3[,c(3,8:13)])
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Outcome),"Outcome"]<-"SAE"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Effect.measure),"Effect.measure"]<-"RR"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.studies),"Number.of.studies"]<-"Indirect evidence only"
glen_direct_evidence_rating4[is.na(glen_direct_evidence_rating4$Number.of.participant),"Number.of.participant"]<-"-"


glen_direct_evidence_rating5<-cbind(glen_direct_evidence_rating4,
                                    num1=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat1,mtxIR_SAE_subl)),
                                    num2=as.numeric(changing_treatment_name2num(glen_direct_evidence_rating4$treat2,mtxIR_SAE_subl)))
switch_func2<-function(x){
  if(which.max(x)==1){
    return(sort(x))
  }else{
    return(x)
  }
}


##ratings with comb
glen_direct_evidence_rating5[,14:15]<-t(apply(glen_direct_evidence_rating5[,14:15],1,switch_func2))
glen_direct_evidence_rating6<-cbind(glen_direct_evidence_rating5,
                                    comb=paste(glen_direct_evidence_rating5[,14],
                                               glen_direct_evidence_rating5[,15],sep = ":"))
glen_direct_evidence_rating6_onlyDirect<-glen_direct_evidence_rating6[!is.na(glen_direct_evidence_rating6$Certainty.of.evidence..ignores.imprecision.),]

##estimates with comb
grade_table1_with_comb<-cbind(grade_table1,num1=as.numeric(as.vector(grade_table1$treat1)),num2=as.numeric(as.vector(grade_table1$treat2)))
#grade_table1_with_comb<-grade_table1.1
grade_table1_with_comb[,12:13]<-t(apply(grade_table1_with_comb[,12:13],1,switch_func2))
View(grade_table1_with_comb)
grade_table1_with_comb1<-cbind(grade_table1_with_comb,
                               comb=paste(grade_table1_with_comb[,12],
                                          grade_table1_with_comb[,13],sep = ":"))
#grade_table1_with_comb1<-grade_table1_with_comb1[,-12]
View(grade_table1_with_comb1)
glen_direct_evidence_rating5$comb
grade_table1_with_comb1$comb

##Temporary GRADE look
Temp_GRADE<-merge(glen_direct_evidence_rating6[,-c(14,15)],grade_table1_with_comb1[,c(6:11,14)],
                  by=c(colnames(glen_direct_evidence_rating6)[16]),all.x=T,all.y=T
)

Temp_GRADE2<-cbind(Temp_GRADE[,c(2:14)],
                   Indirect=paste(Temp_GRADE[,15]," (",Temp_GRADE[,16],", ",Temp_GRADE[,17],")",sep = ""),
                   NMA=paste(Temp_GRADE[,18]," (",Temp_GRADE[,19],", ",Temp_GRADE[,20],")",sep = ""),
                   comb=Temp_GRADE[,c(1)]
)

Temp_GRADE3<-as.data.frame(t(apply(Temp_GRADE2, 1,replace_fun )))
Temp_GRADE3[is.na(Temp_GRADE3)]<-"-"
Temp_GRADE4<-cbind(Temp_GRADE3[,c(1:14)],QoE=" ",Temp_GRADE3[,15],QoE=" ",Temp_GRADE3[,16])
colnames(Temp_GRADE4)[c(15:18)]<-c("Quality of evidence","NMA","Quality of evidence","comb")
#Temp_GRADE4[Temp_GRADE4$comb%in%c("2:11","2:4"),"Indirect"]<-"NE"
write.csv(Temp_GRADE4,"Initial GRADE table for SAE MTX-IR_AUG2026.csv",row.names = F)
##Modified from Glens's suggestion
Temp_GRADE5<-cbind(Temp_GRADE4[,-c(17,18)],QoE=" ",Temp_GRADE4[,c(17,18)])
colnames(Temp_GRADE5)[17]<-"Quality of evidence without Imprecision"
write.csv(Temp_GRADE5,"Initial GRADE table for SAE MTX-IR_AUG2026 modified.csv",row.names = F)


#initial_gradeTable_haq<-read.csv("Initial GRADE table for HAQ modified.csv")
initial_gradeTable_sae<-Temp_GRADE5
estimate_contribution_matrix<-dat.net.contrib$random
all_nma_estimates<-row.names(estimate_contribution_matrix)

all_nma_estimates<-flip_to_smaller_to_larger(all_nma_estimates)


#ordering row and column names of contribution matrix
row.names(estimate_contribution_matrix)<-all_nma_estimates
colnames(estimate_contribution_matrix)<-flip_to_smaller_to_larger(colnames(estimate_contribution_matrix))


initial_gradeTable_sae2<-initial_gradeTable_sae[match(all_nma_estimates,initial_gradeTable_sae$comb),]

rating_scale<-data.frame(original_mednum=c(5:1),medname=c("HIGH","MODERATE","LOW","VERY LOW","EXTREMELY LOW"),
                         new_mednum=c(5:1))

#Creating contribution for the top two contributors, hence input number is 2
first_order_loop_contri_colnames <- generate_loop_labels(2)

indirect_rating<-list()

#Adding first order loop contribution
first_order_loop_contri<-as.data.frame(matrix("",nrow(initial_gradeTable_sae2),
                                              length(first_order_loop_contri_colnames)))

colnames(first_order_loop_contri)<-first_order_loop_contri_colnames
initial_gradeTable_sae2.1<-cbind(initial_gradeTable_sae2[,c(1:14)],first_order_loop_contri,
                                 initial_gradeTable_sae2[,-c(1:14)])

initial_gradeTable_sae2.1<-as.data.frame(contri_matrix_fxn2(initial_gradeTable_sae2.1,
                                                            all_nma_estimates,
                                                            estimate_contribution_matrix,
                                                            glen_direct_evidence_rating6_onlyDirect,
                                                            top_n = 2))

# Standardize the list to keep only the first column if multiple columns exist
standardized_list <- lapply(initial_gradeTable_sae2.1, function(x) {
  if (is.data.frame(x)) {
    return(x[, 1, drop = FALSE])  # Keep only the first column, ensure it's still a data frame
  } else {
    return(data.frame(x))         # Convert other types of list elements to data frame
  }
})
# Get original column names for first column of each list element
column_names <- names(sapply(standardized_list, function(x) colnames(x)[1]))

# Combine all elements of the list into a single data frame
result_df <- do.call(cbind, standardized_list)
result_df <-result_df[,-24] 
# Assign the original column names to the result data frame
colnames(result_df) <- column_names
initial_gradeTable_sae2.1<-result_df
initial_gradeTable_sae2.1$combI<-str_replace(initial_gradeTable_sae2.1$comb, ":", " vs ")

initial_gradeTable_sae2.1<-process_dataset(initial_gradeTable_sae2.1, 2)
initial_gradeTable_sae3<-initial_gradeTable_sae2.1[match(initial_gradeTable_sae$comb,initial_gradeTable_sae2.1$comb),]
estimate_contribution_matrix2 <- estimate_contribution_matrix
rescaled_estimate_contribution_matrix <- estimate_contribution_matrix2 / apply(estimate_contribution_matrix2, 1, sum) * 100
write.csv(rescaled_estimate_contribution_matrix,"Contribution matrix for SAE_MTXIR outcome_2026.csv",row.names = T)
#
###########Assessing Intransitivity
sub_data_intransitivity <- initial_gradeTable_sae3[initial_gradeTable_sae3$Indirect != "-", ]
top_n=2

sub_data_intransitivity2 <- sub_data_intransitivity
for (i in 1:top_n) {
  sub_data_intransitivity2[[paste0("dir_contri", i)]] <- gsub(" vs ", ":", sub_data_intransitivity[[paste0("First order loop ", i)]])
}

for (i in 1:nrow(sub_data_intransitivity2)) {
  direct_contri_indirectness_ratings <- initial_gradeTable_sae3[initial_gradeTable_sae3$comb %in% unlist(sub_data_intransitivity2[i, paste0("dir_contri", 1:top_n)]), ]
  
  intransitivity_vector <- check_vector(direct_contri_indirectness_ratings$Indirectness)
  
  rating_indir_num <- sum(
    as.numeric(changing_treatment_name2num(sub_data_intransitivity2$`Quality of evidence`[i], rating_scale)),
    as.numeric(intransitivity_vector[1])
  )
  rating_indir <- selected_treatment_name(rating_indir_num, rating_scale)
  
  sub_data_intransitivity2$Intransitivity[i] <- intransitivity_vector[2]
  sub_data_intransitivity2$`Certainty of evidence (ignores imprecision)-Indirect`[i] <- rating_indir
}



initial_gradeTable_sae3.1<-cbind(initial_gradeTable_sae3[,c(1:23)],
                                 Intransitivity="",
                                 `Certainty of evidence (ignores imprecision)-Indirect`="",
                                 initial_gradeTable_sae3[,-c(1:23)])
initial_gradeTable_sae3.1[initial_gradeTable_sae3.1$Indirect!="-",c("Intransitivity",
                                                                    "Certainty of evidence (ignores imprecision)-Indirect")]<-sub_data_intransitivity2[,c("Intransitivity",
                                                                                                                                                          "Certainty of evidence (ignores imprecision)-Indirect")]



initial_gradeTable_sae4<-initial_gradeTable_sae3.1


#######Creating percentage contribution of both direct and indirect evidence
####
colnames_vector<-colnames(estimate_contribution_matrix)
grade_temp<-cbind(initial_gradeTable_sae4[,c(1:25)],DirectPerContri="",IndirectPerContri="",initial_gradeTable_sae4[,-c(1:25)])
###Applying the P-values
P_values_nodesplit_saemtxIR2
grade_temp2<-cbind(grade_temp[,c(1:29)],p_value="-",grade_temp[,-c(1:29)])

drt<-grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),]
P_values_nodesplit_saemtxIR3 <- P_values_nodesplit_saemtxIR2[match(drt$comb, P_values_nodesplit_saemtxIR2$comb), , drop = FALSE]
grade_temp2[(grade_temp2$Direct!="-")&((grade_temp2$Indirect!="-")&(grade_temp2$Indirect!="NE")),"p_value"]<-round(P_values_nodesplit_saemtxIR3$p,2)

subSAESOF<-grade_temp2[grade_temp2$p_value!="-",]
subcontriMatSAE<-estimate_contribution_matrix[row.names(estimate_contribution_matrix)%in%subSAESOF$comb,]


# Initializing dirPer and IndirPer vectors
dirPer <- numeric(nrow(subcontriMatSAE))
IndirPer <- numeric(nrow(subcontriMatSAE))

# Loop to calculate dirPer and IndirPer
for (i in 1:nrow(subcontriMatSAE)) {
  dirPer[i] <- subcontriMatSAE[i, colnames(subcontriMatSAE) == row.names(subcontriMatSAE)[i]]
  IndirPer[i] <- 1- dirPer[i]
}

# Print the results
cat("dirPer:", dirPer, "\n")
cat("IndirPer:", IndirPer, "\n")

grade_temp2[match(row.names(subcontriMatSAE), grade_temp2$comb),"DirectPerContri"]<-round(dirPer,2)
grade_temp2[match(row.names(subcontriMatSAE), grade_temp2$comb),"IndirectPerContri"]<-round(IndirPer,2)



##Determining which evidence (direct or indirect) that contributes the most the NMA evidence
grade_temp2.1<-cbind(grade_temp2[,1:28],initialqoe="",grade_temp2[,-c(1:28)])
#When only direct or only indirect estimate is present
grade_temp2.1$initialqoe<-ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.!="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`==""),
                                 grade_temp2.1$Certainty.of.evidence..ignores.imprecision.,
                                 ifelse((grade_temp2.1$Certainty.of.evidence..ignores.imprecision.=="-")&(grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`!="-"),
                                        grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`,""))
#When both direct and indirect are present
fg<-c()
sub_grade_temp2.1<-grade_temp2.1[(grade_temp2.1$initialqoe==""),]
for(i in 1:nrow(sub_grade_temp2.1)){
  if(any(sub_grade_temp2.1[i,c(13,25)]=="Flag for manual review")){
    fg[i]<-"Flag for manual review"
  }else{
    if(length(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60))==0){ #Changing from the 70% to 60%
      if(which.max(c(as.numeric(changing_treatment_name2num(sub_grade_temp2.1$Certainty.of.evidence..ignores.imprecision.[i],rating_scale)),as.numeric(changing_treatment_name2num(sub_grade_temp2.1$`Certainty of evidence (ignores imprecision)-Indirect`[i],rating_scale))))==1
      ){
        fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
      }else{
        fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
        
      }
    }else if(which(apply(sub_grade_temp2.1[i,c(26,27)],1,as.numeric) > .60)==1){
      fg[i]<-sub_grade_temp2.1[i,"Certainty.of.evidence..ignores.imprecision."]
    }else{
      fg[i]<-sub_grade_temp2.1[i,"Certainty of evidence (ignores imprecision)-Indirect"]
    }
  }
  
  
}


grade_temp2.1[(grade_temp2.1$initialqoe==""),"initialqoe"]<-fg

####Evaluating evidence at the NMA phase
###########Assessing Incoherence
grade_temp2.2<-cbind(grade_temp2.1[,c(1:31)],Incoherence="",grade_temp2.1[,-c(1:31)])
incho<-c()
sub_grade_temp2.2<-grade_temp2.2[(grade_temp2.2$p_value!="-"),]
for(i in 1:nrow(sub_grade_temp2.2)){
  if(length(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60))==0){#Changing from the 70% to 60%
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Serious"
    }else if(as.numeric(sub_grade_temp2.2[i,"p_value"])>0.2){
      incho[i]<-"Not Serious"
      
    }else{
      incho[i]<-"Some Concerns"
    }
  }else if(which(apply(sub_grade_temp2.2[i,c(26,27)],1,as.numeric) > .60)==1){
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }else{
    if(as.numeric(sub_grade_temp2.2[i,"p_value"])<0.05){
      incho[i]<-"Some Concerns"
    }else{
      incho[i]<-"Not Serious"
    } 
  }
  
  
  
}


grade_temp2.2[(grade_temp2.2$p_value!="-"),"Incoherence"]<-incho


initial_gradeTable_sae5<-grade_temp2.2
initial_gradeTable_sae5<-initial_gradeTable_sae5[,-c(30,33)]
colnames(initial_gradeTable_sae5)[c(29,30)]<-c("Initial NMA evidence rating",
                                               "P-Value")

###After incoherence evaluation
initial_gradeTable_sae5.1<-cbind(initial_gradeTable_sae5[,c(1:31)],
                                 qoeAfterIncohenrence="",
                                 initial_gradeTable_sae5[,-c(1:31)])
initial_gradeTable_sae5.1[(initial_gradeTable_sae5.1$`P-Value`=="-"),"qoeAfterIncohenrence"]<-initial_gradeTable_sae5.1[(initial_gradeTable_sae5.1$`P-Value`=="-"),"Initial NMA evidence rating"]

sub_data_incoherence<-initial_gradeTable_sae5.1[initial_gradeTable_sae5.1$Incoherence!="",]

for(i in 1:nrow(sub_data_incoherence)){
  
  incoherence_vector<-check_vector(sub_data_incoherence[i,"Incoherence"])
  if(incoherence_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }else{
    rating_incho_num<-sum(as.numeric(changing_treatment_name2num(sub_data_incoherence$`Initial NMA evidence rating`[i],rating_scale)),
                          as.numeric(incoherence_vector[1]))
    rating_incho<-selected_treatment_name(rating_incho_num,rating_scale)
  }
  sub_data_incoherence$qoeAfterIncohenrence[i]<-rating_incho
}

initial_gradeTable_sae5.1[(initial_gradeTable_sae5.1$`P-Value`!="-"),"qoeAfterIncohenrence"]<- sub_data_incoherence$qoeAfterIncohenrence

#Replacing NMA estimates with direct estimates only when only direct evidence is present
initial_gradeTable_sae5.1[(initial_gradeTable_sae5.1$Direct!="-")&(initial_gradeTable_sae5.1$Indirect=="-"),"NMA"]<-initial_gradeTable_sae5.1[(initial_gradeTable_sae5.1$Direct!="-")&(initial_gradeTable_sae5.1$Indirect=="-"),"Direct"]


###Rating NMA after imprecision
initial_gradeTable_sae6<-cbind(initial_gradeTable_sae5.1,
                               p.estimate=apply(as.data.frame(initial_gradeTable_sae5.1$NMA),1,function(x){
                                 as.numeric(sub(" .*", "", x))
                               }),
                               t(apply(as.data.frame(initial_gradeTable_sae5.1$NMA),1,function(x){
                                 as.numeric(unlist(strsplit(gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]),",")))
                                 
                               })))

colnames(initial_gradeTable_sae6)[c(36,37)]<-c("l.bound","u.bound")




to_higher_lower_text <- function(e1) {
  sign_val <- sign(e1)
  abs_val <- abs(round(e1, 2))
  text <- ifelse(sign_val < 0, "lower", "higher")
  return(paste(abs_val, text))
}

####Running Rates NMA model
SAE_data<-current_SAE_dataset_mtxIR_withROB2_indrectness2
SAE_data$treatment<-renumber(SAE_data,"treatment")
SAE_data<-SAE_data[order(SAE_data$study),c(1:6)]


SAE_data_wide2<-reshape(number_of_arms(SAE_data), idvar=c("study","study_name"), timevar = "Arms", direction="wide")
row.names(SAE_data_wide2)<-NULL
xcv<-function(x){
  sum(!is.na(x))
}
numberOfArms<-apply(SAE_data_wide2[,c(3,7,11)],1,xcv)

SAE_data_wide3<-cbind(SAE_data_wide2,na=numberOfArms)
View(SAE_data_wide3)
colnames(SAE_data_wide3)<-c("study", "study_name", "t.1","e.1","n.1", "r.1",
                            "t.2","e.2","n.2","r.2","t.3","e.3","n.3","r.3", "na")
                            #"t.4","e.4","n.4","r.4","na" )
#######
####Running Rates NMA model
# load data for MTC
MTCData_rates <- as.data.frame(SAE_data_wide3)
r <- data.matrix(MTCData_rates[,c(6,10,14)])
n <- data.matrix(MTCData_rates[,c(5,9,13)])
t <- data.matrix(MTCData_rates[,c(3,7,11)])
E <- data.matrix(MTCData_rates[,c(4,8,12)])
na <- data.matrix(MTCData_rates[,"na"])
nt <- max(t, na.rm=TRUE)
ns <- nrow(r)

NTpbo <- length(t[,"t.1"][t[,"t.1"]==1])
xPBO <- r[t[,'t.1']==1,"r.1"]
nPBO <- n[t[,'t.1']==1,"n.1"]
E.PBO <- E[t[,'t.1']==1,"e.1"]

# define initial values 
# (other inits will be automatically generated when calling WinBUGS)

initv1 <- list(  d=c(NA,rep(0,nt-1)), mu=rep(0,ns),
                 .RNG.name="base::Super-Duper", .RNG.seed=1)

#loading rate model txt file
model.file4<-"Rate_nma_analysis.txt"
n.chains=3
model.inits <- initv1 

rates.model3<-jags.model(file=model.file4,
                         data= list("r"=r,"t"=t,"E"=E,"na" = na[,1], "nt" = nt, "ns" = ns,
                                    "NTpbo"= NTpbo,"xPBO"=xPBO,"E.PBO"=E.PBO),
                         n.chains=n.chains,
                         inits = model.inits,
                         n.adapt=5000
)

update(object=rates.model3,n.iter=30000,progress.bar='text')
## specify parameters to be monitored
params4 <- c("RR","illrisk")
## run JAGS and save posterior samples
samps4 <- coda.samples( rates.model3, params4, n.iter=10000 )
summary_samps<-summary(samps4)

#Baseline estimate
baseline_risk <-summary_samps$quantiles[which(rownames(summary_samps$quantiles)=="illrisk[1]"),3]
typeof(estimates_relative_to_baseline)

#Converting estimates relative to baseline to absolutes
#nma_estimate_bio_dataframe3<-as.data.frame(apply(t(do.call(rbind,nma_estimate_bio_dataframe2)),2,as.numeric))
estimates_relative_to_baseline<-initial_gradeTable_sae6[1:(nt-1),c(1,2,35:37)]
#estimates_relative_to_baseline<-initial_gradeTable_sae6[initial_gradeTable_sae6$treat2=="MTX",c(1,2,35:37)]
#absolutes_from_baseline<-cbind(baseline_risk,(estimates_relative_to_baseline[,3:5]*baseline_risk)/((1-baseline_risk)+(estimates_relative_to_baseline[,3:5]*baseline_risk)))
absolutes_from_baseline<-cbind(baseline_risk,1-exp(-estimates_relative_to_baseline[,3:5]*baseline_risk))
colnames(absolutes_from_baseline)[2:4]<-rep("",3)
estimates_relative_to_baseline_combine<-cbind(estimates_relative_to_baseline,absolutes_from_baseline)

#Getting all baseline values
all_baseline<-as.data.frame(c(baseline_risk,estimates_relative_to_baseline_combine[,7]))
#1-exp(-RR*ACR)
colnames(all_baseline)<-"base"
com.list<-list()
unique_treat2<-unique(initial_gradeTable_sae6$treat2)
for(i in 1:(nrow(estimates_relative_to_baseline_combine))){
  base<-all_baseline[i,]
  d.sub<-initial_gradeTable_sae6[initial_gradeTable_sae6$treat2==unique_treat2[i],c(1,2,35:37)]
  absolutes_conversion<-cbind(base,1-exp(-d.sub[,3:5]*base))
  #absolutes_conversion<-cbind(base,(d.sub[,3:5]*base)/((1-base)+(d.sub[,3:5]*base)))
  com.list[[i]]<-cbind(d.sub,absolutes_conversion)
}
com.dataframe<-do.call(rbind,com.list)
colnames(com.dataframe)[7:9]<-rep("",3)

# for(i in 4:11){
#   com.dataframe[com.dataframe$treat2==unique_treat2[i],"base"]<-all_baseline[i+1,]
# }


com.dataframe2<-cbind(com.dataframe,com.dataframe[,6:9]*1000)
more_fewer_text<-function(e1,e2){
  if(e1-e2>0){
    text<-"fewer"
  }else{
    text<-"more"
  }
  return(paste(round(abs(e1-e2),0),text))
}

more_fewer_signed<-function(e1,e2){
  return(round(e2-e1,0))
}

com.dataframe_withtext<-as.data.frame(cbind(com.dataframe,
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_signed,com.dataframe2[,10],com.dataframe2[,13]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,11]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,12]),
                                            mapply(more_fewer_text,com.dataframe2[,10],com.dataframe2[,13])))

prefix.assumedrisk <- ''
prefix.corrrisk<-''
suffix.assumedrisk <- ' per 1000 pt-years'
suffix.corrrisk <- ' per 1000 pt-years'
result_com.dataframe_withtext<-cbind(com.dataframe_withtext[,1:2],
                                     paste(prefix.assumedrisk,round(com.dataframe_withtext[,6]*1000,0), suffix.assumedrisk, sep=''),
                                     paste(prefix.corrrisk, round(com.dataframe_withtext[,7]*1000,0), suffix.corrrisk, ' (', round(com.dataframe_withtext[,8]*1000,0), ' to ',round(com.dataframe_withtext[,9]*1000,0),')', sep=''),
                                     paste(prefix.corrrisk, com.dataframe_withtext[,13], ' (', com.dataframe_withtext[,14], ' to ',com.dataframe_withtext[,15],')', sep=''),
                                     # paste(prefix.corrrisk, com.dataframe_withtext[,3], ' (', com.dataframe_withtext[,4], ' to ',com.dataframe_withtext[,5],')', sep=''),
                                     com.dataframe_withtext[,10:15]
)

colnames(result_com.dataframe_withtext)[c(3:11)]<-c("Assumed risk (comparator)",	"Corresponding risk (intervention)",
                                                    "Rate difference","p.estimate","l.bound","u.bound","","","")



initial_gradeTable_sae7<-cbind(initial_gradeTable_sae6[,-c(35:37)],
                               result_com.dataframe_withtext[,-c(1,2)],pair="")
#colnames(initial_gradeTable_sae7)[32]<-"Quality.of.evidence.after.Imprecision"

initial_gradeTable_sae8<-cbind(`Intervention number`=as.numeric(changing_treatment_name2num(initial_gradeTable_sae7$treat1,mtxIR_SAE_subl)),
                               `Comparator number`=as.numeric(changing_treatment_name2num(initial_gradeTable_sae7$treat2,mtxIR_SAE_subl)),
                               initial_gradeTable_sae7)
colnames(initial_gradeTable_sae8)[3:4]<-c("Intervention",	"Comparator")
initial_gradeTable_sae8<-initial_gradeTable_sae8[,-c(46)]
#write.csv(initial_gradeTable_ACR508,"ACR50 SoF table without imprecision assessment May2023.csv",row.names = F)
write.csv(initial_gradeTable_sae8,"SAE MTX-IR SoF table without imprecision assessment AUG2026.csv",row.names = F)

effect.size.rating_list<-lapply(split(initial_gradeTable_sae7, seq_len(nrow(initial_gradeTable_sae7))), function(x) {
  imprecision_fxn(x[,c("p.estimate","l.bound","u.bound")],MID = 50, large.effect = 100, framework = T)
})

as.numeric(do.call(rbind, effect.size.rating_list))

rating_scale_imprecision<-data.frame(original_mednum=c(0:-3),medname=c("Not Serious","Serious","Very Serious","Extremely Serious"),
                                     new_mednum=c(0:-3))

selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Example dataset rather than using imprecision application
initial_gradeTable_sae9<-initial_gradeTable_sae8
initial_gradeTable_sae9<-cbind(initial_gradeTable_sae9[,c(1:34,37:45,35:36)])
initial_gradeTable_sae9.1<-cbind(initial_gradeTable_sae9[,c(1:37)],
                                 Small.effect.threshold= 50,
                                 Large.effect.threshold= 100,
                                 Framework = "NULL",
                                 Imprecision.1 = "",
                                 Certainty.of.evidence = "",
                                 initial_gradeTable_sae9[,c(44,45)])


initial_gradeTable_sae9.1$Imprecision.1<-selected_treatment_name(as.numeric(do.call(rbind, effect.size.rating_list)),rating_scale_imprecision)

###Final certainty rating of the NMA
detailed_SOF_SAE_mtxIR<-initial_gradeTable_sae9.1
rating_scale2<-data.frame(original_mednum=c(0,0,-1,-2,-3),medname=c("Not Serious","Some Concerns","Serious","Very Serious","Extremely Serious"),
                          new_mednum=c(0,0,-1,-2,-3))
for(i in 1:nrow(detailed_SOF_SAE_mtxIR)){
  
  imprecision_vector<-check_vector(detailed_SOF_SAE_mtxIR[i,"Imprecision.1"])
  if(imprecision_vector[2]=="Some Concerns"){
    #rating_indir<-"Flag for manual review"
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_SAE_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_SAE_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }else{
    rating_impre_num<-sum(as.numeric(changing_treatment_name2num(detailed_SOF_SAE_mtxIR$qoeAfterIncohenrence[i],rating_scale)),
                          as.numeric(changing_treatment_name2num(detailed_SOF_SAE_mtxIR[i,"Imprecision.1"],rating_scale2)))
    rating_impre<-selected_treatment_name(rating_impre_num,rating_scale)
  }
  detailed_SOF_SAE_mtxIR$Certainty.of.evidence[i]<-rating_impre
}


View(detailed_SOF_SAE_mtxIR[is.na(detailed_SOF_SAE_mtxIR$Certainty.of.evidence),])
detailed_SOF_SAE_mtxIR[is.na(detailed_SOF_SAE_mtxIR$Certainty.of.evidence),"Certainty.of.evidence"]<-"EXTREMELY LOW"

#detailed_SOF_SAE_mtxIR$comb<-initial_gradeTable_sae7$comb

write.csv(detailed_SOF_SAE_mtxIR,"SAE_MTX-IR_Detail_SoF_OCT2026.csv", row.names = F)


# ###Combining MTX-IR and PBO stopped_mtxIR
# MTX_IR_PBO_stopped<-as.data.frame(read_xlsx("Combining MTX-IR and PBO stopped_mtxIR_dataset_21July2026.xlsx"))
# 
# library(dplyr)
# for(i in 1:nrow(MTX_IR_PBO_stopped)){
#   ppp<- master.medtable_dmard_ir %>% filter(MTX_IR_PBO_stopped$treatment[i] ==master.medtable_dmard_ir$medname )
#   MTX_IR_PBO_stopped$treatment[i]<-ppp[1,1]
# }
# MTX_IR_PBO_stopped.medlegend <- new.medtable(MTX_IR_PBO_stopped,'treatment')
# MTX_IR_PBO_stopped.medlegend
# 
# ##NMA plot reason for narrowing of network
# nma.plot(MTX_IR_PBO_stopped[,1],MTX_IR_PBO_stopped[,3],data = MTX_IR_PBO_stopped, trtname = as.vector(MTX_IR_PBO_stopped.medlegend$medname))
# nma.plot(MTX_IR_PBO_stopped[,1],MTX_IR_PBO_stopped[,3],data = MTX_IR_PBO_stopped)
# 
# #Narrow 1
# MTX_IR_PBO_stopped_narrow<-narrow.all(MTX_IR_PBO_stopped)
# MTX_IR_PBO_stopped_narrow
# 
# #sub network plot 1
# MTX_IR_PBO_stopped.d<-MTX_IR_PBO_stopped_narrow[[1]]
# length(unique(MTX_IR_PBO_stopped.d$treatment))
# MTX_IR_PBO_stopped_subl<-new.medtable(MTX_IR_PBO_stopped.d,'treatment')
# 
# nma.plot(MTX_IR_PBO_stopped.d[,1],MTX_IR_PBO_stopped.d[,3],data = MTX_IR_PBO_stopped.d, trtname = as.vector(MTX_IR_PBO_stopped_subl$medname), title = "ACR50 Network plot-Combining MTX-IR and PBO stopped_dmardIR")
# 
# 
# ##############################################################################
# #Performing sub network plot 1
# library(rjags)
# library(coda)
# library(gemtc)
# View(MTX_IR_PBO_stopped.d)
# MTX_IR_PBO_stopped_subl
# 
# # network_data_ACR50_dmardIR_1<-dmardIR.d%>%
# #   select(c("study", "study_name","treatment","responders", "sampleSize", "month",
# #            "OVERALL..ALL.OTHER.OUTCOMES", "Indirectness"))
# network_data_ACR50_MTX_IR_PBO_stopped_1<-MTX_IR_PBO_stopped.d%>%
#   select(c("study", "study_name","treatment","responders", "sampleSize", "month"))
# 
# network_data_ACR50_MTX_IR_PBO_stopped_1$treatment<-renumber(network_data_ACR50_MTX_IR_PBO_stopped_1,"treatment")
# ##To download acr50 data
# network_data_ACR50_MTX_IR_PBO_stopped_data_download<-network_data_ACR50_MTX_IR_PBO_stopped_1
# network_data_ACR50_MTX_IR_PBO_stopped_data_download$treatment<-selected_treatment_name(network_data_ACR50_MTX_IR_PBO_stopped_data_download$treatment,MTX_IR_PBO_stopped_subl)
# write.csv(network_data_ACR50_MTX_IR_PBO_stopped_data_download,"ACR50 MTX_IR_PBO_stopped dataset_JUNE2026.csv",row.names = F)
# # 1. Open jpeg file
# jpeg("ACR50 -MTX_IR_PBO_stopped Network plot 1 with named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_MTX_IR_PBO_stopped_1$study,network_data_ACR50_MTX_IR_PBO_stopped_1$treatment,data = network_data_ACR50_MTX_IR_PBO_stopped_1, trtname = as.vector(MTX_IR_PBO_stopped_subl$medname),title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# # 1. Open jpeg file
# jpeg("ACR50 -MTX_IR_PBO_stopped Network plot 1 without named interventions_JUNE2026.jpg", width = 3480, height = 2050, res=300)
# nma.plot(network_data_ACR50_MTX_IR_PBO_stopped_1$study,network_data_ACR50_MTX_IR_PBO_stopped_1$treatment,data = network_data_ACR50_MTX_IR_PBO_stopped_1,title = "ACR50 Network plot")
# # Close the pdf file
# dev.off() 
# 
# network_data_ACR50_MTX_IR_PBO_stopped_1<-mtc.network(network_data_ACR50_MTX_IR_PBO_stopped_1)
# 
# 
# ##RUN the model
# #mtc.hy.empirical.lor("semi-objective", "non-pharma")
# #mtc.hy.prior("std.dev", "dunif", 0, 2)
# #mtc.hy.prior("prec", "dgamma", 0.01, 0.01)
# set.seed(501)
# network_data_ACR50_MTX_IR_PBO_stopped_1.model <- mtc.model(network_data_ACR50_MTX_IR_PBO_stopped_1, likelihood='binom', link='logit',hy.prior=mtc.hy.empirical.lor("subjective", "pharma-pharma"))
# 
# # By default, the model will have 4 chains - generate a seed for each
# seeds <- sample.int(4, n = .Machine$integer.max)
# # Apply JAGS RNG settings to each chain
# network_data_ACR50_MTX_IR_PBO_stopped_1.model$inits<- mapply(c, network_data_ACR50_MTX_IR_PBO_stopped_1.model$inits, list(
#   list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seeds[1]),
#   list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=seeds[2]),
#   list(.RNG.name="base::Super-Duper", .RNG.seed=seeds[3]),
#   list(.RNG.name="base::Mersenne-Twister", .RNG.seed=seeds[4])), SIMPLIFY=FALSE)
# results_MTX_IR_PBO_stopped<-mtc.run(network_data_ACR50_MTX_IR_PBO_stopped_1.model,n.iter = 50000, n.adapt = 30000)
# summary(results_MTX_IR_PBO_stopped)
# results_MTX_IR_PBO_stopped.12<-relative.effect(results_MTX_IR_PBO_stopped,t1="12")
# 
# View(relative.effect.table.OR2(results_MTX_IR_PBO_stopped, comparator = "12"))
# refOR_MTX_IR_PBO_stopped<-as.data.frame(relative.effect.table.OR2(results_MTX_IR_PBO_stopped, comparator = "12"))
# 
# #library for mapvalues
# library(plyr)
# row.names(refOR_MTX_IR_PBO_stopped)<-mapvalues(row.names(refOR_MTX_IR_PBO_stopped),
#                                                from = c(1:(nrow(MTX_IR_PBO_stopped_subl))),
#                                                to = as.vector(MTX_IR_PBO_stopped_subl$medname))
# colnames(refOR_MTX_IR_PBO_stopped)<-mapvalues(colnames(refOR_MTX_IR_PBO_stopped),
#                                               from = c(1:(nrow(MTX_IR_PBO_stopped_subl))),
#                                               to = as.vector(MTX_IR_PBO_stopped_subl$medname))
# View(refOR_MTX_IR_PBO_stopped)
# # which(row(refOR)==col(refOR))
# # refOR[row(refOR)==col(refOR)]
# index_diagonal2<-data.frame(index1=1:nrow(refOR_MTX_IR_PBO_stopped),index2=1:nrow(refOR_MTX_IR_PBO_stopped))
# options(warn=-1)
# for(i in 1:nrow(refOR_MTX_IR_PBO_stopped)){
#   refOR_MTX_IR_PBO_stopped[index_diagonal2$index1[i],index_diagonal2$index2[i]]<-selected_treatment_name(as.numeric(refOR_MTX_IR_PBO_stopped[row(refOR_MTX_IR_PBO_stopped)==col(refOR_MTX_IR_PBO_stopped)])[i],MTX_IR_PBO_stopped_subl)
# }
# options(warn=0)
# View(refOR_MTX_IR_PBO_stopped)
# 
# write.csv(refOR_MTX_IR_PBO_stopped,"ACR50 -MTX_IR_PBO_stopped Network pairwise comparisons_JUNE2026.csv", row.names = T)
# # library(openxlsx)
# # write.xlsx(refOR_dmardIR,"Results for DMARD-IR_JUNE2026.xlsx",sheetName = "ACR50",rowNames = T)
# 
