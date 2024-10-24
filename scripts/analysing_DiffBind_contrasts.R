## ---------------------------
##
##
## Purpose of script: process CSUC results of DiffBind CEBPa
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-18
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(DiffBind)
## ---------------------------

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/")

Mdba <- readRDS("DBA.keptlength.RDS")
# #Establishing a contrast
# Mdba <- dba.contrast(Mdba, categories=DBA_FACTOR, minMembers = 2)
# 
# #Performing the differential enrichment analysis
# dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# dba.show(dbObj, bContrasts=T)	
 Mdba$samples$Condition=factor(Mdba$samples$Condition,levels=c("p42","p30","EV"))
 Mdba$samples$Treatment=factor(Mdba$samples$Treatment,levels=c("UT","LPS"))
# 
Mdba2=Mdba
DBA <- dba.contrast(Mdba, design = '~ Condition * Treatment',reorderMeta=list(Condition=c("p42","p30","EV"),
                                                                              Treatment=c("UT","LPS")))
DBA$DESeq2

library(ExploreModelMatrix)

vd <- VisualizeDesign(sampleData = data.frame(Condition=Mdba$samples$Condition,
                                              Treatment=Mdba$samples$Treatment), 
                      designFormula =  ~ Condition * Treatment,
                      textSizeFitted = 4)
      
cowplot::plot_grid(plotlist = vd$plotlist)

#p42LPSvsp42UT
DBA1 <- dba.contrast(DBA, contrast = list("Treatment_LPS_vs_UT"))
dbObj <- dba.analyze(DBA)
dbObj1 <- dba.analyze(DBA1)
de1=dbObj1$contrasts[[1]]$DESeq2$de
de1$diff=ifelse(de1$padj<=0.05,ifelse(de1$fold>0,"UP","DOWN"),"NO")
#check if manually is the same:
DBA2 <- dba.contrast(DBA, group1 = Mdba$masks$p42&Mdba$masks$LPS,
                     group2 = Mdba$masks$p42&Mdba$masks$UT,
                     name1 = "p42LPS",
                     name2 = "p42UT")
dbObj2 <- dba.analyze(DBA2)
dba.show(dbObj, bContrasts=T)
de2=dbObj2$contrasts[[1]]$DESeq2$de
de2$diff=ifelse(de2$padj<=0.05,ifelse(de2$fold>0,"UP","DOWN"),"NO")
newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
Mdba$samples$Condition=newfeat #this does not work
Mdba$class[DBA_CONDITION,]=newfeat #this does


DBA <- dba.contrast(Mdba, 
                    design = '~ Condition')
DBA$DESeq2
vd <- VisualizeDesign(sampleData = data.frame(Condition=Mdba$samples$Condition),
                                              
                      designFormula =  ~ Condition ,
                      textSizeFitted = 4)

cowplot::plot_grid(plotlist = vd$plotlist)

#p42LPSvsp42UT with this kind of design:
DBA3 <- dba.contrast(DBA, contrast = c("Condition","p42LPS","p42UT"))
DBA4 <- dba.contrast(DBA, contrast = list("Condition_p42LPS_vs_EVLPS",
                                          "Condition_p42UT_vs_EVLPS"))

de4$diff=ifelse(de4$padj<=0.05,ifelse(de4$fold>0,"UP","DOWN"),"NO")
de3$diff=ifelse(de3$padj<=0.05,ifelse(de3$fold>0,"UP","DOWN"),"NO")

#These two ways are equivalent:
#p30LPSvsp42LPS:
#note DBA1 was generated with ~Condition*Treatment / (Condition=EV  EV  p30 p30 p30 p30 p30 p30 p42 p42 p42 p42)
DBA5= dba.contrast(DBA1, contrast = c(0,1,0,0,1,0))

#DBA was generated with ~Condition / Condition="EVLPS"  "EVUT"   "p30LPS" "p30LPS" "p30LPS" "p30UT" 
#"p30UT"  "p30UT"  "p42LPS" "p42LPS" "p42UT"  "p42UT" 
DBA6 <- dba.contrast(DBA, contrast = c("Condition","p30LPS","p42LPS"))
dbObj6
