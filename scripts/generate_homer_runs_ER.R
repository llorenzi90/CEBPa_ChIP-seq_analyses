## ---------------------------
##
##
## Purpose of script: generate homer runs for p30 vs p42 comparisons
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-12
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
# 1)3 modes run background related: 
#         A) against whole genome (i.e no background specified), 
#         B) against all peaks excluding own set, 
#         C) Against specific set (in case of p30, p42 and viceversa, in case of common, individual to one or the other -and in this case this is identical to B-) 
# 
# 2) Generate bed files for all conditions/backgrounds and 
# 
# First I need to do a liftover of all peaks because homer does not support mm39
# 
# 3)
# ##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
#~/homer/bin/findMotifsGenome.pl -h
# Program will find de novo and known motifs in regions in the genome
# 
# Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]
# Example: findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8

#for each set p30, p42, common
#create outdir

# create outdir specific background: background_allgenome,
#   background_all_other_peaks, background_antagonist_peaks

#
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/p30_vs_p42/bed_files/")
list_files=lapply(list.files(),read.table)
names(list_files)=gsub("_ER_UT_peaks.bed","",list.files())

head(list_files$p30)
#generate bed files with all combinations to use as backgrounds:
combs=combn(names(list_files),2)
combs
colnames(combs)=c("p30_and_common","p30_only_and_p42_only","p42_and_common")
for (col in colnames(combs)) {
  write.table(rbind(list_files[[combs[1,col]]],
                    list_files[[combs[2,col]]]),paste0(col,"_ER_UT_peaks.bed"),
              quote = F, col.names = F,row.names = F,sep = "\t")
  
}



##Liftover of all files
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/get_relative_path_Rfunction.R")
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/p30_vs_p42/")
dir.create("bed_files_lifted_over_mm10")
liftover_exe="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/liftOver"
liftover_exe=get_relative_path(liftover_exe)
#chain=import.chain(chain_file)
chain_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/liftover_chains/mm39ToMm10.over.chain"
chain_file=get_relative_path(chain_file)
for(fi in list.files("bed_files",full.names = T)){
  liftoverfname=gsub(".bed",".liftover_mm10.bed",basename(fi))
  liftoverunmapped=gsub(".bed",".liftover_mm10.unmapped.bed",basename(fi))
  
  cmm=paste(liftover_exe,fi,chain_file,paste0("bed_files_lifted_over_mm10/",liftoverfname),paste0("bed_files_lifted_over_mm10/",liftoverunmapped))

  system(cmm)
  }

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/homer/")
#example cmm: 
#"../../../../../scripts/liftOver bed_files/p42_ER_UT_peaks.bed ../../../../../references/mouse/liftover_chains/mm39ToMm10.over.chain bed_files_lifted_over_mm10/p42_ER_UT_peaks.liftover_mm10.bed bed_files_lifted_over_mm10/p42_ER_UT_peaks.liftover_mm10.unmapped.bed"
findmotifsgenome="/home/llorenzi/homer/bin/findMotifsGenome.pl"
bedfilesdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/p30_vs_p42/bed_files_lifted_over_mm10/"
bedfilesdir=get_relative_path(bedfilesdir)
patt="_ER_UT_peaks.liftover_mm10.bed"

in_list=list(p30=list(dname="p30",bgs=c("p42","p42_and_common")),
             p42=list(dname="p42",bgs=c("p30","p30_and_common")),
             common=list(dname="p30_p42",bgs=c("p30_only_and_p42_only")))
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/homer/")

for(set in names(in_list)){
  tli=in_list[[set]]
  inbed=paste0(bedfilesdir,"/",tli$dname,patt)
  sh_file=paste0(set,"_homer_commands.sh")
  dir.create(set)
  wd="background_all_genome"
  cmm1=paste(findmotifsgenome,inbed,"mm10",paste0(set,"/",wd))
  write(cmm1,sh_file)
  for (bg in tli$bgs) {
    wd=paste0("background_",bg)
    cmm2=paste(findmotifsgenome,inbed,"mm10",paste0(set,"/",wd),"-bg",paste0(bedfilesdir,"/",bg,patt))
    write(cmm2,sh_file,append = T)
  }
}
