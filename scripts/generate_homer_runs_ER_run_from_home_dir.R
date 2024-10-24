## ---------------------------
##
##
## Purpose of script:generate homer run to run it from home (independent of internet connection)
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-13
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
## ---------------------------setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/homer/")
#example cmm: 
#"../../../../../scripts/liftOver bed_files/p42_ER_UT_peaks.bed ../../../../../references/mouse/liftover_chains/mm39ToMm10.over.chain bed_files_lifted_over_mm10/p42_ER_UT_peaks.liftover_mm10.bed bed_files_lifted_over_mm10/p42_ER_UT_peaks.liftover_mm10.unmapped.bed"
setwd("/home/llorenzi/tmp_run_homer/")
findmotifsgenome="/home/llorenzi/homer/bin/findMotifsGenome.pl"
bedfilesdir="bed_files_lifted_over_mm10"
#bedfilesdir=get_relative_path(bedfilesdir)
patt="_ER_UT_peaks.liftover_mm10.bed"

in_list=list(p30=list(dname="p30",bgs=c("p42","p42_and_common")),
             p42=list(dname="p42",bgs=c("p30","p30_and_common")),
             common=list(dname="p30_p42",bgs=c("p30_only_and_p42_only")))

for(set in names(in_list)){
  tli=in_list[[set]]
  inbed=paste0(bedfilesdir,"/",tli$dname,patt)
  sh_file=paste0(set,"_homer_commands.sh")
  dir.create(set)
  
  system(paste0("rm ",sh_file))
  for (bg in tli$bgs) {
    wd=paste0("background_",bg)
    cmm1=paste(findmotifsgenome,inbed,"mm10",paste0(set,"/",wd),"-bg",paste0(bedfilesdir,"/",bg,patt))
    write(cmm1,sh_file,append = T)
  }
  wd="background_all_genome"
  cmm2=paste(findmotifsgenome,inbed,"mm10",paste0(set,"/",wd))
  write(cmm2,sh_file,append = T)
}

