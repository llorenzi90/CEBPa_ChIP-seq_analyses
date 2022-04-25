# sampleSheet data frame containing sample sheet, or file name of sample sheet to load (ignored
#                                                                                       if DBA is specified). Columns names in sample sheet may include:
# • SampleID: Identifier string for sample. Must be unique for each sample.
# • Tissue: Identifier string for tissue type
# • Factor: Identifier string for factor
# • Condition: Identifier string for condition
# • Treatment: Identifier string for treatment
# • Replicate: Replicate number of sample
# • bamReads: file path for bam file containing aligned reads for ChIP sample
# • bamControl: file path for bam file containing aligned reads for control
# sample
# dba 5
# • Spikein: file path for bam file containing aligned spike-in reads
# • ControlID: Identifier string for control sample
# • Peaks: path for file containing peaks for sample. Format determined by
# PeakCaller field or caller parameter
# • PeakCaller: Identifier string for peak caller used. If Peaks is not a bed
# file, this will determine how the Peaks file is parsed. If missing, will use
# default peak caller specified in caller parameter. Possible values:
#   – “raw”: text file file; peak score is in fourth column
# – “bed”: .bed file; peak score is in fifth column
# – “narrow”: default peak.format: narrowPeaks file
# – “macs”: MACS .xls file
# – “swembl”: SWEMBL .peaks file
# – “bayes”: bayesPeak file
# – “peakset”: peakset written out using pv.writepeakset
# – “fp4”: FindPeaks v4
# • PeakFormat: string indicating format for peak files; see PeakCaller and
# dba.peakset
# • ScoreCol: column in peak files that contains peak scores
# • LowerBetter: logical indicating that lower scores signify better peaks
# • Counts: file path for externally computed read counts; see dba.peakset
# (counts parameter)
# For sample sheets loaded from a file, the accepted formats are comma-separated
# values (column headers, followed by one line per sample), or Excel-formatted
# spreadsheets (.xls or .xlsx extension). Leading and trailing white space will
# be removed from all values, with a warning

example_sample_sheet <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/DS_AMKL/ChIP-seq/analyses/DiffBind/sample_sheet_Carini_csuc.csv")
countdata=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/merged_peaks_EV.p30.p42_LPS.UT.counts")

colnames(countdata) <- gsub(".sorted.markedDups.bam","",colnames(countdata))
colnames(countdata)

colda <- data.frame(SampleID = colnames(countdata))
colda$SampleID[colda$SampleID=="P42_LPS_R1"] <- "p42_LPS_R1"                    
colda$Condition <- gsub("(.*)_(.*)_(.*)","\\1",colda$SampleID)
#colda$Condition[colda$Condition=="P42"] <- "p42"
colda$Treatment=gsub("(.*)_(.*)_(.*)","\\2",colda$SampleID)
colda$Replicate=gsub("(.*)_(.*)_(.*)","\\3",colda$SampleID)

#bam files
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/")
list_bam_Files <- read.table("list_bam_files.txt")
list_bam_Files <- list_bam_Files$V1
sample_names_bam_Files <- gsub(".sorted.markedDups.bam","",basename(list_bam_Files))
sample_names_bam_Files[sample_names_bam_Files=="P42_LPS_R1"] <- "p42_LPS_R1"                    

table(colda$SampleID%in%sample_names_bam_Files)
colda$bamReads <- paste0("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/",
                         list_bam_Files[match(colda$SampleID,sample_names_bam_Files)])
inp_bams <- grep("INP",list_bam_Files,value = T)
inp_bams
names_inp_bams=gsub("INP_","",gsub(".sorted.markedDups.bam","",basename(inp_bams)))
table(colda$SampleID%in%names_inp_bams)
names_inp_bams[names_inp_bams=="P42_LPS_R1"] <- "p42_LPS_R1"
table(names_inp_bams%in%colda$SampleID)
table(names_inp_bams)
#I have to take new input file for sample INP_p30_LPS_R3
inp_bams[duplicated(names_inp_bams)]
inp_bams <- inp_bams[!duplicated(names_inp_bams)]
names_inp_bams <- names_inp_bams[!duplicated(names_inp_bams)]

colda$bamControl <- paste0("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/",
                           inp_bams[match(colda$SampleID,names_inp_bams)])
colda$ControlID <- gsub(".sorted.markedDups.bam","",basename(colda$bamControl))
colda$ControlID <- gsub("P42","p42",colda$ControlID)


colda$Peaks <- paste0(colda$bamReads,".NoModel_q_0.05_peaks.narrowPeak")
colda$PeakCaller="macs"

write.csv(colda,"sample_sheet_csuc.csv",row.names = F)
