setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/")
list.files("analyses/DiffBind/")
sample_sheet=read.csv("analyses/DiffBind/sample_sheet_csuc.csv")
sequenced_pairs=read.table("data/total_sequenced_pairs_H3K27ac.txt",sep = "\t",header = T)

sample_sheet$Sequenced_pairs=sequenced_pairs$Total.reads[match(tolower(sample_sheet$SampleID),
                                                               tolower(sequenced_pairs$Sample))]
sequenced_pairs_input=sequenced_pairs[sequenced_pairs$Type.of.data=="input",]
sequenced_pairs_input$Sample <- gsub("INP_","",sequenced_pairs_input$Sample)
sample_sheet$Input_Sequenced_pairs=sequenced_pairs_input$Total.reads[match(tolower(sample_sheet$SampleID),
                                                               tolower(sequenced_pairs_input$Sample))]

countdata=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/merged_peaks_EV.p30.p42_LPS.UT.counts")
colnames(countdata) <- gsub(".sorted.markedDups.bam","",colnames(countdata))
reads_in_peaks=colSums(countdata)
sample_sheet$Pair_reads_in_peks=reads_in_peaks[match(tolower(sample_sheet$SampleID),
                                                     tolower(names(reads_in_peaks)))]
sample_sheet$FRiP=sample_sheet$Pair_reads_in_peks/sample_sheet$Sequenced_pairs
write.csv(sample_sheet[,c(1:4,10:12)],"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/sample_sheet_H3K27ac.csv",row.names = F)

