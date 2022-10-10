library(GGally)
library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(biomaRt)
library(limma)
library(VennDiagram)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)

###Data processing for downstream analysis
noncoding <- openxlsx::read.xlsx("~/Project/WJN/noncoding.xlsx") %>% 
  .[grep("lcl",.$Accession),]
colnames(noncoding) <- gsub("[Abundances.(Normalized):.]","",colnames(noncoding))
rownames(noncoding) <- noncoding[,1]
noncoding <- noncoding[,-1]
for (i in colnames(noncoding)) {
  noncoding[,i] <- as.numeric(noncoding[,i])
}

colnames(noncoding) <- gsub("w1","Hup-T3",colnames(noncoding)) %>% 
  gsub("w2","Hup-T4",.) %>% gsub("w3","Patu8988T",.) %>% 
  gsub("w4","QGP-1",.) %>% gsub("w5","Panc1",.)

noncoding$na_num <- apply(noncoding, 1, function(x){table(is.na(x))["TRUE"]}) %>% ifelse(is.na(.),0,.)
noncoding <- noncoding[noncoding$na_num != 15,]

####caculate peptide expression correlation of different cell lines
data_list <- c("Hup-T3","Hup-T4","Patu8988T","QGP-1","Panc1")
names(data_list) <- data_list
data_list <- lapply(data_list, function(x){
  noncoding[,grep(x,colnames(noncoding))]
})

pdf("~/Project/WJN/fig3b.pdf",height = 6,width = 7)
for (i in names(data_list)) {
  p <- ggpairs(data_list[[i]],
               lower = list(continuous = wrap(ggscatter), method="pearson"),
               diag = list(continuous = wrap(ggdehist)))+
    theme_bw()+theme(panel.grid = element_blank())
  print(p)
}
dev.off()

###count of different initiate codon type

data_list <- lapply(data_list, function(x){
  x$na_num <- apply(x, 1, function(xx){table(is.na(xx))["TRUE"]}) %>% ifelse(is.na(.),0,.)
  x <- x[x$na_num < 2,]
  x <- x[,-4]
  x$Mean <- apply(x,1,function(xx)mean(xx,na.rm=T))
  x$ORF_ID <- rownames(x)
  return(x)
})

noncoding_sub <- data_list[[1]][,c(1:3,5)]
for (i in 2:5) {
  noncoding_sub <- merge(noncoding_sub,data_list[[i]][,c(1:3,5)],by="ORF_ID",all=T)
}

openxlsx::write.xlsx(noncoding_sub,"~/Project/WJN/filtered_exp.xlsx")

rnacentral <- readBStringSet("~/Project/ProDatabase/reference/RNAcentral/RNAcentral_homo_sapiens_GRCh38.fa",
                             format = "fasta",nrec = -1L,skip=0L, seek.first.rec=TRUE, use.names=TRUE) %>% 
  as.data.frame()
rnacentral$RNAcentraID <- rownames(rnacentral)
colnames(rnacentral)[1] <- "RNA_Sequence"

seqdata <- readBStringSet("~/Project/ProDatabase/reference/RNAcentral_CCLE_uniprot/combine_CCLE_cpat_cpc_noncoding_RNAcentral_uniprot.fasta",
                          format = "fasta",nrec = -1L,skip=0L, seek.first.rec=TRUE, use.names=TRUE) %>% as.data.frame()
rownames(seqdata) <- rownames(seqdata) %>% gsub("[|]","_",.) %>% 
  gsub(" unnamed protein product","",.) %>% gsub(", partial","",.)
seqdata_sub <- seqdata[rownames(seqdata) %in% noncoding_sub$ORF_ID,,drop=F]

seqdata_sub$ORF_ID <- rownames(seqdata_sub)
colnames(seqdata_sub)[1] <- "Pep_sequence"

first_char <- lapply(gregexpr("[:]",rownames(seqdata_sub)), function(x){(x[1])+1}) %>% unlist()
sec_char <- lapply(gregexpr("[:]",rownames(seqdata_sub)), function(x){(x[2])-1}) %>% unlist()
seqdata_sub$start <- substr(rownames(seqdata_sub),first_char,sec_char) %>% as.numeric()
seqdata_sub$start <- (seqdata_sub$start)+1
seqdata_sub$end <- substring(rownames(seqdata_sub),sec_char+2) %>% as.numeric()
seqdata_sub$end <- (seqdata_sub$end)+1
seqdata_sub$RNAcentraID <- substr(rownames(seqdata_sub),
                                  regexpr("URS",rownames(seqdata_sub)),
                                  regexpr(":",rownames(seqdata_sub))-1)
seqdata_sub$ORF_ID <- rownames(seqdata_sub)
colnames(seqdata_sub)[1] <- "Pep_sequence"

seqdata_sub <- merge(seqdata_sub,rnacentral,by="RNAcentraID")
seqdata_sub$ini_codon <- substr(seqdata_sub$RNA_Sequence,seqdata_sub$start,(seqdata_sub$start)+2)

ini_codon <- as.data.frame(table(seqdata_sub$ini_codon))
colnames(ini_codon) <- c("ini_codon","Num")
ini_codon <- ini_codon[order(ini_codon$Num),]
ini_codon$ini_codon <- gsub("T","U",ini_codon$ini_codon)
myLabel = as.vector(ini_codon$ini_codon)
myLabel = paste0(myLabel,"(",ini_codon$Num,", ",round(ini_codon$Num/sum(ini_codon$Num) * 100,2),"%)")
myLabel = myLabel[c(grep("UUG",myLabel),grep("CUG",myLabel),grep("AUG",myLabel))]
openxlsx::write.xlsx(ini_codon,"~/Project/WJN/ini_codong_number_total.xlsx")

pdf("~/Project/WJN/fig4a.pdf",height = 6,width = 5)
ggplot(ini_codon, aes(x = "", y = Num, fill = ini_codon)) +
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid=element_blank(),
        panel.background = element_blank()) + 
  scale_fill_discrete(breaks = ini_codon$ini_codon, labels = myLabel) +
  geom_text(aes(y = Num/2 + c(0, cumsum(Num)[-length(Num)]), x = 1, label = myLabel), size = 5)
dev.off()

###Peptide length
seqdata_sub$Pep_length <- nchar(seqdata_sub$Pep_sequence)
count_stat <- as.data.frame(table(cut(seqdata_sub$Pep_length,breaks = c(0,25,50,75,100,max(seqdata_sub$Pep_length)))))
names(count_stat) <- c("PEP_length","count")
count_stat$PEP_length <- factor(x=c("0-25","25-50","50-75","75-100",">100"),
                                levels = c("0-25","25-50","50-75","75-100",">100"))
openxlsx::write.xlsx(count_stat,"~/Project/WJN/Pep_length_Freq.xlsx")
pdf("~/Project/WJN/fig4b1.pdf",height = 5,width = 6)
ggplot(data = count_stat,mapping = aes(x=PEP_length,y=count,color=PEP_length,fill=PEP_length))+
  geom_bar(stat = "identity") + 
  theme_classic()+
  geom_text(aes(y = count-2, x = PEP_length, label = count),color="black", size = 5)
dev.off()

####Expression of peptide with different length

rownames(noncoding_sub) <- noncoding_sub$ORF_ID
noncoding_sub <- noncoding_sub[,-1]
noncoding_sub$Mean <- apply(noncoding_sub,1,function(x)mean(x,na.rm=T)) %>% log2

Mean <- data.frame(ORF_ID=rownames(noncoding_sub),Mean=noncoding_sub$Mean)
Mean <- left_join(Mean,seqdata_sub,by="ORF_ID")
Mean <- mutate(Mean,seq_length=nchar(Pep_sequence))
Mean <- mutate(Mean,length_section=cut(seq_length,scientific = FALSE,
                                       breaks = c(0,25,50,75,100,max(seq_length))))
Mean$length_section <- Mean$length_section  %>% 
  gsub("[(]","",.) %>% 
  gsub("[]]","",.) %>% 
  gsub("0,25","0-25",.) %>% 
  gsub("25,50","25-50",.) %>% 
  gsub("50,75","50-75",.)%>% 
  gsub("75,100","75-100",.)%>% 
  gsub("[+]","",.) %>%
  gsub("100,2.56e03",">100",.)
Mean$length_section <- factor(Mean$length_section,
                              levels = c("0-25","25-50","50-75","75-100",">100"))

compaired_length <- list(c("0-25",">100"),c("25-50",">100"),
                         c("50-75",">100"),c("75-100",">100"))

pdf("~/Project/WJN/fig4c.pdf",height = 9,width = 8)
ggplot(Mean,aes(x=length_section,y=Mean,color=length_section))+
  geom_boxplot()+geom_jitter(position = position_jitterdodge()) +
  theme_bw() + 
  ylab("log2(Peptiede Expression)")+
  theme(axis.text.x=element_text(size=15), 
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=16), 
        axis.title.y=element_text(size = 20),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1), 
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_signif(comparisons = compaired_length,
              step_increase = 0.1,
              map_signif_level = F,
              test = "t.test",
              color="black")
dev.off()



###count of different codon type in different peptide length section 

ini_codon_section <- as.data.frame(table(Mean$length_section,Mean$ini_codon))
colnames(ini_codon_section) <- c("length_section","ini_codon","Num")
ini_codon_section$ini_codon <- gsub("T","U",ini_codon_section$ini_codon)

length_section <- as.character(unique(ini_codon_section$length_section))
names(length_section) <- length_section
length_section <- lapply(length_section, function(x){
  ini_codon_section[ini_codon_section$length_section %in% x,]
})

openxlsx::write.xlsx(length_section,"~/Project/WJN/ini_codong_number_length_section.xlsx")

pdf("~/Project/WJN/fig4b2.pdf",height = 6,width = 5)
for (i in c("0-25","25-50","50-75","75-100",">100")) {
  length_section[[i]] <- length_section[[i]][order(length_section[[i]]$Num),]
  myLabel = as.vector(length_section[[i]]$ini_codon)
  myLabel = paste0(myLabel,"(",length_section[[i]]$Num,", ",round(length_section[[i]]$Num/sum(length_section[[i]]$Num) * 100,2),"%)")
  myLabel = myLabel[c(grep("UUG",myLabel),grep("CUG",myLabel),grep("AUG",myLabel))]
  p <- ggplot(length_section[[i]], aes(x = "", y = Num, fill = ini_codon)) +
    geom_bar(stat = "identity", width = 1) +  
    coord_polar(theta = "y") + 
    labs(x = "", y = "", title = "") + 
    theme(axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          panel.grid=element_blank(),
          panel.background = element_blank()) + 
    scale_fill_discrete(breaks = length_section[[i]]$ini_codon, labels = myLabel) +
    geom_text(aes(y = Num/2 + c(0, cumsum(Num)[-length(Num)]), x = 1, label = myLabel), size = 5)+
    ggtitle(i)
  print(p)
}
dev.off()
