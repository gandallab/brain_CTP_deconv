# Julia Sidorenko
# 04/07/19
# This script uses dbSNP build 144 (dbSNP144.GRCh37) to rename genotyped SNPs on Illumina array to rsIDs

##PLEASE NOTE 1: the small insertions and deletions (I/D) are left unchanged in this pipeline 
##PLEASE NOTE 2: If the alleles for insertions and deletions (>1bp change) are other than "I" / "D" ,
# they won't be filtered out and can be replaced by SNPs (1bp change) at that specific position 

# UPDATES by Chloe:
# 191206: modified to avoid changes in bimfile order when chrX is added
# Correctly labelling chr23 and chr25 creates issues with ordering
# Preserve initial order by dealing with X somewhat separately when sorting files
# Shift chrX formatting before removing I/D (or else lose I/D variants on chrX)

print("Reading libraries")
library(plyr)
library(data.table)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

ref <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snpcount(ref)
seqlevels(ref)

print("Reading .bim file")
# Read in bim file
bimname = commandArgs(trailingOnly = TRUE)[1]
bimfile <- fread(bimname, colClasses = list(character = 1))

print("Deal with chrX")
## change 23 and 25 to X
## CAUTION: if there are other chr than 1~23,25, remove them use plink. Don't remove them only in bim file.
bimfile$V1<-ifelse(bimfile$V1=="23", "X",bimfile$V1)
bimfile$V1<-ifelse(bimfile$V1=="25", "X",bimfile$V1)
bimfile$V1<-ifelse(bimfile$V1=="XY", "X",bimfile$V1)

print("Separate the indels")
#separte the insertions/deletions
d<-subset(bimfile,V5=="D" | V5=="I" | V6=="D" | V6=="I" )
bim<-subset(bimfile,V5!="D" & V5!="I" & V6!="D" & V6!="I" )
bim$index <- 1:nrow(bim)
head(bim)
tail(bim)

print("Per chromosome, get SNP locations and match to bim file")
# For each chromosome download the SNP locations and match to bim file
a <- ddply(bim, .(V1), .progress="text", function(x)
{
  x <- mutate(x)
  chr <- x$V1[1]
  snps <- data.frame(snpsBySeqname(ref, chr))
  snps <- subset(snps, pos %in% x$V4, select=-c(alleles_as_ambig, strand))
  snps <- subset(snps, !duplicated(pos))
  snps <- subset(snps, !duplicated(RefSNP_id))
  x <- merge(x, snps, by.x="V4", by.y="pos", all.x=T)
  x <- x[order(x$index), ]
  index <- !is.na(x$RefSNP_id)
  x$V2[index] <- x$RefSNP_id[index]
  x <- subset(x, select=c(V1, V2, V3, V4, V5, V6))
  return(x)
})

print("Assign IDs if not present in chr:bp format")
# If the SNP has no ID, assign chr:bp12345 to avoid lots of variants with "."
miss <- which(a[,2] == ".")
a[miss,2] <- paste(a[miss,1], a[miss,4], sep = ":")

print("Label duplicates")
# If there are duplicated SNPs for any reason then label them as _2, _3 etc
temp <- rle(a$V2)
dup_out <- temp[which(duplicated(temp))]
temp2 <- paste0(rep(temp$values, times = temp$lengths), "_", unlist(lapply(temp$lengths, seq_len)))
temp2 <- gsub("_1", "", temp2)
a$V2 <- temp2

print("Clean output")
# Rename X and XY to 23
a[a[,1]=="X",1] <- 23
d[which(d[,1]=="X"),1] <- 23

# bind the insertions/deletions back
c <- rbind(a,d)

# add PAR boundaries boundaries 'b37'/'hg19': GRCh37/UCSC human genome 19, boundaries 2699520 and 154931044 
c$V1<-ifelse(c$V1==23 & c$V4<2699520, "25", c$V1)
c$V1<-ifelse(c$V1==23 & c$V4>154931044, "25", c$V1)

# reorder the SNPs
c_aut <- c[which(c$V1 %in% 1:22),]
c_x <- c[which(c$V1 %in% c(23,25)),]

sorted.aut <- c_aut[with(c_aut, order(as.numeric(V1),V4)),]
sorted.x <- c_x[which(c_x$V1 %in% c(23, 25)),]
sorted.x <- c_x[with(c_x, order(V4)),]
sorted.all <- rbind(sorted.aut, sorted.x)

#check the order
head(sorted.x)
tail(sorted.x)
length(which(bimfile$V4!=sorted.all$V4))

print("Write outputs")
# Save file
write.table(sorted.all, file=bimname, row=F, col=F, qu=F)
#write.table(c, file=paste(bimname, "_rsed_dbSNP144.GRCh37", sep="" ), row=F, col=F, qu=F)
write.table(dup_out, file=paste(bimname, "snpdup", sep="" ), row=F, col=F, qu=F)

