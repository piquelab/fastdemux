## 

library(tidyverse)

library(data.table)

md <- fread("../build/out.SCAIP10-PHA.full.2024-04-19.info.txt.gz")


##dlda <- fread("../build/out.SCAIP10-PHA.full.2024-04-19.dlda.txt.gz")
##dim(dlda)

md$BARCODE.old <- md$BARCODE
bclist <- fread("../../scALOFT_2024/counts_cellranger_hg38/SCAIP10-PHA/raw_feature_bc_matrix/barcodes.tsv.gz",header=FALSE)
md$BARCODE <- bclist$V1[md$bcnum+1]

sum(md$BARCODE==md$BARCODE.old)

demuxlet <- fread("/rs/rs_grp_scaloft/scALOFT_2024/counts_cellranger_hg38/demuxlet/demuxlet/SCAIP10-PHA.out.best")

head(demuxlet)

md2 <- inner_join(md,demuxlet)

table(md2$DROPLET.TYPE,md2$SNG.BEST.GUESS==md2$BestSample)

md2 <- md2 %>% 
  filter(Nsnp>100,Numi>100)

table(md2$DROPLET.TYPE,md2$SNG.BEST.GUESS==md2$BestSample)


##dlda2 <- dlda[md2$Bcnum+1,]/sqrt(md2$Numi)
##cs <- colMeans(dlda2)
##hist(cs,breaks=50)
##sum(cs>1.8)

table(md2$DROPLET.TYPE,md2$SNG.BEST.GUESS==md2$BestSample)

tt<- table(md2$BestSample[md2$DropType==1])
inbatch <- names(which(tt>200))

it <- which(tt>500)

table(md2$DropType==1 ,md2$SNG.BEST.GUESS==md2$BestSample)

table(md2$DropType,md2$DROPLET.TYPE)

tt <- table(md2$DROPLET.TYPE,md2$BestSample %in% inbatch)
tt
tt/rowSums(tt)


tt <- table(md2$DropType==1,md2$BestSample %in% inbatch)
tt
tt[2,]/sum(tt[2,])

tt <-table(md2$DropType,md2$BestSample %in% inbatch)
tt
tt/rowSums(tt)


hist(log10(md2$BestScore),breaks=100)

ggplot(md2, aes(x=log10(BestScore),fill=DROPLET.TYPE)) + 
  geom_histogram() +
  theme_bw()






