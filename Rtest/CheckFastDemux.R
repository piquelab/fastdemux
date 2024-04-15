## 

library(tidyverse)

library(data.table)

md <- fread("../build/out.test.info.txt.gz")
colnames(md) <- c("Bcnum","Nsnp","Numi","BestScore","BestSample","SecondBestScore","SecondBestSample")
head(md)




dlda <- fread("../build/out.test.dlda.txt.gz")
dim(dlda)

bclist <- fread("../../scALOFT_2024/counts_cellranger_hg38/SCAIP10-PHA/raw_feature_bc_matrix/barcodes.tsv.gz",header=FALSE)
md$BARCODE <- bclist$V1
rm(bclist)
head(md)

demuxlet <- fread("/rs/rs_grp_scaloft/scALOFT_2024/counts_cellranger_hg38/demuxlet/demuxlet/SCAIP10-PHA.out.best")

head(demuxlet)

md2 <- inner_join(md,demuxlet)

table(md2$DROPLET.TYPE,md2$SNG.BEST.GUESS==md2$BestSample)

md2 <- md2 %>% 
  mutate(BestScore=BestScore/sqrt(Numi),
         SecondBestScore=SecondBestScore/sqrt(Numi),
         DiffScore=BestScore-SecondBestScore,
         DoubScore=sqrt(2)/2*(BestScore+SecondBestScore)-BestScore) 

md2 <- md2 %>% 
  filter(Nsnp>100,Numi>100)

dlda2 <- dlda[md2$Bcnum+1,]/sqrt(md2$Numi)

cs <- colMeans(dlda2)
hist(cs,breaks=50)

sum(cs>1.8)

table(md2$DROPLET.TYPE,md2$SNG.BEST.GUESS==md2$BestSample)

tt<- table(md2$BestSample[md2$DiffScore>5])
inbatch <- names(which(tt>500))

it <- which(tt>500)

table(md2$DiffScore>6 & md2$DoubScore< 0,md2$SNG.BEST.GUESS==md2$BestSample)


table(md2$DROPLET.TYPE,md2$BestScore>5 & md2$DoubScore>0)

table(md2$DROPLET.TYPE,md2$DiffScore>5 )


table(md2$DROPLET.TYPE,md2$DiffScore>5 & md2$DoubScore< 0)

table(md2$DROPLET.TYPE,md2$BestSample %in% inbatch)

table(md2$DiffScore>5 & md2$DoubScore<0 & md2$BestScore > 5,md2$BestSample %in% inbatch)

hist(unlist(dlda2[,10]),breaks=100)

md2$DiffScore[md2$DiffScore>6]=6;

ggplot(md2, aes(x=DiffScore,fill=DROPLET.TYPE)) + 
  geom_histogram() +
  theme_bw()

ggplot(md2, aes(x=DiffScore,y=BestScore,color=DROPLET.TYPE)) + 
  geom_point() +
  theme_bw()






