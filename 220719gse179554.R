#220719 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179554
#https://github.com/kpatel427/YouTubeTutorials/blob/main/dataManipulation.R 

# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# get metadata 
# get metadata --------
gse <- getGEO(GEO = 'GSE179554', GSEMatrix = TRUE)
# Error: The size of the connection buffer (131072) was not large enough                                          0s
# to fit a complete line:
#   * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

gse
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
# select, mutate, rename ------------
metadata.modified <- metadata %>%
  select(1,8,47,52) 
setwd("~/Desktop")
save(metadata.modified,file = '220719metadata_gse.rda')

  rename(tissue = source_name_ch1) %>%
  rename(condition = characteristics_ch1.2) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(condition = gsub("treatment: ", "", condition))

#https://mp.weixin.qq.com/s/uNvcIuhFrVY46Mh44GvLww
fs <- list.files('GSE179554_RAW/',
                 #pattern = '_Col_',
                 full.names = T)
fs
fs2 <- list.files('GSE179554_RAW/')
                 #pattern = '_Col_',
                 #full.names = T)
fs2
substring(fs2,1,10)
library(data.table)
fread(file.path('GSE179554_RAW/',fs[1]))


gid <- fread(fs[1],data.table = F)[,1]
# 1 v1 ENSMUSG00000000001 v2 Gnai3 v3 986
gid2 <- fread(fs[1],data.table = F)
#
dim(gid2)
colnames(gid2)
length(gid)
#54843
rawcount <- do.call(cbind,
                    lapply(fs,function(x){
                      fread(x,data.table = F)[,3]
                    }))

rawcount[1:4,1:4]
head(rawcount) 
dim(rawcount)
#54843 104 - > after keep, turned into 35866 104
#18977 were excluded

colnames(rawcount) <- substring(fs2,1,10)
rownames(rawcount) <- gid2$V2
save(rawcount,file = '220721rawcount.rda')
#without gene names error! error! error!
#read.table()
keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
#>=floor(0,75*ncol(rawcount))
table(keep)
length(keep)
View(keep)
dim(keep)
summary(keep)
rawcount2 <- rawcount[keep,]
#14870 * 104, 39973 were excluded
colnames(rawcount2) <- substring(fs2,1,10)
rawcount2[1:4,1:4]
head(rawcount2) 
dim(rawcount2)
rawcount2 <- as.data.frame(rawcount2)
save(rawcount2,file = '220723rawcount2.rda')

colnames(keep)
gid <- gid[keep]
gid2 <- gid2[keep]
#error! error! error! 
head(gid)
length(gid)
#from 54843 to 14870

exprset <- read.table(fs[1],header = T)[,1]
#error! error! error!

for (n in 2:104){
  rawcount=read.table(fs[n],header = T)
  exprset=cbind(exprset,rawcount[,3])
}
#no usage
{
library(AnnoProbe)
ids=annoGene(gid ,
              ID_type = 'ENSEMBL',
             species = 'mouse')
#0.34% of input IDs are fail to annotate... 
head(ids)
#14819 * 6 
#including "SYMBOL"   "biotypes" "ENSEMBL" 
#"chr"      "start"    "end"
colnames(ids)
ids=ids[!duplicated(ids$SYMBOL),]
rawcount=rawcount[match(ids$ENSEMBL, gid),]
rownames(rawcount) = ids$SYMBOL
}

rawcount2[1:4,1:4]
rawcount2['GAPDH',]
#Error in rawcount2["GAPDH", ] : subscript out of bounds
rawcount2['CXCR2',]
cg = log2(edgeR::cpm(rawcount)+1)['GAPDH',]
#error! error! error!
gene <- rownames(rawcount2)
'GAPDH' in gene
#Error: unexpected 'in' in "'GAPDH' in"
library(stringr)
str_detect('GAPDH',gene)
#error! error! error!
rawcount2 <- as.data.frame(rawcount2)
t(rawcount2)#error
rawcount3 <- data.frame(t(rawcount2))
exprset <- cbind(metadata.modified,rawcount3)
#104 * 14874 (14870+4)
save(exprset,file = 'exprset220724.rda')
#geneialized for all genes
colnames(exprset)
table(exprset$`drugs:ch1`)
exprset2 <- subset(exprset, `drugs:ch1` == 'vehicle' |
                     `drugs:ch1` == 'no treatment')

exprset_cxc <- exprset2[,c("title", "source_name_ch1", "drugs:ch1", "treatment:ch1",
                                 'Gapdh','Cxcl1','Cxcl2','Cxcr2')]
save(exprset_cxc,file = 'exprset_cxc220724.rda')
#specialized for cxc genes expression
table(exprset2$source_name_ch1)
exprset_liver <- subset(exprset2, source_name_ch1 == 'liver')
exprset_liver$Gapdh
exprset_liver$Cxcl1
exprset_liver$Cxcl2
exprset_liver$Cxcr2
exprset_lcxc <- exprset_liver[,c("title", "source_name_ch1", "drugs:ch1", "treatment:ch1",
                 'Gapdh','Cxcl1','Cxcl2','Cxcr2')]

library(edgeR)
cg <- log2(edgeR::cpm(exprset_lcxc[,-(1:4)])+1)[,'Gapdh']
gp <- exprset_lcxc$`treatment:ch1`
oneway.test(cg~gp)
#p-value = 0.3112

#https://statsandr.com/blog/how-to-do-a-t-test-or-anova-for-many-variables-at-once-in-r-and-communicate-the-results-in-a-better-way/
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("sham-operated", "CLP-6h"), 
                       c("sham-operated", "CLP-12h"), 
                       c("CLP-6h", "CLP-12h")) 
# comparisons for post-hoc tests
# Edit until here
exprset_lcxc2 <- log2(edgeR::cpm(exprset_lcxc[,-(1:4)])+1)
exprset_lcxc2 <- cbind(exprset_lcxc[,1:4],exprset_lcxc2)
head(exprset_lcxc2)
names(exprset_lcxc2)
# Edit from here
exprset_lcxc2$`treatment:ch1` <- factor(exprset_lcxc2$`treatment:ch1`,
                                        levels = c("sham-operated",
                                                   "CLP-6h",
                                                   "CLP-12h"))
save(exprset_lcxc2,file = 'exprset_lcxc2220724.rda')
x <- which(names(exprset_lcxc2) == "treatment:ch1") # name of grouping variable
y <- which(names(exprset_lcxc2) == "Gapdh" # names of variables to test
           | names(exprset_lcxc2) == "Cxcl1" |
             names(exprset_lcxc2) == "Cxcl2" |
             names(exprset_lcxc2) == "Cxcr2")
# Edit at your own risk
library(ggplot2)
library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(exprset_lcxc2,
                   x = colnames(exprset_lcxc2[j]), 
                   y = colnames(exprset_lcxc2[i]),
                   color = colnames(exprset_lcxc2[j]),
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., 
                                                ", p-value = ",
                                                ..p.format..)),
                             method = method1, 
                             label.y = max(exprset_lcxc2[, i], 
                                           na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, 
                           method = method2, 
                           label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}


exprset_cxc$`treatment:ch1` <- factor(exprset_cxc$`treatment:ch1`,
                                      levels = c("sham-operated",
                                                 "CLP-6h",
                                                 "CLP-12h"))
exprset_lung <- subset(exprset_cxc, source_name_ch1 == 'lung')

method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("sham-operated", "CLP-6h"), 
                       c("sham-operated", "CLP-12h"), 
                       c("CLP-6h", "CLP-12h")) 
# comparisons for post-hoc tests
# Edit until here
exprset_lung2 <- log2(edgeR::cpm(exprset_lung[,-(1:4)])+1)
exprset_lung2 <- cbind(exprset_lung[,1:4],exprset_lung2)
head(exprset_lung2)
names(exprset_lung2)
# Edit from here
exprset_lung2$`treatment:ch1` <- factor(exprset_lung2$`treatment:ch1`,
                                        levels = c("sham-operated",
                                                   "CLP-6h",
                                                   "CLP-12h"))
save(exprset_lung2,file = 'exprset_lung2_2220724.rda')
x <- which(names(exprset_lung2) == "treatment:ch1") # name of grouping variable
y <- which(names(exprset_lung2) == "Gapdh" # names of variables to test
           | names(exprset_lung2) == "Cxcl1" |
             names(exprset_lung2) == "Cxcl2" |
             names(exprset_lung2) == "Cxcr2")
# Edit at your own risk
library(ggplot2)
library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(exprset_lung2,
                   x = colnames(exprset_lung2[j]), 
                   y = colnames(exprset_lung2[i]),
                   color = colnames(exprset_lung2[j]),
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., 
                                                ", p-value = ",
                                                ..p.format..)),
                             method = method1, 
                             label.y = max(exprset_lung2[, i], 
                                           na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, 
                           method = method2, 
                           label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}

exprset_lung3 <- gather(exprset_lung2[,-(1:3)],
                        'gene','expression' ,-"treatment:ch1")
exprset_lung3$gene <- factor(exprset_lung3$gene,
                             levels = c("Gapdh","Cxcr2",
                                        "Cxcl1","Cxcl2"))
save(exprset_lung3,file = 'exprset_lung3_220724.rda')
    p <- ggboxplot(exprset_lung3,
                   "treatment:ch1", 
                   'expression',
                   color = 'gene',
                   #legend = "none",
                   palette = "npg",
                   #add = "jitter"
    )

p
    print(
      p + stat_compare_means(aes(label = paste0(..method.., 
                                                ", p-value = ",
                                                ..p.format..)),
                             method = method1, 
                             label.y = max(exprset_lung3[, i], 
                                           na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, 
                           method = method2, 
                           label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )

load("/Users/nanaliu/Desktop/exprset_lcxc2220724.rda")
exprset_lcxc2$`treatment:ch1`#right order / level
exprset_lcxc3 <- gather(exprset_lcxc2[,-(1:3)],
                        'gene','expression' ,-"treatment:ch1")
exprset_lcxc3$gene <- factor(exprset_lcxc3$gene,
                             levels = c("Gapdh","Cxcr2",
                                        "Cxcl1","Cxcl2"))
save(exprset_lcxc3,file = 'exprset_lcxc3_220724.rda')
p <- ggboxplot(exprset_lcxc3,
               "treatment:ch1", 
               'expression',
               color = 'gene',
               #legend = "none",
               palette = "npg",
               #add = "jitter"
)

p







