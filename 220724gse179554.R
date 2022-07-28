#
load("/Users/nanaliu/Desktop/exprset_cxc220724.rda")
exprset_cxc2 <- log2(edgeR::cpm(exprset_cxc[,-(1:4)])+1)
exprset_cxc2 <- cbind(exprset_cxc[,c(2,4)],exprset_cxc2)
save(exprset_cxc2,file = 'exprset_cxc2_220724.rda')

exprset_cxc3 <- exprset_cxc2[-c(49:60),]
#to exclude ex_vivo samples
exprset_cxc3$`treatment:ch1` <- factor(exprset_cxc3$`treatment:ch1`,
                                        levels = c("sham-operated",
                                                   "CLP-6h",
                                                   "CLP-12h"))
exprset_cxc3$`treatment:ch1`

library(dplyr)
library(tidyverse)
exprset_cxc4 <- gather(exprset_cxc3,
                        'gene','expression' ,
                       -c("source_name_ch1","treatment:ch1"))
exprset_cxc4$gene <- factor(exprset_cxc4$gene,
                             levels = c("Gapdh","Cxcr2",
                                        "Cxcl1","Cxcl2"))
save(exprset_cxc4,file = 'exprset_cxc4_220724.rda')
#long format of cxc gene expression with all organs expected for in vitro samples
library(ggpubr)
p <- ggboxplot(exprset_cxc4,
               "treatment:ch1", 
               'expression',
               color = 'gene',
               #legend = "none",
               palette = "npg",
               facet.by = 'source_name_ch1',
               #nrow=1,
               add = "point",
               xlab = 'Treatment',
               ylab = 'Expression'
)
p
#for all organs
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("sham-operated", "CLP-6h"), 
                       c("sham-operated", "CLP-12h"), 
                       c("CLP-6h", "CLP-12h")) 
# comparisons for post-hoc tests
x <- which(names(exprset_cxc3) == "treatment:ch1") # name of grouping variable
y <- which(names(exprset_cxc3) == "Gapdh" # names of variables to test
           | names(exprset_cxc3) == "Cxcr2" |
             names(exprset_cxc3) == "Cxcl2" |
             names(exprset_cxc3) == "Cxcl1")
# Edit at your own risk
library(ggplot2)
library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(exprset_cxc3,
                   x = colnames(exprset_cxc3[j]), 
                   y = colnames(exprset_cxc3[i]),
                   color = colnames(exprset_cxc3[j]),
                   facet.by = 'source_name_ch1',
                   legend = "none",
                   palette = "npg",
                   add = "point",
                   xlab = 'Treatment'
                   )
    print(p + stat_compare_means(aes(label = paste0(..method.., 
                                                ", p-value = ",
                                                ..p.format..)),
                             method = method1, 
                             size=3,
                             label.x = 2.4,
                             label.y = min(exprset_cxc3[, i], 
                                           na.rm = TRUE))+ #
        stat_compare_means(comparisons = my_comparisons, 
                           method = method2, 
                           size=3,
                           label = "p.format"))
  }
}

load("/Users/nanaliu/Desktop/gse17955/exprset_cxc2_220724.rda")
exprset_cxc5 <- exprset_cxc2[c(49:60),]
#to exclude ex_vivo samples
exprset_cxc5$`treatment:ch1` <- factor(exprset_cxc5$`treatment:ch1`,
                                       levels = c("Control",
                                                  "LPS 1h",
                                                  "LPS 4h",
                                                  "LPS 24h"))
exprset_cxc5$`treatment:ch1`
library(dplyr)
library(tidyverse)
exprset_cxc6 <- gather(exprset_cxc5,
                       'gene','expression' ,
                       -c("source_name_ch1","treatment:ch1"))
exprset_cxc6$`treatment:ch1`
exprset_cxc6$gene <- factor(exprset_cxc6$gene,
                            levels = c("Gapdh","Cxcl1",
                                       "Cxcl2","Cxcr2"))
save(exprset_cxc6,file = 'exprset_cxc6_220724.rda')

library(ggpubr)
p <- ggboxplot(exprset_cxc6,
               "treatment:ch1", 
               'expression',
               color = 'gene',
               #legend = "none",
               palette = "npg",
               #facet.by = 'source_name_ch1',
               #nrow=1,
               add = "point",
               xlab = 'Treatment',
               ylab = 'Expression'
)
p
#for all organs
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Control","LPS 1h"),
                       c("Control","LPS 4h"),
                       c("Control","LPS 24h"),
                       c("LPS 1h","LPS 4h"), 
                       c("LPS 1h","LPS 24h"), 
                       c("LPS 4h","LPS 24h")) 
# comparisons for post-hoc tests
x <- which(names(exprset_cxc5) == "treatment:ch1") # name of grouping variable
y <- which(names(exprset_cxc5) == "Gapdh" # names of variables to test
           | names(exprset_cxc5) == "Cxcr2" |
             names(exprset_cxc5) == "Cxcl2" |
             names(exprset_cxc5) == "Cxcl1")
# Edit at your own risk
library(ggplot2)
library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(exprset_cxc5,
                   x = colnames(exprset_cxc5[j]), 
                   y = colnames(exprset_cxc5[i]),
                   color = colnames(exprset_cxc5[j]),
                   #facet.by = 'source_name_ch1',
                   legend = "none",
                   palette = "npg",
                   add = "point",
                   xlab = 'Treatment'
    )
    print(p + stat_compare_means(aes(label = paste0(..method.., 
                                                    ", p-value = ",
                                                    ..p.format..)),
                                 method = method1, 
                                 size=3,
                                 label.x = .6,
                                 label.y = max(exprset_cxc5[, i], 
                                               na.rm = TRUE))+ #
            stat_compare_means(comparisons = my_comparisons, 
                               method = method2, 
                               size=3,
                               label = "p.format"))
  }
}
