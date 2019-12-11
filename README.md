# my_TCGA
TCGA data analysis with R
#Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.1 (2019-07-05)


# Load the bioconductor installer. 
# Install the main RTCGA package
# https://bioconductor.org/install

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

library("BiocManager")


# Install the clinical and mRNA gene expression data packages

BiocManager::install("RTCGA")   #RTCGA_1.16.0
BiocManager::install("RTCGA.clinical")  #RTCGA.clinical_20151101.16.0
BiocManager::install('RTCGA.rnaseq')  #RTCGA.rnaseq_20151101.16.0
BiocManager::install("RTCGA.mRNA")   #RTCGA.mRNA_1.14.0
BiocManager::install("RTCGA.mutations")   
BiocManager::install("RTCGA.miRNASeq")
BiocManager::install("RTCGA.methylation")  
BiocManager::install("RTCGA.PANCAN12")
BiocManager::install("RTCGA.CNV")
BiocManager::install("RTCGA.RPPA")

library("RTCGA")
library("RTCGA.clinical") 
library('RTCGA.rnaseq')
library("RTCGA.mRNA")
library('RTCGA.mutations') 
library("RTCGA.miRNASeq")
library("RTCGA.methylation")
library("RTCGA.PANCAN12")
library("RTCGA.CNV")
library("RTCGA.RPPA")



BiocManager::install("RTCGA.clinical.20160128")
BiocManager::install("RTCGA.mutations.20160128")
BiocManager::install("RTCGA.rnaseq.20160128")
BiocManager::install("RTCGA.CNV.20160128")
BiocManager::install("RTCGA.RPPA.20160128")
BiocManager::install("RTCGA.mRNA.20160128")
BiocManager::install("RTCGA.miRNASeq.20160128")
BiocManager::install("RTCGA.methylation.20160128")

library("RTCGA.clinical.20160128")
library("RTCGA.mutations.20160128")
library("RTCGA.rnaseq.20160128")
library("RTCGA.CNV.20160128")
library("RTCGA.RPPA.20160128")
library("RTCGA.mRNA.20160128")
library("RTCGA.miRNASeq.20160128")
library("RTCGA.methylation.20160128")


#Refer from: https://cloud.tencent.com/developer/article/1168376

#获取TCGA-cancer数据
library(RTCGA)
all_TCGA_cancers=infoTCGA()
DT::datatable(all_TCGA_cancers)
colnames(all_TCGA_cancers)
rownames(all_TCGA_cancers)

save(all_TCGA_cancers,file="all_TCGA_cancers.results20191211.Rdata")



#指定任意基因从任意癌症里面获取【芯片表达】数据
#这里我们拿下面3种癌症做示范：
#Breast invasive carcinoma (BRCA)
#Ovarian serous cystadenocarcinoma (OV)
#Lung squamous cell carcinoma (LUSC)

rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA.mRNA)
expr <- expressionsTCGA(BRCA.mRNA, OV.mRNA, LUSC.mRNA,
                        extract.cols = c("GATA3", "PTEN", "XBP1","ESR1", "MUC1"))


expr


nb_samples <- table(expr$dataset)
nb_samples


#其中要注意的是mRNA并不是rnaseq，两者不太一样，具体样本数量，可以看最前面的表格。
#下面简化一下标识，方便可视化展现
expr$dataset <- gsub(pattern = ".mRNA", replacement = "",  expr$dataset)
expr$bcr_patient_barcode <- paste0(expr$dataset, c(1:590, 1:561, 1:154))
expr


#绘制指定基因在不同癌症的表达量区别boxplot
library(ggpubr)


## GATA3
ggboxplot(expr, x = "dataset", y = "GATA3",
          title = "GATA3", ylab = "Expression",
          color = "dataset", palette = "jco")


## PTEN
ggboxplot(expr, x = "dataset", y = "PTEN",
          title = "PTEN", ylab = "Expression",
          color = "dataset", palette = "jco")

### 注意这个配色可以自选的： RColorBrewer::display.brewer.all()  
#这里选择的是 ggsci 包的配色方案，包括： “npg”, “aaas”, “lancet”, “jco”, “ucscgb”, “uchicago”, 
#“simpsons” and “rickandmorty”，针对常见的SCI杂志的需求开发的。

my_comparisons <- list(c("BRCA", "OV"), c("OV", "LUSC"))
ggboxplot(expr, x = "dataset", y = "GATA3",
          title = "GATA3", ylab = "Expression",
          color = "dataset", palette = "jco")+
  stat_compare_means(comparisons = my_comparisons)

#这些统计学检验，也是被包装成了函数：
compare_means(c(GATA3, PTEN, XBP1) ~ dataset, data = expr)


#更多boxplot参数
label.select.criteria <- list(criteria = "`y` > 3.9 & `x` %in% c('BRCA', 'OV')")
ggboxplot(expr, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          combine = TRUE,
          color = "dataset", palette = "jco",
          ylab = "Expression", 
          label = "bcr_patient_barcode",              # column containing point labels
          label.select = label.select.criteria,       # Select some labels to display
          font.label = list(size = 9, face = "italic"), # label font
          repel = TRUE                                # Avoid label text overplotting
)

#其中 combine = TRUE 会把多个boxplot并排画在一起，其实没有ggplot自带的分面好用。
#还可以使用 merge = TRUE or merge = “asis” or merge = "flip" 来把多个boxplot 合并，效果不一样。


#还有翻转，如下：
ggboxplot(expr, x = "dataset", y = "GATA3",
          title = "GATA3", ylab = "Expression",
          color = "dataset", palette = "jco",
          rotate = TRUE)

#更多可视化详见： http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/



#指定任意基因从任意癌症里面获取【测序表达】数据
#还是同样的3种癌症和5个基因做示范，这个时候的基因ID稍微有点麻烦，不仅仅是要symbol还要entrez的ID，
#具体需要看 https://wiki.nci.nih.gov/display/TCGA/RNASeq+Version+2 的解释

#如下：
rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA)
library(RTCGA.rnaseq)
expr <- expressionsTCGA(BRCA.rnaseq, OV.rnaseq, LUSC.rnaseq,
                        extract.cols = c("GATA3|2625", "PTEN|5728", "XBP1|7494","ESR1|2099", "MUC1|4582"))
expr

nb_samples <- table(expr$dataset)
nb_samples

library(ggpubr)
# ESR1|2099
ggboxplot(expr, x = "dataset", y = "`PTEN|5728`",
          title = "ESR1|2099", ylab = "Expression",
          color = "dataset", palette = "jco")


#更多可视化见：http://rtcga.github.io/RTCGA/articles/Visualizations.html




#用全部的rnaseq的表达数据来做主成分分析
## RNASeq expressions
rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA.rnaseq)
library(dplyr) 

## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union


expressionsTCGA(BRCA.rnaseq, OV.rnaseq, HNSC.rnaseq) %>%
  dplyr::rename(cohort = dataset) %>%  
  filter(substr(bcr_patient_barcode, 14, 15) == "01") -> BRCA.OV.HNSC.rnaseq.cancer
pcaTCGA(BRCA.OV.HNSC.rnaseq.cancer, "cohort") -> pca_plot
plot(pca_plot)

#因为是全部的表达数据，所以非常耗时，但是可以很明显看到乳腺癌和卵巢癌关系要近一点，头颈癌症就要远一点。

#用5个基因在3个癌症的表达量做主成分分析

expr %>%  
  filter(substr(bcr_patient_barcode, 14, 15) == "01") -> rnaseq.5genes.3cancers
DT::datatable(rnaseq.5genes.3cancers)
#pcaTCGA(rnaseq.5genes.3cancers, "dataset") -> pca_plot
#plot(pca_plot)

#该包里面的pcaTCGA函数不好用，其实可以自己做PCA分析。

#用突变数据做生存分析
rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA.mutations)
# library(dplyr) if did not load at start
library(survminer)
mutationsTCGA(BRCA.mutations, OV.mutations) %>%
  filter(Hugo_Symbol == 'TP53') %>%
  filter(substr(bcr_patient_barcode, 14, 15) ==
           "01") %>% # cancer tissue
  mutate(bcr_patient_barcode =
           substr(bcr_patient_barcode, 1, 12)) ->
  BRCA_OV.mutations
library(RTCGA.clinical)
survivalTCGA(
  BRCA.clinical,
  OV.clinical,
  extract.cols = "admin.disease_code"
) %>%
  dplyr::rename(disease = admin.disease_code) ->
  BRCA_OV.clinical
BRCA_OV.clinical %>%
  left_join(
    BRCA_OV.mutations,
    by = "bcr_patient_barcode"
  ) %>%
  mutate(TP53 =
           ifelse(!is.na(Variant_Classification), "Mut","WILDorNOINFO")) ->
  BRCA_OV.clinical_mutations
BRCA_OV.clinical_mutations %>%
  select(times, patient.vital_status, disease, TP53) -> BRCA_OV.2plot
kmTCGA(
  BRCA_OV.2plot,
  explanatory.names = c("TP53", "disease"),
  break.time.by = 400,
  xlim = c(0,2000),
  pval = TRUE) -> km_plot

## Scale for 'colour' is already present. Adding another scale for
## 'colour', which will replace the existing scale.

## Scale for 'fill' is already present. Adding another scale for 'fill',
## which will replace the existing scale.

print(km_plot)



#多个基因在多种癌症的表达量热图
rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA.rnaseq)
# perfrom plot
# library(dplyr) if did not load at start
expressionsTCGA(
  ACC.rnaseq,
  BLCA.rnaseq,
  BRCA.rnaseq,
  OV.rnaseq,
  extract.cols = 
    c("MET|4233",
      "ZNF500|26048",
      "ZNF501|115560")
) %>%
  dplyr::rename(cohort = dataset,
                MET = `MET|4233`) %>%
  #cancer samples
  filter(substr(bcr_patient_barcode, 14, 15) ==
           "01") %>%
  mutate(MET = cut(MET,
                   round(quantile(MET, probs = seq(0,1,0.25)), -2),
                   include.lowest = TRUE,
                   dig.lab = 5)) -> ACC_BLCA_BRCA_OV.rnaseq
ACC_BLCA_BRCA_OV.rnaseq %>%
  select(-bcr_patient_barcode) %>%
  group_by(cohort, MET) %>%
  summarise_each(funs(median)) %>%
  mutate(ZNF500 = round(`ZNF500|26048`),
         ZNF501 = round(`ZNF501|115560`)) ->
  ACC_BLCA_BRCA_OV.rnaseq.medians

## `summarise_each()` is deprecated.
## Use `summarise_all()`, `summarise_at()` or `summarise_if()` instead.
## To map `funs` over all variables, use `summarise_all()`

heatmapTCGA(ACC_BLCA_BRCA_OV.rnaseq.medians,
            "cohort", "MET", "ZNF500",
            title = "Heatmap of ZNF500 expression")


#一个R包不仅仅是提供一个数据下载接口，更重要的是里面封装了一些便于使用的统计分析函数。
#更多参考 生信技能树GATK4系列教程 ：GATK4的gvcf流程


#http://rtcga.github.io/RTCGA/articles/Visualizations.html
#Facet examples with mutations datasets

## facet example
rm(list=ls())
load("all_TCGA_cancers.results20191211.Rdata")

library(RTCGA.mutations)
# library(dplyr) if did not load at start
mutationsTCGA(
  BRCA.mutations,
  OV.mutations,
  ACC.mutations,
  BLCA.mutations
) %>%
  filter(Hugo_Symbol == 'TP53') %>%
  filter(substr(bcr_patient_barcode, 14, 15) ==
           "01") %>% # cancer tissue
  mutate(bcr_patient_barcode =
           substr(bcr_patient_barcode, 1, 12)) ->
  ACC_BLCA_BRCA_OV.mutations

mutationsTCGA(
  BRCA.mutations,
  OV.mutations,
  ACC.mutations,
  BLCA.mutations
) -> ACC_BLCA_BRCA_OV.mutations_all

ACC_BLCA_BRCA_OV.rnaseq %>%
  mutate(bcr_patient_barcode =
           substr(bcr_patient_barcode, 1, 15)) %>%
  filter(bcr_patient_barcode %in%
           substr(ACC_BLCA_BRCA_OV.mutations_all$bcr_patient_barcode, 1, 15)) %>% 
  # took patients for which we had any mutation information
  # so avoided patients without any information about mutations
  mutate(bcr_patient_barcode =
           substr(bcr_patient_barcode, 1, 12)) %>%
  # strin_length(ACC_BLCA_BRCA_OV.mutations$bcr_patient_barcode) == 12
  left_join(ACC_BLCA_BRCA_OV.mutations,
            by = "bcr_patient_barcode") %>% #joined only with tumor patients
  mutate(TP53 = 
           ifelse(!is.na(Variant_Classification), "Mut", "WILD")
  ) %>%
  select(-bcr_patient_barcode, -Variant_Classification,
         -dataset, -Hugo_Symbol) %>% 
  group_by(cohort, MET, TP53) %>% 
  summarise_each(funs(median)) %>% 
  mutate(ZNF501 = round(`ZNF501|115560`)) ->
  ACC_BLCA_BRCA_OV.rnaseq_TP53mutations_ZNF501medians

heatmapTCGA(
  ACC_BLCA_BRCA_OV.rnaseq_TP53mutations_ZNF501medians,
  "cohort",
  "MET",
  fill = "ZNF501",
  facet.names = "TP53",
  title = "Heatmap of ZNF501 expression"
)


heatmapTCGA(
  ACC_BLCA_BRCA_OV.rnaseq_TP53mutations_ZNF501medians,
  "TP53",
  "MET",
  fill = "ZNF501",
  facet.names = "cohort",
  title = "Heatmap of ZNF501 expression"
)

heatmapTCGA(
  ACC_BLCA_BRCA_OV.rnaseq_TP53mutations_ZNF501medians,
  "TP53",
  "cohort",
  fill = "ZNF501",
  facet.names = "MET",
  title = "Heatmap of ZNF501 expression"
)

