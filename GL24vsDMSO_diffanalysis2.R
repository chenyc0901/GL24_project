library(tidyverse)
library(edgeR)
library(biomaRt)
library(gplots)
library(org.Hs.eg.db)
library(clusterProfiler)
library(statmod)
library(EnhancedVolcano)

rm(list=ls())
set.seed(100)
# readcount data from RSEM (expected counts)
data  <-  read.delim(file="gene_readcount_157_231.txt", header = TRUE, row.names = 1)
groups <-  factor(rep(c("MDAMB157_DMSO", "MDAMB157_GL24", "MDAMB231_DMSO", "MDAMB231_GL24"), each = 5))
cells <- factor(rep(c("MDAMB157", "MDAMB231"), each=10))
treat <- factor(rep(c("DMSO", "GL24"), each=5, times=2), levels=c("DMSO","GL24"))

# To create DGEList object
y <- DGEList(counts = data,group=groups)
#To filter condition use default method filterByExpr
keep  <- filterByExpr(y)
#To change library size after filter
y <-  y[keep, ,keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~cells+treat)
rownames(design) <- colnames(y)
y <- estimateDisp(y,design, robust=TRUE)


fit <- glmQLFit(y,design)
lrt <- glmQLFTest(fit)

Treatment_compare<- topTags(lrt,n=nrow(lrt)) %>%
  .$table %>% 
  rownames_to_column("gene_id")

#convert transcript id to gene symbol
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",version = 104)
geneid <- Treatment_compare$gene_id
G_list <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","gene_biotype", "hgnc_symbol","entrezgene_id"),
                values = geneid,
                mart = mart)

#transform read counts data for merge
data_readcounts  <- data %>% rownames_to_column("gene_id")
# merge all data files
mergedata <- left_join(Treatment_compare, G_list, by = c("gene_id"="ensembl_gene_id")) %>% 
             left_join(., data_readcounts,by = "gene_id")
#Output file: mergedata_difftreatment.csv
write.csv(mergedata, file = "mergedata_difftreatment.csv", row.names = FALSE)
#select protein coding gene for following process
proteincoding_data <- mergedata %>% 
                      filter(mergedata$gene_biotype=="protein_coding")
proteincoding_data %>% dim() #14059 29

##### To analyze result #######
#To draw volcano plot
#filter condition: logFC >3, and FDR<10e-7
filtercondition <- abs(proteincoding_data$logFC)>4 & proteincoding_data$FDR<10e-7

pdf("difftreat_normalizedCell_volcano.pdf",width=12, height=12)
EnhancedVolcano(proteincoding_data,
                lab=proteincoding_data$hgnc_symbol,
                selectLab=proteincoding_data[which(filtercondition),]$hgnc_symbol,
                xlab = bquote(~Log[2]~ 'fold change'),
                x="logFC",
                y="FDR",
                legendPosition = 'top',
                legendLabSize = 10,
                FCcutoff = 2,
                pCutoff= 10e-5,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colConnectors = 'black',
                max.overlaps = Inf)
dev.off()



#####GSEA analysis#####
#load hallmark gene set
file_gmt_entrez =  'h.all.v2022.1.Hs.entrez.gmt'
gs_hallmark_entrez <-  read.gmt(file_gmt_entrez)
geneList_diffTreatment <-  proteincoding_data %>%
                           distinct(entrezgene_id, .keep_all=TRUE) %>% 
                           arrange(desc(logFC)) %>%
                           na.omit()
geneList <- geneList_diffTreatment$logFC %>% 
            setNames(geneList_diffTreatment$entrezgene_id)

gsea_hallmark <-  GSEA(geneList, TERM2GENE = gs_hallmark_entrez)
gsea_hallmark_table <- gsea_hallmark@result %>% 
                       dplyr::as_tibble() %>% 
                       arrange(desc(NES))
gsea_hallmark_table %>% glimpse()


###GO analysis###
geneset <- mergedata %>% 
  dplyr::filter(FDR < 0.01,abs(logFC)>1, gene_biotype == "protein_coding" & hgnc_symbol!="")
# dplyr::filter(FDR_MDAMB231 < 0.01 & FDR_MDAMB157 < 0.01 & abs(logFC_MDAMB231)> 2 & abs(logFC_MDAMB157) > 2 & gene_biotype == "protein_coding" & hgnc_symbol!="")
geneset %>% dim() # 2342 genes for GO analysis
ora_GO <-  enrichGO(
  gene = geneset$gene_id,
  OrgDb = org.Hs.eg.db,
  keyType= "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
)
# GO result :ora_GO@result

### KEGG analysis###
library(R.utils)# This package is used for download kegg
R.utils::setOption("clusterProfiler.download.method","auto")
# To get ENTRENSID id for enrichKEGG
ENSEMBL  <- geneset$gene_id
ENSEMBL  %>% str() #2342
egIDs <- stack(mget(ENSEMBL, org.Hs.egENSEMBL2EG, ifnotfound = NA))

kegg_GO  <- enrichKEGG(
  gene=egIDs$values,
  organism = 'hsa',
  pAdjustMethod = "BH",
)
#KEGG_result: kegg_GO@result
