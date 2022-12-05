library(tidyverse)
library(edgeR)
library(biomaRt)
library(enrichplot)
library(corrplot) 
library(VennDiagram)
library(gplots)
library(factoextra) # for PCA plot
library(EnhancedVolcano)
library(heatmap3)

rm(list=ls())
set.seed(100)
# readcount data from RSEM (expected counts)
data  <-  read.delim(file="gene_readcount_raw.tsv", header = TRUE, row.names = 1)
groups = factor(rep(c("Hs578T_DMSO","Hs578T_GL24","MDAMB157_DMSO", "MDAMB157_GL24","MDAMB231_DMSO","MDAMB231_GL24"), each = 5))
# To create DGEList object
y <- DGEList(counts = data, group = groups)
#To filter condition use default method filterByExpr
keep  <- filterByExpr(y)
#To change library size after filter
y <-  y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+groups)
y <- estimateDisp(y,design)

#log transformation
#logCPM = log2(CPM +2/L):L is the average library size in millions
lcpm <- cpm(y,log=TRUE)
#normalization plot
boxplot(lcpm,las=2, col=ncol(y), main="normalization result", cex.axis=0.6)
#correlation
corrplot.mixed(cor(lcpm), tl.col = "black", upper = "pie", tl.pos = "lt", number.cex=0.4)

#PCA plot
pca <- prcomp(t(lcpm),scale=TRUE)
pdf("PCAplot1.pdf")
grp <- factor(y$samples$group)
eig.val <- get_eigenvalue(pca)
head(eig.val, n = 10L)
fviz_pca_ind(pca,
             geom = c("point"),
             habillage = grp,
             pointsize = 3,
             #addEllipses = TRUE,
             #ellipse.level=0.95,
             mean.point = FALSE,
             repel = TRUE     # Avoid text overlapping
)+scale_color_brewer(palette="Set1") +
  theme_minimal()
dev.off()

# fit model
fit <- glmQLFit(y,design)

# GL21 vs DMSO expression in Hs578T cell line
qlf_Hs578T <- glmQLFTest(fit, contrast=c(-1,1,0,0,0,0))
dgetable_Hs578T<- topTags(qlf_Hs578T,n = nrow(qlf_Hs578T))  %>% 
                  .$table  %>% 
                   rownames_to_column('gene_id')

# GL21 vs DMSO expression in MDAMB157 cell line 
qlf_MDAMB157 <- glmQLFTest(fit, contrast=c(0,0,-1,1,0,0))
dgetable_MDAMB157 <- topTags(qlf_MDAMB157,n=nrow(qlf_MDAMB157)) %>% 
                    .$table %>% 
                    rownames_to_column('gene_id')

# GL21 vs DMSO expression in MDAMB231 cell line
qlf_MDAMB231 <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,1))
dgetable_MDAMB231 <- topTags(qlf_MDAMB231,n=nrow(qlf_MDAMB231)) %>% 
                    .$table %>% 
                    rownames_to_column('gene_id')

#convert transcript id to gene symbol
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",version = 104)
geneid <- dgetable_Hs578T$gene_id
G_list <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","gene_biotype", "hgnc_symbol","entrezgene_id"),
                values = geneid,
                mart = mart)

# To transform readcounts data for merge
data_readcounts  <- data %>% rownames_to_column("gene_id")
# To merge all datafile
mergedata <- left_join(dgetable_Hs578T, G_list, by = c("gene_id"="ensembl_gene_id")) %>% 
             left_join(., dgetable_MDAMB157, by = "gene_id", suffix = c("","_MDAMB157"))  %>% 
             left_join(., dgetable_MDAMB231, by = "gene_id", suffix = c("_Hs578T","_MDAMB231")) %>% 
             left_join(., data_readcounts,by = "gene_id")


##### To analysize the result ######
#To draw volcanoplot 
proteincoding_data <- mergedata %>% filter(mergedata$gene_biotype=="protein_coding")
#overlap(112 DE gene overlap, up=42, down=70 )
overlap_up <- proteincoding_data$logFC_MDAMB157 > 2 & 
  proteincoding_data$logFC_MDAMB231 > 2 & 
  proteincoding_data$FDR_MDAMB157 < 10e-5 & 
  proteincoding_data$FDR_MDAMB231 < 10e-5 
overlap_down <-  proteincoding_data$logFC_MDAMB157 < -2 & 
  proteincoding_data$logFC_MDAMB231 < -2 & 
  proteincoding_data$FDR_MDAMB157 < 10e-5 & 
  proteincoding_data$FDR_MDAMB231 < 10e-5 
#add color panel
proteincoding_data <- proteincoding_data %>% mutate(color=ifelse(overlap_up,'red',ifelse(overlap_down,'blue','grey')))

draw_volcanoplot <- function(cell){
                    logFC <- paste0("logFC_",cell)
                    FDR <- paste0("FDR_",cell)
                    output <- paste0(cell,".pdf")
                    keyvals.color <- proteincoding_data$color
                    names(keyvals.color)[keyvals.color == 'red'] <- 'overlap_up'
                    names(keyvals.color)[keyvals.color == 'blue'] <- 'overlap_down'
                    names(keyvals.color)[keyvals.color == 'grey'] <- 'non-overlap'
                    volcanoplot <- EnhancedVolcano(proteincoding_data,
                                   lab=NA,
                                    #lab=proteincoding_data$hgnc_symbol,
                                    #selectLab=proteincoding_data[which(filtercondition),]$hgnc_symbol,
                                    xlab = bquote(~Log[2]~ 'fold change'),
                                    x=logFC,
                                    y=FDR,
                                    colCustom=keyvals.color,
                                    legendPosition = 'top',
                                    legendLabSize = 10,
                                    FCcutoff = 2,
                                    pCutoff= 10e-5,
                                    colAlpha = 0.3,
                                    drawConnectors = TRUE,
                                    widthConnectors = 0.75,
                                    maxoverlapsConnectors = Inf,
                                    colConnectors = 'black',
                                    max.overlaps = Inf) 
                   ggsave(volcanoplot, file=output, width=17 ,height=17, units="cm", device='pdf')
}
draw_volcanoplot("Hs578T")
draw_volcanoplot("MDAMB231")
draw_volcanoplot("MDAMB157")

#To draw heatmap by logCPM
logCPM  <- cpm(y, log=TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 
logCPM <- left_join(logCPM, G_list, by = c("gene_id"="ensembl_gene_id")) %>% 
  filter(gene_biotype=="protein_coding") %>% 
  distinct(gene_id,.keep_all=TRUE)
rownames(logCPM) <- logCPM %>%
  mutate(name = paste(.$gene_id, .$hgnc_symbol, sep=";")) %>% 
  .$name

logCPM_mt <- logCPM %>%
  dplyr::select(-c(gene_id,entrezgene_id,gene_biotype,hgnc_symbol)) %>% 
  as.matrix()

#filter condition(MDAMB231,MDAMB157): FDR <10e-5, logFC > 3
gene.list <- proteincoding_data %>% 
                   dplyr::filter(FDR_MDAMB231 < 10e-5 & 
                                 FDR_MDAMB157 < 10e-5 & 
                                 abs(logFC_MDAMB231)> 3 & 
                                 abs(logFC_MDAMB157) > 3 ) %>% 
                   mutate(name=paste(.$gene_id, .$hgnc_symbol, sep=";")) %>%
                   .$name %>% 
                   as.character()

gene.list %>% str()
logCPM_mt <- logCPM_mt[gene.list,]
logCPM_mt%>% head()
logCPM_mt_transpose <- t(scale(t(logCPM_mt)))
col.pan <- colorpanel(100, "blue", "white", "red")

pdf("heatmap.pdf")
ColSideColors <- cbind(Treatment=rep(c("orange", "green"), each=5, times=3),Cell_type=c(rep("brown1",10),rep("mediumpurple2",20)))
heatmap3(logCPM_mt_transpose, 
          col=col.pan, 
          Rowv=TRUE,
          Colv=TRUE,
          scale="none",
          ColSideColors = ColSideColors, 
          margin=c(10,9),
          cexRow=0.7,
          cexCol=0.8)
dev.off()
