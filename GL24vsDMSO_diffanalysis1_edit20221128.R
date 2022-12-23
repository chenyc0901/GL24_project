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

proteincoding_data %>% filter(abs(logFC_Hs578T)>1 & FDR_Hs578T<0.01) %>% glimpse()

draw_volcanoplot <- function(cell){
                    logFC <- paste0("logFC_",cell)
                    FDR <- paste0("FDR_",cell)
                    output <- paste0(cell,".pdf")
                    Up_regulated <- proteincoding_data[logFC] > 1 & proteincoding_data[FDR] < 1e-2 
                    Down_regulated <- proteincoding_data[logFC] < -1 & proteincoding_data[FDR] < 1e-2 
                    proteincoding_data <- proteincoding_data %>% mutate(color = ifelse(Up_regulated, 'red', ifelse(Down_regulated, 'blue', 'grey')))
                    keyvals.color <- proteincoding_data$color
                    names(keyvals.color)[keyvals.color == 'red'] <- 'Up-regulated'
                    names(keyvals.color)[keyvals.color == 'blue'] <- 'Down-regulated'
                    names(keyvals.color)[keyvals.color == 'grey'] <- 'Not-significant'
                    volcanoplot <- EnhancedVolcano(proteincoding_data,
                                   lab=NA,
                                   title=cell,
                                    xlab = bquote(~Log[2]~ 'fold change'),
                                    ylab = bquote(~-Log[10]~ '(adjusted p value)'),
                                    x=logFC,
                                    y=FDR,
                                    colCustom=keyvals.color,
                                    legendPosition = 'right',
                                    legendLabSize = 10,
                                    FCcutoff = 1,
                                    pCutoff= 1e-2,
                                    colAlpha = 1,
                                    drawConnectors = TRUE,
                                    widthConnectors = 0.75,
                                    maxoverlapsConnectors = Inf,
                                    colConnectors = 'black',
                                    #col=c("blue", "grey", "grey", "red"),
                                    max.overlaps = Inf) +
                      scale_x_continuous(breaks = seq(-10,10,5), limits = c(-10,10))+
                      scale_y_continuous(breaks = seq(0,30,5), limits = c(0,30))
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

#filter condition(MDAMB231,MDAMB157): FDR <0.01, abs(logFC) > 1, n=1047
gene.list <- proteincoding_data %>% 
                   dplyr::filter(FDR_MDAMB231 < 0.01 & 
                                 FDR_MDAMB157 < 0.01 & 
                                 abs(logFC_MDAMB231)> 1 & 
                                 abs(logFC_MDAMB157) > 1 ) %>% 
                   mutate(name=paste(.$gene_id, .$hgnc_symbol, sep=";")) %>%
                   .$name %>% 
                   as.character()

gene.list %>% str()
logCPM_mt <- logCPM_mt[gene.list,]
logCPM_mt%>% head()
logCPM_mt_transpose <- t(scale(t(logCPM_mt)))
col.pan <- colorpanel(100, "blue", "white", "red")

pdf("heatmap_1047.pdf")
ColSideColors <- cbind(Treatment=rep(c("orange", "green"), each=5, times=3),Cell_type=c(rep("brown1",10),rep("mediumpurple2",20)))
heatmap3(logCPM_mt_transpose, 
          col=col.pan, 
          labRow=FALSE,
          Rowv=TRUE,
          Colv=TRUE,
          scale="none",
          ColSideColors = ColSideColors, 
          margin=c(10,9),
          cexRow=0.7,
          cexCol=0.8)
dev.off()



#filter condition(MDAMB231,MDAMB157): FDR <0.01, abs(logFC) > 1,same direction, n=959

gene.list <- proteincoding_data %>% 
             dplyr::filter(FDR_MDAMB231 < 0.01 & 
             FDR_MDAMB157 < 0.01 & 
             abs(logFC_MDAMB231)> 1 & 
             abs(logFC_MDAMB157) > 1 &
             logFC_MDAMB231*logFC_MDAMB157>0) %>% 
             mutate(name=paste(.$gene_id, .$hgnc_symbol, sep=";")) %>%
             .$name %>% 
             as.character()

gene.list %>% str()
logCPM_mt <- logCPM_mt[gene.list,]
logCPM_mt%>% head()
logCPM_mt_transpose <- t(scale(t(logCPM_mt)))
col.pan <- colorpanel(100, "blue", "white", "red")

pdf("heatmap_same_direction_959.pdf")
ColSideColors <- cbind(Treatment=rep(c("orange", "green"), each=5, times=3),Cell_type=c(rep("brown1",10),rep("mediumpurple2",20)))
heatmap3(logCPM_mt_transpose, 
         col=col.pan, 
         labRow=FALSE,
         Rowv=TRUE,
         Colv=TRUE,
         scale="none",
         ColSideColors = ColSideColors, 
         margin=c(10,9),
         cexRow=0.7,
         cexCol=0.8)
dev.off()
