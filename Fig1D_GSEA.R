library(clusterProfiler)
library(msigdbr)
library(tidyverse)

# Download Hallmark gene sets from MSigDB
hs_gsea_H <- msigdbr(species = "Homo sapiens", 
                     category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Pull out just the columns corresponding to gene symbols and LogFC 
mydata.df.sub <- data.frame(geneID = res$symbol, LogFC = res$log2FoldChange)

# Construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res.H <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_H, verbose=FALSE)
myGSEA.df.H <- as_tibble(myGSEA.res.H@result)

# select the top 10 enriched pathways
myGSEA.df.H.top10 <- myGSEA.df.H %>% 
  top_n(10, NES)

# plot
ggplot(myGSEA.df.H.top10, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(y = "Normalized enrichment score", x="") +
  coord_flip() + # Flip the coordinates to make it horizontal
  guides(fill = guide_colorbar(reverse = TRUE)) +
  theme_bw()