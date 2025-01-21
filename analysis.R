library(DESeq2)
library(org.Hs.eg.db)

# load in the data
sample_info <- read.table("sample_info.txt", row.names=1)
counts <- read.table("raw_feature_counts.txt", row.names = 1)

# perform differential gene expression with DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

# filter to remove low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# choose reference level; I want to identify genes that are differentially expressed in the M003.1 condition compared to the mCherry negative control condition 
dds$condition <- relevel(dds$condition, ref = "cherry")

# differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# create symbol column for gene names
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                      keys=ens.str,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")


