library(EnhancedVolcano)
library(ggrepel)

# select genes to show labels for
labels = c("IL6","TNF", "IL10", "IL1B", "PIM1", "PIM2", "NFKB2", "CD274", "NFKBIZ", "CD80", "ACOD1", "CCL8", "CCL4", "TNFSF18", "IFNB1", "TNFAIP3" , "RELB", "ZC3H12A", "CCL2", "MMP9", "CCL4L2", "CXCL13", "NFKBIA", "SOCS3", "EBI3", "STARD13", "FAXDC2", "BMP6", "PDK4", "PIK3IP1", "BCL2A1", "TNFRSF4", "CXCL10")

# make volcano plot
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 10e-10,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                selectLab = labels,
                gridlines.minor = FALSE,
                gridlines.major = FALSE)
