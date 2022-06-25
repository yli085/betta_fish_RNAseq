# script to extract counts for only Ovary/Testis groups, and Long fin/Short fin groups, respectively.
# Perform normalization and wilcoxon test on gene expression of a given set of genes, and
# make a boxplot with significance.

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(writexl)


# input data --------------------------------------------------------------

## all data
data0 <- read_delim("data/counts_Rmatrix.txt", delim="\t", col_types=list(Geneid="c", .default="i"))
meta0 <- read_excel("data/meta/betta_fish_meta.xlsx")


# Ovary vs Testis --------------------------------------------------------------

## Input data: ovary and testis
meta <- meta0[-1]
meta <- meta %>% filter(Group %in% c("Testis", "Ovary"))
meta <- meta %>% column_to_rownames("Sample")
meta$Group <- factor(meta$Group, levels = c("Ovary", "Testis"))
data <- data0[c("Geneid", rownames(meta))]
data <- data %>% column_to_rownames("Geneid")

## lcpm
group <- meta$Group
x <- DGEList(counts=data, group=group)
x <- calcNormFactors(x, method = "TMM")
# x$samples$norm.factors
lcpm_df <- cpm(x, normalized.lib.sizes=TRUE, log = T)

# extract igenes
id_tab <- read_excel("data/ids.xlsx", sheet = 1, range="A2:B10")
id_tab$gene_id <- str_replace(id_tab$`gene_id in BS_25104_geneID.csv`, "model.chr2", "TU.Chr2")
# make df
df <- lcpm_df[match(id_tab$gene_id, row.names(lcpm_df)),] %>% as.data.frame()
rownames(df) <- id_tab$gene_name 
df <- df %>% rownames_to_column("gene")
# transform data
ngenes = dim(id_tab)[1]
df <- df %>% pivot_longer(2:ncol(df), names_to = "sample", values_to = "expr")
df <- df %>% mutate(group=rep(meta$Group, ngenes))
# factorize variables
df$group <- factor(df$group, levels=c("Ovary", "Testis")) 
df$gene <- as_factor(df$gene)


## plot
# p value
pValue <- df %>% group_by(gene) %>%
  do(test = wilcox.test(expr~group, data=., paired=FALSE)) %>%
  summarise(gene, pval = test$p.value)

# multiple comparison adjust
pValue <- p.adjust(pValue$pval, method="fdr")

my_annotations <- as.character(symnum(pValue, #corr = FALSE, na = FALSE,
                                      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***",'**','*','ns')))

exprMax <- df %>% group_by(gene) %>% summarise(maxExpr=max(expr))

anno <- data.frame(x=exprMax$gene,y=exprMax$maxExpr+1,label=my_annotations,gene=exprMax$gene)
anno_sig <- anno[which(anno$label!='ns'),]
anno_ns <-  anno[which(anno$label=='ns'),]

# plot
# line width (box, axis), text/label character size (x, y axis, legend, annotation text)
# customized: y axis breaks, legend labels, legend positions
bxplt <- ggplot(df, aes(x=gene, y=expr, group = interaction(group, gene)))

bxplt+geom_boxplot(aes(fill=group),
                   width=0.5, lwd=1, fatten=1) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  scale_fill_manual(values=c("dodgerblue", "tan2"))+
  scale_y_continuous(breaks=c(0, 4, 8, 12))+
  geom_text(data =anno_sig, aes(x, y, label=label, group=NULL),
            size=14, color="red", fontface="bold") +
  geom_text(data =anno_ns, aes(x, y, label=label, group=NULL),
            size=10, color="black") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = c(0.9, 0.9),
        legend.box.background = element_rect(colour = "black", size=1)) +
  theme(axis.title=element_text(size=24),
        axis.text = element_text(color='black', size=20),
        axis.text.x = element_text(face="italic", size =24, angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.line = element_line(colour = "black", size=1),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

# saved as 12 x 6 inches pdf.



# Long fin vs Short fin ---------------------------------------------------

## Input data: Long-fin & Short-fin
meta <- meta0[-1]
meta <- meta %>% filter(Group %in% c("Short_tail", "Long_tail"))
meta <- meta %>% column_to_rownames("Sample")
meta$Group <- factor(meta$Group, levels = c("Long_tail", "Short_tail"))
data <- data0[c("Geneid", rownames(meta))]
data <- data %>% column_to_rownames("Geneid")

## lcpm 
group <- meta$Group
x <- DGEList(counts=data, group=group)
x <- calcNormFactors(x, method = "TMM")
# x$samples$norm.factors
lcpm_df <- cpm(x, normalized.lib.sizes=TRUE, log = T)

# extract igenes
id_tab <- read_excel("data/ids.xlsx", sheet = 2, range="A2:B10")
id_tab$gene_id <- str_replace(id_tab$`gene_id in BS_25104_geneID.csv`, "model.chr17", "TU.Chr17")
# make df
df <- lcpm_df[match(id_tab$gene_id, row.names(lcpm_df)),] %>% as.data.frame()
rownames(df) <- id_tab$gene_name 
df <- df %>% rownames_to_column("gene")
# transform data
ngenes = dim(id_tab)[1]
df <- df %>% pivot_longer(2:ncol(df), names_to = "sample", values_to = "expr")
df <- df %>% mutate(group=rep(meta$Group, ngenes))
# factorize variables
df$group <- factor(df$group, levels=c("Long_tail", "Short_tail")) 
df$gene <- as_factor(df$gene)

## plot
# p value
pValue <- df %>% group_by(gene) %>%
  do(test = wilcox.test(expr~group, data=., paired=FALSE)) %>%
  summarise(gene, pval = test$p.value)

# multiple comparison adjust
pValue <- p.adjust(pValue$pval, method="fdr")

my_annotations <- as.character(symnum(pValue, #corr = FALSE, na = FALSE,
                                      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***",'**','*','ns')))

exprMax <- df %>% group_by(gene) %>% summarise(maxExpr=max(expr))

anno <- data.frame(x=exprMax$gene,y=exprMax$maxExpr+1,label=my_annotations,gene=exprMax$gene)
anno_sig <- anno[which(anno$label!='ns'),]
anno_ns <-  anno[which(anno$label=='ns'),]

# plot
# line width (box, axis), text/label character size (x, y axis, legend, annotation text)
# customized: y axis breaks, legend labels, legend positions

bxplt <- ggplot(df, aes(x=gene, y=expr, group = interaction(group, gene)))

bxplt+geom_boxplot(aes(fill=group),
                   width=0.5, lwd=1, fatten=1) +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  scale_fill_manual(values=c("dodgerblue", "tan2"),
                    labels = c("Halfmoon", "HMPK"))+
  scale_y_continuous(breaks=c(0, 4, 8))+
  geom_text(data =anno_sig, aes(x, y, label=label, group=NULL),
            size=14, color="red", fontface="bold") +
  geom_text(data =anno_ns, aes(x, y, label=label, group=NULL),
            size=10, color="black") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = c(0.92, 0.93),
        legend.box.background = element_rect(colour = "black", size=1)) +
  theme(axis.title=element_text(size=24),
        axis.text = element_text(color='black', size=20),
        axis.text.x = element_text(face="italic", size =24, angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.line = element_line(colour = "black", size=1),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

# saved in 12 X 6 inches pdf
