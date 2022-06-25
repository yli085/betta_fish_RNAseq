# Betta fish RNA-seq project

## Questions
-	Study sex determination by comparing a given set of genes’ expression levels of the Ovary tissue and Testis samples.
-	Study tail length determination by comparing a given set of genes’ expression levels of the long fin and short fin samples.
-	PCA analysis of the RNA-seq samples.

## RNA-seq sample information and library construction
For the RNA-Seq analysis, two-month-old fish individuals were used. Four independent RNA libraries were prepared and sequenced: 
1. library of pooled RNA samples of 11 tissues (brain, liver, muscle, eye, skin, scale, fin, intestine, testis, ovary, embryo) of one HMPK individual was constructed and subjected to both full-length transcriptome sequencing using the PacBio Sequel platform and short reads sequencing with the Illumina HiSeq 2500 platform. 
2. Libraries of RNA samples of the caudal fin tissues from HMPK (n = 5), Halfmoon (n = 5), Crowntail (n = 1), Veiltail (n = 1), wild B. splendens (n = 1), B. imbellis (n = 1), B. mahachaiensis (n = 1), B. smaragdina (n = 1) individuals were prepared and sequenced with the Illumina HiSeq 2500 platform. 
3. Libraries of RNA samples of the pectoral fin from dumbo (n = 5) and non-dumbo (n = 5) individuals were prepared and sequenced with the Illumina HiSeq 2500 platform. 
4. Libraries of ovary RNA samples from female (n = 5) and testis samples from male (n = 5) individuals were prepared and sequenced with the Illumina HiSeq 2500 platform. 

## Data analysis
The raw sequencing data were first processed by removing adapters and low-quality reads. The resulting clean paired-end reads were mapped to our assembly using STAR. Then the read counts were summarized at the gene level with our gene annotation file using featureCounts (v1.5.0-p3). Reads counts were further normalized with the trimmed mean of M-values (TMM) method implemented in the edgeR package (3.18.1). The log CPM (Counts Per Million) values were obtained using the cpm function from edgeR.

**Cite**: Zhang, W., Wang, H., Brandt D.Y.C., Hu B., Sheng J., Wang M., Luo H., **Li Y**., et. al. (2022) The genetic architecture of phenotypic diversity in the betta fish (Betta splendens). *Science Advances* (Accepted)  
