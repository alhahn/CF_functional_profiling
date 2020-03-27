library("ggplot2")
library("phyloseq")
library("DESeq2")
library("vegan")

##Loading OTU and taxa tables
my_data.frame = read.csv(file = "Path_otutable.csv", header = T, sep = ",")
rownames(my_data.frame) <- paste0("OTU", 1:nrow(my_data.frame))
my_otu_table = otu_table(my_data.frame, taxa_are_rows = TRUE)

my_tax_data.frame = read.csv(file = "Path_taxtable.csv", header=F, sep=",")
rownames(my_tax_data.frame) <- paste0("OTU", 1:nrow(my_tax_data.frame))
my_tax_tab = tax_table(my_tax_data.frame, errorIfNULL = TRUE)
rownames(my_tax_tab) <- paste0("OTU", 1:nrow(my_tax_tab))
colnames(my_tax_tab) <- c("Superclass1", "Superclass2", "Superclass3", "Superclass4", "Pathway")

##Create phyloseq class experiment-level object
physeq = phyloseq(my_otu_table, my_tax_tab)

## Loading Sample data and labeling which samples by Clinical state and Pt ID
sampledata = read.csv("meta_Landen_Path.csv", header = F, sep = ",")
sampledata <- as.data.frame(sampledata)
rownames(sampledata) <- sampledata[,1]
sampledata[,1] <- NULL
colnames(sampledata) <- c("ETF", "PtID")
sampledata = sample_data(sampledata)
#```

######DESEQ2###############

##Merge phyloseq and sampledata to create new phylo-class experiment-level object with taxa, samples, and variables
physeq2 = merge_phyloseq(sampledata, physeq)

# Changing phyloseq object into DEseq object
patho.des = phyloseq_to_deseq2(physeq2, ~ ETF)

# calculate geometric means prior to estimate size factors - this is because we have zeros
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(patho.des), 1, gm_mean)

#Size factor and dispersion estimates
patho.des = estimateSizeFactors(patho.des, geoMeans = geoMeans)
patho.des <- estimateDispersions(patho.des)
patho.des.var <- getVarianceStabilizedData(patho.des)
patho.des.deseq = DESeq(patho.des, fitType="local")


# Filering to find the top most differentially abundant OTUs, ordering by p value, and removing taxa with NA value
#contrast - "variable", "numerator", "denominator", so positive log fold change would be higher in numerator 
#and negative log fold change would be higher in the denominator
res = results(patho.des.deseq, cooksCutoff = FALSE, contrast = c("ETF", "T", "F"))
res = res[order(res$padj, na.last=NA), ]
res 
#this is where you can see that the log2 fold change is being calculate between Clinical State
#log2 fold change (MAP): Clinical State
#Wald test p-value: Clinical State
#export res table into excel
write.csv(x=res, file="Pathway_res.csv")


#set significance
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq2)[rownames(sigtab), ], "matrix"))

head(sigtab)

dim(sigtab)

#export significant log2 fold change into excel
write.csv(x=sigtab, file="Pathway_sigtab.csv")

### Looking at the Pathways that are significantly different

x = tapply(sigtab$log2FoldChange, sigtab$Pathway, function(x) max(x))
x = sort(x, TRUE)
sigtab$Pathway = factor(as.character(sigtab$Pathway), levels=names(x))

# Colored by Pathway
patho.sp.genus.l2fc_plot <- ggplot(sigtab, aes(x=log2FoldChange, y=Pathway, color=Superclass2)) + geom_vline(xintercept = 0.0, color = "gray", size = 1) + geom_point(size=3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
patho.sp.genus.l2fc_plot
ggsave("patho.sp.genus.l2fc.pdf", patho.sp.genus.l2fc_plot)
