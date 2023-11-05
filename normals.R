options(stringsAsFactors=F)
Sys.setenv("VROOM_CONNECTION_SIZE"=500072)
library(dplyr)
library(GEOquery)
library(openxlsx)
library(BioNERO)
library(tidyverse)
library(janitor)
library(DESeq2)
library(FCBF)
library(gtools)
library(enrichR)
library(annotables)
library(WGCNA)
library(officer)
library(flextable)

setwd("C:/Users/Manisri/OneDrive - International Institute of Information Technology/mouse_models_hcc/human_transplantation_datasets/")

source("modified_functions.R")

g1 = getGEO(filename="GSE148355_series_matrix.txt.gz", AnnotGPL=F, GSEMatrix=T)

xx = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,-2]

naming = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,c(1,2)]

gene_mapping = naming[,2] %>% setNames(naming[,1])

coding_features = gene_mapping %>% data.frame() %>% setNames("gene_symbol") %>% 
  subset(gene_symbol %in% (grch38 %>% filter(biotype == "protein_coding") %>% pull(symbol) %>% unique   )      ) %>% 
  rownames %>% sort

req_samples = xx %>% names %>% grep("TG|HPC|DH013|N003",.,invert = T,value = T) %>% .[-1] 

fpkm_data = xx %>% column_to_rownames("ID") %>% .[coding_features,req_samples]

req_genes = which(rowSums(fpkm_data == 0) < 54) %>% names

req_exp = log2(fpkm_data[req_genes,] + 1)

all( rownames(fpkm_data)[(fpkm_data %>% rowSums() ) == 0]  == rownames(req_exp)[(req_exp %>% rowSums() ) == 0]  )

top_varying_genes = req_exp %>% 
  as.matrix %>% rowVars() %>% order(decreasing = T) %>% .[1:3000]

datExpr0 = req_exp[top_varying_genes,] %>% t %>% data.frame()

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

hh=req_exp[top_varying_genes,]

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#tiff("coding_sample_b4_filter.tiff", compression = "lzw",width = 7, height = 6, units = "in", res=600)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = 50, col = "red");
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

tiff("coding_sample_after_filter.tiff", compression = "lzw",width = 7, height = 6, units = "in", res=600)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

hh = t(datExpr)

sft <- SFT_fit_1(exp=hh, net_type="signed", cor_method="pearson",rsquared = 0.8)

softpower <- 14

sft_thrsh=sft$dtf %>% mutate(special = ifelse(power == softpower,"yes","no") %>% factor())
mm = round(sft_thrsh$meank[sft_thrsh$power == softpower],1 )
ggpubr::ggarrange(sft$plot1,sft$plot2 + geom_hline(yintercept = mm, color="brown3") + 
                    geom_text(aes(0,mm, label=mm),hjust=0.1,vjust=-0.2,color="blue"))

p1=ggpubr::ggscatter(sft_thrsh, x = "power", y = "fit", ylim = c(0,1), label = "power", size = 1, font.label = 10, 
                     color = "special", font.tickslab = 10, ytickslab.rt = 90,
                     palette = c(yes = "blue",no="black")) + 
  geom_hline(yintercept = 0.8, color = "brown3") + theme_bw() + ggtitle("Scale independence") +
  xlab("Soft threshold (power)") + ylab(expression(paste("Scale-free topology fit - ",R^{2} ) ) ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        panel.grid = element_blank())

p2=ggpubr::ggscatter(sft_thrsh, x = "power", y = "meank",
                     xlab = "Soft threshold (power)", ylab = "Mean connectivity (k)", 
                     title = "Mean connectivity", label = "power", 
                     size = 1, font.label = 10, color = "special", font.tickslab = 10, 
                     ytickslab.rt = 90, palette = c(yes = "blue",no="black")) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", panel.grid = element_blank())+
  geom_hline(yintercept = mm, color="brown3") +
  geom_text(aes(0,mm, label=mm),hjust=0.1,vjust=-0.2,color="blue")

scale_free_plot = ggpubr::ggarrange(p1, p2) 

ggpubr::ggarrange(p1, p2)

ggsave("coding_scale_free_thph.tiff", height=4, width = 7, device = "tiff", dpi = 600, plot = scale_free_plot)

#net <- exp2gcn(hh, net_type="signed",  SFTpower=power, cor_method="pearson")
net <- exp2gcn_1(hh,include_TOM = T,module_merging_threshold = 1,
                 net_type="signed",  SFTpower=softpower, cor_method="pearson", min_module_size = 50)

tiff("coding_dendrogram_thph.tiff", compression = "lzw",width = 5, height = 5, units = "in", res=600)
plot_dendro_and_colors(net)
dev.off()

normal = rep(1, rownames(datExpr) %>% grep("N",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("N",.,value=T))

fl = rep(2, rownames(datExpr) %>% grep("FL",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("FL",.,value=T))

fh = rep(3, rownames(datExpr) %>% grep("FH",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("FH",.,value=T))

cs = rep(4, rownames(datExpr) %>% grep("CS",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("CS",.,value=T))

dl = rep(5, rownames(datExpr) %>% grep("DL",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("DL",.,value=T))

dh = rep(6, rownames(datExpr) %>% grep("DH",.,value=T) %>% length) %>% 
  setNames(rownames(datExpr) %>% grep("DH",.,value=T))


traitData = list(normal,fl,fh,cs,dl,dh) %>% lapply(function(x)  x %>% data.frame %>% 
            setNames("pathology")) %>% do.call(rbind,.) %>% .[colnames(hh),,drop=F]

tiff("coding_survival_module_thph.tiff", compression = "lzw",width = 7, height = 6, units = "in", res=600)
MEtrait <- module_trait_cor(exp=hh, MEs=net$MEs, cor_method="pearson", continuous_trait = T,
                            metadata = traitData ,
                            cex.text = 0.8, cex.lab.x = 0.8, cex.lab.y = 0.8)
dev.off()


tiff("coding_eigen_gene_thph.tiff", compression = "lzw",width = 5, height = 7, units = "in", res=600)
plot_eigengene_network(net)
dev.off()

tiff("coding_module_genes_thph.tiff", compression = "lzw",width = 6, height = 5, units = "in", res=600)
plot_ngenes_per_module(net)
dev.off()

################# Eigen genes to plot
for (i in net$MEs %>% names){
toplot = net$MEs[,i,drop=F] %>% cbind(traitData[rownames(.),,drop=F]) %>% 
  arrange(pathology)
correlation = MEtrait$cor[MEtrait$ME == i] %>% signif(2)
pval = MEtrait$pvalue[MEtrait$ME == i] %>% signif(2)

bp = barplot( toplot %>% .[,1] , 
              col=rep(c("green","yellow","orange","royalblue","red","brown"), toplot[,2] %>% table),
             ylab = "", main = 
               paste0(gsub("ME","",i),"\nEigenGene with pathology(cor: ",correlation,",pvalue: ",pval,")"))

legend("bottom",legend = c("N","FL","FH","CS","DL","DH"),bty = "n",
       fill = c("green","yellow","orange","royalblue","red","brown"),
       border = "black", horiz = T)

}

#### Module wise KEGG pathway enrichment

###### pathway enrichment

net$genes_and_modules %>% group_split(.$Modules,.keep = F) %>% 
  lapply(function(x) x %>% mutate(symbol = gene_mapping[Genes]) %>% as.data.frame) %>%
  setNames(net$genes_and_modules %>% group_split(.$Modules,.keep = F) %>%
             lapply(function(z) z$Modules %>% unique ) %>% unlist ) %>%
  write.xlsx(file = "normals_wgcna_modules.xlsx")
  
