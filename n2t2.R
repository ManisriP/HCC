##### WGCNA for paired normal and tumour samples of PH with clinical information using coding genes

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


source("modified_functions.R")

aa = read.xlsx("korean_cohort_table_1.xlsx", startRow=3, sheet=1)
g1 = getGEO(filename="GSE148355_series_matrix.txt.gz", AnnotGPL=F, GSEMatrix=T)

# ''' Of the 15 normal samples, one was collected from a metastatic cancer patient whounderwent PH, and the other 14 were collected from patients with
# cholangiocarcinoma or polycystic liver disease after histological confirmation. '''

sample_data = pData(g1) %>% mutate(title = title %>% gsub(".*Cat_","",.) %>% gsub("]", "",.) ) %>% 
  dplyr::select(names(.)[c(1,2,10,46,48,49,51)]) %>% setNames(c("Sample_ID","Sample_GEO","Patient_ID","Surgery","Age","Gender","Stage") )

bb = read.xlsx("korean_cohort_table_2.xlsx",sheet=1)

sample_data = pData(g1) %>% mutate(title = title %>% gsub(".*Cat_","",.) %>% gsub("]", "",.) ) %>% 
  dplyr::select(names(.)[c(1,2,10,46,48,49,51)]) %>% setNames(c("Sample_ID","Sample_GEO","Patient_ID","Surgery","Age","Gender","Stage") )

### Sample data2
sample_data2=inner_join(sample_data, aa,by = "Sample_ID") %>% mutate(Category = ifelse(grepl("Normal",.$Stage),
"Normal" , ifelse( grepl("cancer",.$Stage) ,   "Tumour", "Non-tumour")   ))





##### Reading gene expression data and extracting coding features

xx = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,-2]

naming = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,c(1,2)]

gene_mapping = naming[,2] %>% setNames(naming[,1])

coding_features = gene_mapping %>% data.frame() %>% setNames("gene_symbol") %>% 
  subset(gene_symbol %in% (grch38 %>% filter(biotype == "protein_coding") %>% pull(symbol) %>% unique   )      ) %>% 
  rownames %>% sort

fpkm_data = xx %>% column_to_rownames("ID") %>% setNames(names(.) %>% gsub("Cat_|r","",.) %>% gsub("\\.","-",.))

##### WGCNA for normal and tumour samples of PH with clinical information using coding genes

th_sampledata = bb %>%  inner_join(.,sample_data2 %>% dplyr::select(Sample_ID,Sample_GEO, Category), by = "Sample_ID") %>% 
  mutate(Grade = "G1")

tumour_stage = c(1:2, rep(3,3),4) %>% setNames(th_sampledata$Tumor.stage %>% unique %>% mixedsort())

tumour_stage_colors = c("green","yellow",rep("orange",3),"red") %>% setNames(names(tumour_stage))

th_sampledata$Tumor.stage.colors = tumour_stage_colors[th_sampledata$Tumor.stage]

th_sampledata$Tumor.stage = tumour_stage[th_sampledata$Tumor.stage]

th_sampledata[pData(g1)[,1][pData(g1)[,15]  == "tumor stage: G2"] %>% gsub(".*\\Cat_","",.) %>% gsub("]","",.) %>% 
                intersect(rownames(th_sampledata)),"Grade"] = "G2"

th_sampledata[pData(g1)[,1][pData(g1)[,15]  == "tumor stage: G3"] %>% gsub(".*\\Cat_","",.) %>% gsub("]","",.) %>% 
                intersect(rownames(th_sampledata)),"Grade"] = "G3"

th_sampledata = aa %>% filter(grepl("Partial",Operational_method)) %>% 
  mutate(Category = ifelse(grepl("ES|HCC",`Pathology(b)`) , "Tumour","Normal")) %>% 
  dplyr::select(Sample_ID,Category) %>% column_to_rownames("Sample_ID")

#th_sampledata["DH013","Category"] = "Tumour"

th_sampledata = th_sampledata[c(aa %>% filter(grepl("Partial",Operational_method ) & 
                                                grepl("Paired",`Paired_sample(a)`) & !grepl("ES|HCC",`Pathology(b)`)   ) %>% .[,1] , 
                                intersect(rownames(th_sampledata),bb$Sample_ID)),,drop=F]

th_sampledata = th_sampledata %>% mutate(grouping = ifelse(Category == "Tumour",2,1))

th_fpkm_data = log2(fpkm_data[coding_features,th_sampledata %>% rownames()] + 1) 

top_varying_genes = th_fpkm_data %>% 
  as.matrix %>% rowVars() %>% order(decreasing = T) %>% .[1:3500]

thph_tumour_varying_genes = top_varying_genes

table(names(th_fpkm_data) == rownames(th_sampledata))

datExpr0 = th_fpkm_data[top_varying_genes,] %>% t %>% data.frame()

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


hh=th_fpkm_data[top_varying_genes,]

sft <- SFT_fit_1(exp=hh, net_type="signed", cor_method="pearson",rsquared = 0.8)

softpower <- 14

sft_thrsh=sft$dtf %>% mutate(special = ifelse(power == softpower,"yes","no") %>% factor())
mm = round(sft_thrsh$meank[sft_thrsh$power == softpower],1 )
ggpubr::ggarrange(sft$plot1,sft$plot2 + geom_hline(yintercept = mm, color="brown3") + 
                    geom_text(aes(0,mm, label=mm),hjust=0.1,vjust=-0.2,color="blue"))

if(!dir.exists("ph_normal_tumour_samples/")){
  dir.create("ph_normal_tumour_samples/")
}

setwd("ph_normal_tumour_samples/")

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
net <- exp2gcn_1(hh,include_TOM = T,
                 net_type="signed",  SFTpower=softpower, cor_method="pearson", min_module_size = 100)

tiff("coding_survival_module_thph.tiff", compression = "lzw",width = 7, height = 6, units = "in", res=600)
MEtrait <- module_trait_cor(exp=hh, MEs=net$MEs, cor_method="pearson", continuous_trait = T,
                            metadata = th_sampledata %>% .[colnames(hh),] %>% 
                              dplyr::select(c("grouping")  ),
                            cex.text = 0.8, cex.lab.x = 0.8, cex.lab.y = 0.8)
dev.off()


tiff("coding_dendrogram_thph.tiff", compression = "lzw",width = 5, height = 5, units = "in", res=600)
plot_dendro_and_colors(net)
dev.off()

tiff("coding_eigen_gene_thph.tiff", compression = "lzw",width = 5, height = 7, units = "in", res=600)
plot_eigengene_network(net)
dev.off()

tiff("coding_module_genes_thph.tiff", compression = "lzw",width = 6, height = 5, units = "in", res=600)
plot_ngenes_per_module(net)
dev.off()

#### Extracting module-wise genes and their KEGG Pathway enrichment
hubs <- get_hubs_gcn(hh, net)
hubs = hubs %>% mutate(Symbol = gene_mapping[.$Gene])
hubs %>% write.table("hub_genes.txt",sep="\t",col.names = T,row.names = F,quote = F)

net$moduleColors %>% unique %>% setNames(net$moduleColors %>% unique) %>% 
  lapply(function(z)  net$genes_and_modules %>% subset(Modules == z) %>% mutate(symbol = gene_mapping[Genes])  %>% 
           dplyr::select(-Modules) %>% qpcR:::cbind.na(.,enrichR::enrichr(.$symbol, 
                                                                          databases = "KEGG_2021_Human")$KEGG_2021_Human ) )  %>%  
  append(list(hub_genes = hubs %>% qpcR:::cbind.na(.,enrichR::enrichr(.$Symbol, 
                                                                      databases = "KEGG_2021_Human")$KEGG_2021_Human ) ) ) %>% 
  write.xlsx("coding_thph_wgcna_module_genes.xlsx")





