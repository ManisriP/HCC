##### WGCNA for 54 TH and PH tumour samples with clinical information using coding genes

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
sample_data2=inner_join(sample_data, aa,by = "Sample_ID") %>% mutate(Category = ifelse(grepl("Normal",.$Stage), "Normal" , 
                                                                                       ifelse( grepl("cancer",.$Stage) ,   "Tumour", "Non-tumour")   ))





##### Reading gene expression data and extracting coding features

xx = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,-2]

naming = read.csv("GSE148355_1.catholic.fpkm.RNAseqId.txt", sep="\t", header=T)[,c(1,2)]

gene_mapping = naming[,2] %>% setNames(naming[,1])

coding_features = gene_mapping %>% data.frame() %>% setNames("gene_symbol") %>% 
  subset(gene_symbol %in% (grch38 %>% filter(biotype == "protein_coding") %>% pull(symbol) %>% unique   )      ) %>% 
  rownames %>% sort

fpkm_data = xx %>% column_to_rownames("ID") %>% setNames(names(.) %>% gsub("Cat_|r","",.) %>% gsub("\\.","-",.))

##### WGCNA for 54 TH and PH tumour samples with clinical information using coding genes

thph_sampledata = bb %>%  inner_join(.,sample_data2 %>% dplyr::select(Sample_ID,Sample_GEO, Category), by = "Sample_ID") %>% 
  column_to_rownames("Sample_ID") %>%  mutate(OS = ifelse(.$`OS(months)` > 80, 2, 1)) %>%
  mutate(Surgery = ifelse(Operation == "TH", 2, 1)) %>% mutate(Grade = "G1")

tumour_stage = c(1:2, rep(3,3),4) %>% setNames(thph_sampledata$Tumor.stage %>% unique %>% mixedsort())

tumour_stage_colors = c("green","yellow",rep("orange",3),"red") %>% setNames(names(tumour_stage))

thph_sampledata$Tumor.stage.colors = tumour_stage_colors[thph_sampledata$Tumor.stage]

thph_sampledata$Tumor.stage = tumour_stage[thph_sampledata$Tumor.stage]

thph_sampledata[pData(g1)[,1][pData(g1)[,15]  == "tumor stage: G2"] %>% gsub(".*\\Cat_","",.) %>% gsub("]","",.) %>% 
                  intersect(rownames(thph_sampledata)),"Grade"] = "G2"

thph_sampledata[pData(g1)[,1][pData(g1)[,15]  == "tumor stage: G3"] %>% gsub(".*\\Cat_","",.) %>% gsub("]","",.) %>% 
                  intersect(rownames(thph_sampledata)),"Grade"] = "G3"


thph_fpkm_data = log2(fpkm_data[coding_features,thph_sampledata %>% rownames()] + 1) 

top_varying_genes = thph_fpkm_data %>% 
  as.matrix %>% rowVars() %>% order(decreasing = T) %>% .[1:5000]

thph_tumour_varying_genes = top_varying_genes

table(names(thph_fpkm_data) == rownames(thph_sampledata))

datExpr0 = thph_fpkm_data[top_varying_genes,] %>% t %>% data.frame()

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

hh=thph_fpkm_data[top_varying_genes,]

sft <- SFT_fit_1(exp=hh, net_type="signed", cor_method="pearson",rsquared = 0.8)

softpower <- sft$power

sft_thrsh=sft$dtf %>% mutate(special = ifelse(power == softpower,"yes","no") %>% factor())
mm = round(sft_thrsh$meank[sft_thrsh$power == softpower],1 )
ggpubr::ggarrange(sft$plot1,sft$plot2 + geom_hline(yintercept = mm, color="brown3") + 
                    geom_text(aes(0,mm, label=mm),hjust=0.1,vjust=-0.2,color="blue"))

if(!dir.exists("thph_tumour_samples/")){
  dir.create("thph_tumour_samples/")
}

setwd("thph_tumour_samples/")

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
                            metadata = thph_sampledata %>% .[colnames(hh),] %>% 
                              dplyr::select(c("OS(months)", "DFS(months)","AFP(ng/ml)", 
                                              "Maximum.size(cm)","Tumor.stage", "OS","Surgery")  ),
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

tiff("coding_hist_OS_thph.tiff", compression = "lzw",width = 4, height = 4, units = "in", res=600)
hist(thph_sampledata %>% .[colnames(hh),] %>% .$"OS(months)", main="", xlab = "Overall Survival of Samples(in months)")
dev.off()

tiff("coding_hist_DS_thph.tiff", compression = "lzw",width = 4, height = 4, units = "in", res=600)
hist(thph_sampledata %>% .[colnames(hh),] %>% .$"DFS(months)", main="", xlab = "Disease-free Survival of Samples(in months)")
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



### Survival plots

module_sheet_names=c("magenta","pink","yellow","black","brown","tan","blue")

dpi=600
w = 3
conversion=1/524*(dpi*w)/dpi*96

ee = vector("list",length(module_sheet_names)) %>% setNames(module_sheet_names)

for (i in module_sheet_names){

black_module = read.xlsx("thph_tumour_samples/coding_thph_wgcna_module_genes.xlsx",sheet=i)$Genes %>%
  discard(is.na)

traitData = thph_sampledata %>%dplyr::select("DFS(months)","Disease-free") %>% 
  setNames(c("DFS","DFS_status")) %>% mutate(DFS_status = ifelse(DFS_status == "Yes",0,1))

req_exp = thph_fpkm_data[black_module, rownames(traitData)]

eigen_exp = moduleEigengenes(expr = req_exp %>% t,
                             colors = rep("test", nrow(req_exp)) ) %>% .$eigengenes

survival_input = traitData %>% cbind(eigen_exp[rownames(.),,drop=F]) %>%
  mutate(MEtest = ifelse(MEtest < median(MEtest), 1,2 ) )

tiff(paste0("paper_figures/tumour_eigen_korean_",i,".tiff"), width = 3, height = 3, units = "in", res = 600,
     compression = "lzw")
p=survfit(Surv(DFS,DFS_status) ~ MEtest, data = survival_input) %>% 
  ggsurvplot(pval = T,pval.method = F, surv.median.line = "hv",pval.size = 4,pval.coord= c(17,0.27),
             censor.size = 3,
             legend.title= element_blank(),
             # title = paste0("Eigen gene (module ",module_nomen[i],")"),
             title = paste0("Module ",module_nomen[i]),
             xlab = "DFS(months)",
             legend = c(0.3,0.13),
             legend.labs =c(paste0("Low(n=",sum(survival_input$MEtest == 1),")"),
                            paste0("High(n=",sum(survival_input$MEtest == 2),")")),
             break.time.by = 10,
             ggtheme = theme_bw()+
               theme(panel.grid = element_blank(), 
                     axis.text.y = element_text(size=20*conversion, color="black"),
                     axis.text.x = element_text(size=19*conversion, color="black",angle=90,vjust = 0.5),
                     axis.title  = element_text(size=22*conversion, color="black"),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 20*conversion, color="black"),
                     plot.title = element_text(size = 22*conversion, color="black", hjust = 0.5),
                     panel.border = element_rect(color="black",linewidth = 0.7),
                     legend.key.size = unit(1,"pt"), legend.spacing.y = unit(-1,"pt")),
             palette = c("#00BFC4","#F8766D"),risk.table = F, cumevents = F )
print(p)
dev.off()
ee[[i]] = p
}