##if (!requireNamespace("BiocManager", quietly=TRUE))
##install.packages("BiocManager")
##BiocManager::install(c("DEP", "pathwayPCA", "clusterProfiler", 
##                       "pathview", "KEGGgraph", "dplyr", "limma"))
BiocManager::install("topGO")
##install.packages("dplyr")


library(DEP)
library(dplyr)
library(limma)
library(SummarizedExperiment)
library(Glimma)
library(sva)
library(gplots)
library(statmod)
library(ggplot2)
library(pathwayPCA)
library(tidyverse)

setwd("/Users/kankanitdoungkamchan/Documents/Boar_proteomics/Summer_boars_rep1/Proteomics")

data <- read.delim("proteinGroups.txt") 
colnames(data)
data <- filter(data, Reverse != "+", Potential.contaminant != "+")
dim(data)

data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

data_unique <- make_unique(data, "Protein.IDs", "Majority.protein.IDs", delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers

#get name column and data columns

data1_1 <- data_unique[,c(315:316, LFQ_columns)] 
data1 <- data1_1[,!grepl("3131",colnames(data1_1))]

  
colnames(data1_1)

temp2 <- colnames(data1)[-c(1:2)] %>% strsplit(".", fixed = TRUE) 

#make targets 

targets <- data.frame(PigID = sapply(temp2, function(x) x[3]) %>% as.numeric(),
                      day = sapply(temp2, function(x) x[5]) %>% factor(levels = c(0,6,17,24,60)))

#read in treatment and sib info

# pig_info <- read.csv("/Users/kankanitdoungkamchan/Downloads/Label_pig.csv")
pig_info <- read.csv("Label_pig.csv")

pig_info$replicate <- c(1,1,2,4,2,3, 3)

#join together

targets <- left_join(targets, pig_info)

#addin label, condition columns

targets$condition <- paste(targets$Treatment, targets$day, sep ="_")
targets$label <- paste(targets$condition, targets$PigID, sep = "_")

colnames(data1)[-c(1:2)] <- targets$label

#filter data if not found in at least 3 samples
i.filter <- rowSums(data1[,-c(1:2)] >0)
table(i.filter)

data2 <- data1[i.filter >= 3, ] 
summary(data2)

temp <- data2[,-c(1:2)]

min(data2[,-c(1:2)]>0)

data_se <- make_se(data2, columns = c(3:ncol(data2)), expdesign =  targets)
plot_frequency(data_se)

#data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_se)
plot_missval(data_se)

data_filt <- filter_missval(data_se, thr = 0)

plot_coverage(data_filt)

data_norm <- normalize_vsn(data_filt)

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

#Do Glimma MDS clustering


glMDSPlot(assay(data_imp), top = 1500, labels = data_imp$label,
          groups = targets)

contrasts <- c("heat_0_vs_ctrl_0", "heat_6_vs_ctrl_6", "heat_17_vs_ctrl_17", 
               "heat_24_vs_ctrl_24", "heat_60_vs_ctrl_60",
               
               "heat_60_vs_heat_0", "heat_24_vs_heat_0", "heat_17_vs_heat_0",
               "heat_6_vs_heat_0",
               
               "ctrl_60_vs_ctrl_0", "ctrl_24_vs_ctrl_0", "ctrl_17_vs_ctrl_0",
               "ctrl_6_vs_ctrl_0")

data_diff <- test_diff(data_imp, type = "manual", test=contrasts)

# adding in rejection criteria
dep <- add_rejections(data_diff, alpha = 0.10, lfc = log2(1.20))

# plot top two pcas.
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 3, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))


# volcano plots of contrasts
for (contrast in contrasts) {
  print(plot_volcano(dep, contrast = contrast, label_size = 2, add_names = TRUE))
}

# summary object
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()



library(KEGGREST)
#### Obtain list of S. scrofa pathways
ssc_kegg_descriptions <- keggList("pathway", "ssc")

# Shorten descriptions
ssc_kegg_descriptions <- sub(" - Sus scrofa \\(pig\\)", "", ssc_kegg_descriptions)

# Turn the identifiers into KEGG-style pathway identifiers (org_id#####)
names(ssc_kegg_descriptions) <- sub("path:", "", names(ssc_kegg_descriptions))

#### Obtain and parse genes per each pathway
ssc_kegg_genes <- sapply(names(ssc_kegg_descriptions), function(pwid){
  pw <- keggGet(pwid)
  pw <- pw[[1]]$GENE[c(FALSE, TRUE)] # get gene symbols, not descriptions
  pw <- sub(";.+", "", pw) # discard any remaining description
  pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] # remove any mistaken lines that cannot be gene symbols
  pw <- unique(pw) # keep unique symbols
  return(pw)
})

#### Filter terms to exclude those with 0 genes (metabolic pathways)
ssc_kegg_genes <- ssc_kegg_genes[sapply(ssc_kegg_genes, length) != 0]
ssc_kegg_descriptions <- ssc_kegg_descriptions[names(ssc_kegg_descriptions) %in% 
                                                 names(ssc_kegg_genes)]

## Downloading the STRING PIN file to tempdir
url <- "https://stringdb-static.org/download/protein.links.v11.0/9823.protein.links.v11.0.txt.gz"
path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
download.file(url, path2file)

## read STRING pin file
ssc_string_df <- read.table(path2file, header = TRUE)

## filter using combined_score cut-off value of 800
ssc_string_df <- ssc_string_df[ssc_string_df$combined_score >= 800, ]

## fix ids
ssc_string_pin <- data.frame(Interactor_A = sub("^9823\\.", "", ssc_string_df$protein1),
                             Interactor_B = sub("^9823\\.", "", ssc_string_df$protein2))
head(ssc_string_pin, 2)

library(biomaRt)
ssc_ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")
converted <- getBM(attributes = c("ensembl_peptide_id", "uniprot_gn_id"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(ssc_string_pin)),
                   mart = ssc_ensembl)
ssc_string_pin$Interactor_A <- converted$uniprot_gn_id[match(ssc_string_pin$Interactor_A, converted$ensembl_peptide_id)]
ssc_string_pin$Interactor_B <- converted$uniprot_gn_id[match(ssc_string_pin$Interactor_B, converted$ensembl_peptide_id)]
ssc_string_pin <- ssc_string_pin[!is.na(ssc_string_pin$Interactor_A) & !is.na(ssc_string_pin$Interactor_B), ]
ssc_string_pin <- ssc_string_pin[ssc_string_pin$Interactor_A != "" & ssc_string_pin$Interactor_B != "", ]

head(ssc_string_pin, 2)


# remove self interactions
self_intr_cond <- ssc_string_pin$Interactor_A == ssc_string_pin$Interactor_B
ssc_string_pin <- ssc_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
ssc_string_pin <- unique(t(apply(ssc_string_pin, 1, sort))) # this will return a matrix object

ssc_string_pin <- data.frame(A = ssc_string_pin[, 1],
                             pp = "pp",
                             B = ssc_string_pin[, 2])

path2SIF <- file.path(tempdir(), "sscrofaPIN.sif")
write.table(ssc_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
path2SIF <- normalizePath(path2SIF)
head(ssc_string_pin, n = 3)

colnames(ssc_string_pin) <- NULL
rownames(ssc_string_pin) <- NULL

head(ssc_string_pin)

### run pathfind R 
### loading data 

data2[, c(2:3)]

library(pathfindR)
heat17_vs_heat0_output <- run_pathfindR(input = data_results[, c(2,9)],
                                convert2alias = FALSE,
                                gene_sets = "Custom",
                                custom_genes = ssc_kegg_genes,
                                custom_descriptions = ssc_kegg_descriptions,
                                pin_name_path = path2SIF)

heat60_vs_heat0_output <- run_pathfindR(input = data_results[, c(2,28)],
                                        convert2alias = FALSE,
                                        gene_sets = "Custom",
                                        custom_genes = ssc_kegg_genes,
                                        custom_descriptions = ssc_kegg_descriptions,
                                        pin_name_path = path2SIF)


heat6_vs_heat_0_output <- run_pathfindR(input = data_results[, c("ID","heat_6_vs_heat_0_p.val")],
                                        convert2alias = FALSE,
                                        gene_sets = "Custom",
                                        custom_genes = ssc_kegg_genes,
                                        custom_descriptions = ssc_kegg_descriptions,
                                        pin_name_path = path2SIF)


ctrl60_vs_ctrl0_output <- run_pathfindR(input = data_results[, c("ID","ctrl_60_vs_ctrl_0_p.val")],
                                      convert2alias = FALSE,
                                      gene_sets = "Custom",
                                      custom_genes = ssc_kegg_genes,
                                      custom_descriptions = ssc_kegg_descriptions,
                                      pin_name_path = path2SIF)

ctrl6_vs_heat6_output <- run_pathfindR(input = data_results[, c("ID","heat_6_vs_ctrl_6_p.val")],
                                        convert2alias = FALSE,
                                        gene_sets = "Custom",
                                        custom_genes = ssc_kegg_genes,
                                        custom_descriptions = ssc_kegg_descriptions,
                                        pin_name_path = path2SIF)

ctrl0_vs_heat0_output <- run_pathfindR(input = data_results[, c("ID","heat_0_vs_ctrl_0_p.val")],
                                        convert2alias = FALSE,
                                        gene_sets = "Custom",
                                        custom_genes = ssc_kegg_genes,
                                        custom_descriptions = ssc_kegg_descriptions,
                                        pin_name_path = path2SIF)

ctrl17_vs_heat17_output <- run_pathfindR(input = data_results[, c("ID","heat_17_vs_ctrl_17_p.val")],
                                         convert2alias = FALSE,
                                         gene_sets = "Custom",
                                         custom_genes = ssc_kegg_genes,
                                         custom_descriptions = ssc_kegg_descriptions,
                                         pin_name_path = path2SIF)

ctrl24_vs_heat24_output <- run_pathfindR(input = data_results[, c("ID","heat_24_vs_ctrl_24_p.val")],
                                         convert2alias = FALSE,
                                         gene_sets = "Custom",
                                         custom_genes = ssc_kegg_genes,
                                         custom_descriptions = ssc_kegg_descriptions,
                                         pin_name_path = path2SIF)
