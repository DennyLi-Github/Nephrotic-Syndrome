
library(tidyverse)


library(org.Hs.eg.db)

library(annotate)

library(DESeq2)
library(gProfileR)
library(knitr)
library(stats)
library(gprofiler2)


library(GO.db)

BiocManager::install("GOstats")
library(GOstats)

library(pathview)
library(gage)

BiocManager::install("gageData")
library(gageData)


library(clusterProfiler)

load("Three_groups.RData")

#FSGS <- read.delim("GSE200828_fsgs.top.table.tsv")
#MCD <- read.delim("GSE200828_mcd.top.table.tsv")
#MN <- read.delim("GSE200828_mn.top.table.tsv")


FSGS_0.01 <- FSGS %>% filter(adj.P.Val < 0.01)
MCD_0.01 <- MCD %>% filter(adj.P.Val < 0.01)
MN_0.01 <- MN %>% filter(adj.P.Val < 0.01)


FSGS_MCD <- inner_join(FSGS_0.01, MCD_0.01, by = "ID")
Final <- inner_join(FSGS_MCD, MN_0.01, by = "ID")

final_separate <- Final %>% separate(ID, into = c("Entrez_ID", "AT"))

final_separate$minFC <- apply(final_separate[, c(7, 12, 17)], 1, min, na.rm = TRUE)

#final_separate1 <- final_separate[, c(1, 13, 17)]


sig_lfc=0.75




selectGenesUp <- final_separate[final_separate$minFC>sig_lfc, 'Entrez_ID']
selectGenesDown <-final_separate[final_separate$minFC<(-sig_lfc), 'Entrez_ID']

universeGenes <- unique(final_separate$Entrez_ID)
cutoff <- 0.01


upParams <- new("GOHyperGParams",
                geneIds=selectGenesUp,
                universeGeneIds=universeGenes,
                annotation="org.Hs.eg.db",
                ontology="CC",
                pvalueCutoff=cutoff,
                conditional=FALSE,
                testDirection="over")


downParams <- new("GOHyperGParams",
                geneIds=selectGenesDown,
                universeGeneIds=universeGenes,
                annotation="org.Hs.eg.db",
                ontology="CC",
                pvalueCutoff=cutoff,
                conditional=FALSE,
                testDirection="over")


upCC <- hyperGTest(upParams)
summary(upCC)[1:10,]

downCC <- hyperGTest(downParams)
summary(downCC)[1:10,]






upParams <- new("GOHyperGParams",
                geneIds=selectGenesUp,
                universeGeneIds=universeGenes,
                annotation="org.Hs.eg.db",
                ontology="MF",
                pvalueCutoff=cutoff,
                conditional=FALSE,
                testDirection="over")

upMF <- hyperGTest(upParams)
summary(upMF)[1:10,]

downParams <- new("GOHyperGParams",
                geneIds=selectGenesDown,
                universeGeneIds=universeGenes,
                annotation="org.Hs.eg.db",
                ontology="MF",
                pvalueCutoff=cutoff,
                conditional=FALSE,
                testDirection="over")

downMF <- hyperGTest(downParams)
summary(downMF)[1:10,]




foldchanges <- final_separate$minFC
names(foldchanges) <- final_separate$Entrez_ID
head(foldchanges)


data("go.sets.hs")
data("go.subs.hs")

goccsets <- go.sets.hs[go.subs.hs$CC]

goccres <- gage(exprs = foldchanges, gsets = goccsets, same.dir = TRUE)


data("kegg.sets.hs")
data("sigmet.idx.hs")

kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

keggres <- gage(exprs = foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)

head(keggres$greater) 

keggrespathways <- data.frame(id = rownames(keggres$greater), keggres$greater) %>% 
  tibble::as_tibble() %>% 
  filter(row_number() <= 20) %>% 
  .$id %>% 
  as.character()


keggrespathways

keggids <- substr(keggrespathways, start = 1, stop = 8)
keggids

tmp <- sapply(keggids, function(pid) pathview(gene.data = foldchanges, pathway.id = pid, species = "hsa"))







library(pathfindR)

df <- as.data.frame(final_separate)
colnames(df)[1] <- "Entrez_ID"


### Convert Entrez to Ensembl ids and Gene names
Ensembl_GeneNames <- bitr(df$Entrez_ID, fromType = "ENTREZID",
                          toType = c("SYMBOL"),
                          OrgDb = org.Hs.eg.db)

#names(final_separate)

final_separate1 <- final_separate %>% rename(ENTREZID = Entrez_ID)

Merge <- inner_join(final_separate1, Ensembl_GeneNames, by = "ENTREZID")


#enrichment_chart()

Merge1 <- Merge[, -c(1)]

str(Merge1)

col_order <- c("SYMBOL", "minFC", "adj.P.Val")
           
Merge2 <- Merge1[, col_order]

output_df <- run_pathfindR(Merge2, gene_sets = "GO-CC")



output_df <- run_pathfindR(Merge2, pin_name_path = "KEGG")




#############################################################

#Finding the key nodes for string analysis


library(tidyverse)

load("Three_groups.RData")

FSGS_0.01 <- FSGS %>% filter(adj.P.Val < 0.01)
MCD_0.01 <- MCD %>% filter(adj.P.Val < 0.01)
MN_0.01 <- MN %>% filter(adj.P.Val < 0.01)


FSGS_0.01<- FSGS_0.01 %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(200)
MCD_0.01<- MCD_0.01 %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(200)
MN_0.01<- MN_0.01 %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(200)



FSGS_MCD <- inner_join(FSGS_0.01, MCD_0.01, by = "ID")
common <- inner_join(FSGS_MCD, MN_0.01, by = "ID")


common_separate <- common %>% separate(ID, into = c("Entrez_ID", "AT"))


Membrane <- read.csv("Supplementary Table 5_Membrane_David.csv")



Membrane$Entrez_ID <- as.character(Membrane$Entrez_ID)

String <- inner_join(common_separate, Membrane, by = "Entrez_ID")




#########################################################################














