
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








#################################################################
dataset1 <- Merge %>% arrange(adj.P.Val1) %>% head(100)
dataset2 <- Merge %>% arrange(adj.P.Val2) %>% head(100)
dataset3 <- Merge %>% arrange(adj.P.Val3) %>% head(100)

merge_f <- inner_join(dataset1, dataset2, by = "SYMBOL")
merge_final <- inner_join(merge_f, dataset3, by = "SYMBOL")


Symbol <- merge_final[, c(10)]

Symbol_100 <- merge_final[, c(10)]

Symbol_100 <- as.data.frame(Symbol_100)

write.csv(Symbol, file = "Symbol.csv")


write.csv(Symbol_100, file = "Symbol_top100.csv")


################################################
library(tidyverse)

FSGS <- read.delim("GSE200828_fsgs.top.table.tsv")
MCD <- read.delim("GSE200828_mcd.top.table.tsv")
MN <- read.delim("GSE200828_mn.top.table.tsv")

library(tidyverse)

FSGS_positive <- FSGS %>% filter(logFC>0) %>% head(200)
MCD_positive <- MCD %>% filter(logFC>0) %>% head(200)
MN_positive <- MN %>% filter(logFC>0) %>% head(200)


Merge1 <- inner_join(FSGS_positive, MCD_positive, by = "ID")
Final <- inner_join(Merge1, MN_positive, by = "ID")


FSGS_0.05 <- FSGS %>% filter(adj.P.Val < 0.05)
MCD_0.05 <- MCD %>% filter(adj.P.Val < 0.05)
MN_0.05 <- MN %>% filter(adj.P.Val < 0.05)

#FSGS_0.05 <- FSGS_0.05 %>% arrange(desc(abs(logFC))) %>% filter(abs(logFC) > 2)
#MCD_0.05 <- MCD_0.05 %>% arrange(desc(abs(logFC))) %>% filter(abs(logFC) > 2)
#MN_0.05 <- MN_0.05 %>% arrange(desc(abs(logFC))) %>% filter(abs(logFC) > 2)


FSGS_0.05 <- FSGS_0.05 %>% arrange(desc(logFC)) %>% filter(logFC > 0) %>% head(200)
MCD_0.05 <- MCD_0.05 %>% arrange(desc(logFC)) %>% filter(logFC > 0) %>% head(200)
MN_0.05 <- MN_0.05 %>% arrange(desc(logFC)) %>% filter(logFC > 0) %>% head(200)



FSGS_MCD <- inner_join(FSGS_0.05, MCD_0.05, by = "ID")
Final <- inner_join(FSGS_MCD, MN_0.05, by = "ID")

write.csv(Final, file = "Final.csv")

final_separate <- Final %>% separate(ID, into = c("Entrez_ID", "AT"))



library(clusterProfiler)
library(org.Hs.eg.db)
### Convert Entrez to Ensembl ids and Gene names

df <- as.data.frame(final_separate)
colnames(df)[1] <- "Entrez_IDs"

Ensembl_GeneNames <- bitr(df$Entrez_IDs, fromType = "ENTREZID",
                          toType = c("SYMBOL"),
                          OrgDb = org.Hs.eg.db)

write.csv(Ensembl_GeneNames, file = "Ensembl_GeneNames.csv")

#Ensembl_GeneNames <- Ensembl_GeneNames %>% rename(Entrez_ID = ENTREZID)

#########################################################################


library(tidyverse)


Membrane <- read.csv("Membrane.csv")
#Top_200 <- read.csv("Ensembl_GeneNames_top200.csv")
Ensembl_GeneNames <- read.csv("Ensembl_GeneNames.csv")

Membrane <- Membrane %>% rename(ENTREZID = ID)

Merge <- inner_join(Membrane, Ensembl_GeneNames, by = "ENTREZID")




write.csv(Merge, file = "Merge.csv")






library(GO.db)

summary(GOBPANCESTOR)

GOTERM$"GO:0003700"

GOMFPARENTS$"GO:0003700"

GOMFCHILDREN$"GO:0003700"



GOTERM$"GO:0003677"

GOMFPARENTS$"GO:0003677"

GOMFCHILDREN$"GO:0003677"


GOMFPARENTS$"GO:0016021"

GOMFCHILDREN$"GO:0016020"


library(hgu95av2GO)

ll1 = hgu95av2GO[["39613_at"]]












