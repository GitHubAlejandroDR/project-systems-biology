# #!/usr/bin/env Rscript
# # prepara la red de STRING para network propagation
# 
# args = commandArgs(trailingOnly=TRUE)
# library(dplyr)
# 
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("Necesita el archivo de red como argumento", call.=FALSE)
# } else if (length(args)==1) {
#   stop("Necesita umbral de score como argumento", call.=FALSE)
# }
setwd(WD_data)
string_file = "9606.protein.links.v11.5.txt"
score_threshold = 700

# string_file = "network/testSTRING_small.txt"
# score_threshold = 800
# 
# 
# ## funciones ####
# 
# # convierte entre ids de genes. Devuelve el id convertido de un gen
# id.converter <- function(from_id, from_type, to_type) {
#   require(clusterProfiler)
#   require(org.Hs.eg.db)
#   gnm <-tryCatch(
#     bitr(
#       from_id
#       , fromType=from_type
#       , toType=to_type
#       , OrgDb="org.Hs.eg.db"
#     )
#     ,error=function(e){NA}
#   )
#   return(gnm[1,2])
# }
# 
# 
# ensp2entrezid <- function(gene_uniprot) {
#   require(clusterProfiler)
#   require(org.Hs.eg.db)
#   gnm <-tryCatch(bitr(gene_uniprot, fromType="ENSEMBLPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db"),error=function(e){NA})
#   gnm <- tryCatch(if(!is.na(gnm)){gnm <- gnm$ENTREZID},error=function(e){NA})
#   
#   return(as.character(gnm[1]))
# }

############### funciones ##################


df = read.csv(string_file, sep = " ", header=TRUE)
# df = read.csv("network/9606.protein.links.v11.5.txt", sep = " ", header=TRUE)

# df2 = df[sample(1:nrow(df), 500),]
df2= dplyr::filter(df, combined_score >=score_threshold)
df2$p1 = sapply(df2$protein1, function(i) gsub("9606.", "", i))
df2$p2 = sapply(df2$protein2, function(i) gsub("9606.", "", i))
#borramos las 3 primeras columnas
borrar <- c("protein1","protein2","combined_score")
df2 <- df2[ , !(names(df2) %in% borrar)]
write.table(df2, file = "network_diamond_ensp.txt", col.names = F, quote = F, row.names = F, sep = ",")
df_full_h <- read.delim(file="string_protein_annotations_hashimoto.tsv")
df_full_t <- read.delim(file="string_protein_annotations_thyroiditis.tsv")
h_ensp <- as.vector(sapply(df_full_h$identifier, function(i) gsub("9606.", "", i)))
write.table(h_ensp, file = "hashimoto_ensp.txt", col.names = F, quote = F, row.names = F)
t_ensp <- as.vector(sapply(df_full_t$identifier, function(i) gsub("9606.", "", i)))
write.table(t_ensp, file = "thyroiditis_ensp.txt", col.names = F, quote = F, row.names = F)

h_symbol <- df_full_h$X.node #hashimoto en symbol
t_symbol <- df_full_t$X.node #tiroiditis en symbol
# 
# ## Conversion SYMBOL a ENTREZID
# library(org.Hs.eg.db)
# hs <- org.Hs.eg.db
# my.symbols <- as.vector(h_ensp)
# h_entrez <- mapIds(org.Hs.eg.db, keys = h_symbol, keytype="SYMBOL", column = "ENTREZID")
# 
# #convertir de ensp a symbol (pierde tambien mucha info)
# select(org.Hs.eg.db, keys=vector, columns="SYMBOL", keytype="ENSEMBLPROT")
# 
# #convertir de ssymbol a entrezid
# v_prueba <- c("ABL1","AKT1","CDKN1B","FLT3","GRB2","JAK2","LYN","MAPK1","RAC1","RHOA","SRC")
# mapIds(org.Hs.eg.db, keys = v_prueba, keytype="SYMBOL", column = "ENTREZID")
# 
# #convertir de ensp a symbol con biomart (hashimoto)
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list_hashimoto <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","entrezgene_id"),values=df_h$DIAMOnD_node,mart= mart)

# df2$entrezid_1 = sapply(df2$p1, ensp2entrezid)
# df2$entrezid_1 = as.character(df2$entrezid_1)
# df2$entrezid_2 = sapply(df2$p2, ensp2entrezid)
# df2$entrezid_2 = as.character(df2$entrezid_2)

# df2 = dplyr::filter(df2, (entrezid_1 != "character(0)" & entrezid_2 != "character(0)"))
# 
# out_DIAMOnD.df = dplyr::select(df2, entrezid_1, entrezid_2)
# write.table(out_DIAMOnD.df, "network_diamond.txt", sep = ",", quote = F, row.names = F, col.names = F)
# 
# out_guild.df = dplyr::select(df2, entrezid_1, entrezid_2)
# out_guild.df$int = rep("1", nrow(out_guild.df))
# out_guild.df = dplyr::select(out_guild.df, entrezid_1, int, entrezid_2)
# 
# write.table(out_guild.df, "network_guild.txt", sep = " ", quote = F, row.names = F, col.names = F)



