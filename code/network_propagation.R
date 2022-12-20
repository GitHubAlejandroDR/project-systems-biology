# cargamos los .R necesarios
source("./variables.R")

# # creamos txt de los genes
# setwd(WD_data)

# ## hashimoto
# data_hashimoto <- read.csv("genes_for_HP_0000872(hashimoto).csv")
# data_hashimoto <- as.numeric(data_hashimoto$GENE_ENTREZ_ID)
# setwd(WD_results)
# lapply(data_hashimoto,write,"genes_hashimoto.txt",append=TRUE,ncolumns=1000) # lo guardamos en txt
# 
# ## tiroiditis
# setwd(WD_data)
# data_tiroiditis <- read.csv("genes_for_HP_0100646(tiroiditis).csv")
# data_tiroiditis <- as.numeric(data_tiroiditis$GENE_ENTREZ_ID)
# setwd(WD_results)
# lapply(data_hashimoto,write,"genes_tiroiditis.txt",append=TRUE,ncolumns=1000) # lo guardamos en txt
# 

setwd(WD_data)
string_file = "9606.protein.links.v11.5.txt"
score_threshold = 700

df = read.csv(string_file, sep = " ", header=TRUE)
# df = read.csv("network/9606.protein.links.v11.5.txt", sep = " ", header=TRUE)

# df2 = df[sample(1:nrow(df), 500),]
df2= dplyr::filter(df, combined_score >=score_threshold)
df2$p1 = sapply(df2$protein1, function(i) gsub("9606.", "", i))
df2$p2 = sapply(df2$protein2, function(i) gsub("9606.", "", i))
#borramos las 3 primeras columnas
borrar <- c("protein1","protein2","combined_score")
df2 <- df2[ , !(names(df2) %in% borrar)]
setwd(WD_results)
write.table(df2, file = "network_diamond_ensp.txt", col.names = F, quote = F, row.names = F, sep = ",")
setwd(WD_data)
df_full_h <- read.delim(file="string_protein_annotations_hashimoto.tsv")
df_full_t <- read.delim(file="string_protein_annotations_thyroiditis.tsv")
h_ensp <- as.vector(sapply(df_full_h$identifier, function(i) gsub("9606.", "", i)))
setwd(WD_results)
write.table(h_ensp, file = "hashimoto_ensp.txt", col.names = F, quote = F, row.names = F)
t_ensp <- as.vector(sapply(df_full_t$identifier, function(i) gsub("9606.", "", i)))
write.table(t_ensp, file = "thyroiditis_ensp.txt", col.names = F, quote = F, row.names = F)

# h_symbol <- df_full_h$X.node #hashimoto en symbol
# t_symbol <- df_full_t$X.node #tiroiditis en symbol


# DIAMOnD
setwd(WD_code)

## hashimoto
system("python3 DIAMOnD.py data/network_diamond_ensp.txt data/hashimoto_ensp.txt 200 ../results/hashimoto_propagated_200.txt")

## tiroiditis
system("python3 DIAMOnD.py data/network_diamond_ensp.txt data/thyroiditis_ensp.txt 200 ../results/thyroiditis_propagated_200.txt")

# Filtramos aquellos nodos añadidos que tengan un p-valor < 0.03
setwd(WD_results)
df_h <- read.delim(file="hashimoto_propagated_200.txt")

df_t <- read.delim(file="thyroiditis_propagated_200.txt")


