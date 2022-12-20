# cargamos los .R necesarios
source("./variables.R")

# creamos txt de los genes
setwd(WD_data)

## hashimoto
data_hashimoto <- read.csv("genes_for_HP_0000872(hashimoto).csv")
data_hashimoto <- as.numeric(data_hashimoto$GENE_ENTREZ_ID)
setwd(WD_results)
lapply(data_hashimoto,write,"genes_hashimoto.txt",append=TRUE,ncolumns=1000) # lo guardamos en txt

## tiroiditis
setwd(WD_data)
data_tiroiditis <- read.csv("genes_for_HP_0100646(tiroiditis).csv")
data_tiroiditis <- as.numeric(data_tiroiditis$GENE_ENTREZ_ID)
setwd(WD_results)
lapply(data_hashimoto,write,"genes_tiroiditis.txt",append=TRUE,ncolumns=1000) # lo guardamos en txt



# DIAMOnD
setwd(WD_code)

## hashimoto
system("python3 DIAMOnD.py data/network_diamond_ensp.txt data/hashimoto_ensp.txt 200 ../results/hashimoto_propagated_200.txt")

## tiroiditis
system("python3 DIAMOnD.py data/network_diamond_ensp.txt data/thyroiditis_ensp.txt 200 ../results/thyroiditis_propagated_200.txt")

# Filtramos aquellos nodos añadidos que tengan un p-valor < 0.03
setwd(WD_results)
df_h <- read.delim(file="hashimoto_propagated_200.txt")
df_h <- df_h[df_h$p_hyper<=0.03,]
df_t <- read.delim(file="thyroiditis_propagated_200.txt")
df_t<- df_t[df_t$p_hyper<=0.03,]