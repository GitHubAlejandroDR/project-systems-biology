# cargamos los .R necesarios
source("./network_propagation.R")


# Instanciamos el constructor de la clase STRINGdb
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
string.network <- string_db$get_graph()


setwd(WD_results)

# Añadimos genes semilla al resultado de la propagacion
table_propH <- read.delim(file="hashimoto_propagated_200.txt")
list_semH <- read.delim(file="hashimoto_ensp.txt", col.names = F)
df_h <- data.frame(sum_h=c(as.vector(table_propH$DIAMOnD_node),as.vector(list_semH$FALSE.)))

table_propT <- read.delim(file="thyroiditis_propagated_200.txt")
list_semT <- read.delim(file="thyroiditis_ensp.txt", col.names = F)
df_t <- data.frame(sum_t=c(as.vector(table_propT$DIAMOnD_node),as.vector(list_semT$FALSE.)))


# mapeo
tiroiditis_mapped <- string_db$map( df_t, "sum_t", removeUnmappedRows = TRUE )
hashimoto_mapped <- string_db$map( df_h, "sum_h", removeUnmappedRows = TRUE )

#sacamos subred
subred_tiroiditis <- string_db$get_subnetwork(tiroiditis_mapped$STRING_id)
subred_hashimoto <- string_db$get_subnetwork(hashimoto_mapped$STRING_id)


# deteccion de comunidades
df_subred_tiroiditis <- igraph::as_data_frame(subred_tiroiditis)
df_subred_hashimoto <- igraph::as_data_frame(subred_hashimoto)

data_LC_tiroiditis <- getLinkCommunities(df_subred_tiroiditis)
data_LC_hashimoto <- getLinkCommunities(df_subred_hashimoto)


setwd(WD_results)
png("comunidades_tiroiditis.png")
getClusterRelatedness(data_LC_tiroiditis, hcmethod = "ward.D")
dev.off()

png("comunidades_hashimoto.png")
getClusterRelatedness(data_LC_hashimoto, hcmethod = "ward.D")
dev.off()