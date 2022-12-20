# cargamos los .R necesarios
source("./network_propagation.R")


# Instanciamos el constructor de la clase STRINGdb
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
string.network <- string_db$get_graph()


setwd(WD_results)


# mapeo
tiroiditis_mapped <- string_db$map( df_t, "DIAMOnD_node", removeUnmappedRows = TRUE )
hashimoto_mapped <- string_db$map( df_h, "DIAMOnD_node", removeUnmappedRows = TRUE )
#sacamos subred
subred_tiroiditis <- string_db$get_subnetwork(tiroiditis_mapped$STRING_id)
subred_hashimoto <- string_db$get_subnetwork(hashimoto_mapped$STRING_id)

red_hs <- graph_from_data_frame(df2)
df_subred_tiroiditis <- induced_subgraph(red_hs, df_t$DIAMOnD_node, impl = "auto")
df_subred_hashimoto <- induced_subgraph(red_hs, df_h$DIAMOnD_node, impl = "auto")

# deteccion de comunidades
df_subred_tiroiditis <- igraph::as_data_frame(df_subred_tiroiditis)
df_subred_hashimoto <- igraph::as_data_frame(df_subred_hashimoto)

data_LC_tiroiditis <- getLinkCommunities(df_subred_tiroiditis)
data_LC_hashimoto <- getLinkCommunities(df_subred_hashimoto)

## DA ERROR EN LOS PLOT PORQUE HAY 1 CLUSTER SOLO DE CADA FENOTIPO
setwd(WD_results)
png("comunidades_tiroiditis.png")
getClusterRelatedness(data_LC_tiroiditis, hcmethod = "ward.D")
dev.off()

png("comunidades_hashimoto.png")
getClusterRelatedness(data_LC_hashimoto, hcmethod = "ward.D")
dev.off()