# cargamos los .R necesarios
source("./comunidades.R")

### Funciones necesarias  ###

# para pasar de ensp id a entrez id
ensp2entrezid <- function(data){
  require(biomaRt)
  for(i in 1:length(data)){
    if(nchar(data[i])==20){ # eilimina el prefijo de los ENSP 
      data[i] <- substr(data[i],6,20)
    }
  }
 
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","entrezgene_id"),values=data,mart=mart)
  
  return(G_list$entrezgene_id)
}
# 
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","entrezgene_id"),values=getNodesIn(data_LC_hashimoto,clusterids=c(1)),mart= mart)

# para obtener los entrez ids pasando como parametro el num de cluster
# obtenerEntrezID <- function(numcluster,data){
#   ENSP_genes <- substr(getNodesIn(data,clusterids=c(numcluster)),6,20)
#   genes_ENTREZ <- c()
#   for(i in 1:length(ENSP_genes)){
#     if(!is.na(ensp2entrezid(ENSP_genes[i]))){
#       genes_ENTREZ <- c(genes_ENTREZ,as.character(ensp2entrezid(ENSP_genes[i])))
#     }
#   }
#   return(genes_ENTREZ)
# }

# para obtener los clusters con los entrez IDs de los genes que los forman. Obtenemos listas de listas que se usaran como parametro para realizar el analisis funcional
clusters_hpo <- function(data){
  lista <- list()
  for (i in 1:data$numbers[3]){ #recorre el numero de clusters
    name_clus <- as.name(paste('X',i, sep = ""))
    lista[as.character(name_clus)] <- list(c(ensp2entrezid(getNodesIn(data,clusterids=c(i)))))
  }
  return(lista)
}

# tenemos que filtrar los clusters que tengan genes semilla
# se le pasa o data_hashimoto o data_tiroiditis
clustersWithGeneSeed <- function(af,gene_seed){
  clusters_seleccionados <- c()
  for(i in 1:length(af@geneClusters)){
    for(j in 1:length(af@geneClusters[i][[1]])){
      for(l in 1:length(gene_seed)){
        if(gene_seed[l] == af@geneClusters[i][[1]][j]) {
          clusters_seleccionados <- c(clusters_seleccionados,names(af@geneClusters[i]))
          print(as.character(i))
        }
      }
    }
  }
  return(unique(clusters_seleccionados))
}

# ANÁLISIS FUNCIONAL A LOS CLUSTERS
af_hashimoto <- compareCluster(geneCluster = clusters_hpo(data_LC_hashimoto), fun = enrichKEGG, pvalueCutoff= 0.05)
af_tiroiditis <- compareCluster(geneCluster = clusters_hpo(data_LC_tiroiditis),fun=enrichKEGG,pvalueCutoff=0.05)

# seleccionamos las filas que  pertenezcan a un cluster con gene seed
af_tiroiditis@compareClusterResult<- af_tiroiditis@compareClusterResult[af_tiroiditis@compareClusterResult$Cluster ==clustersWithGeneSeed(af_tiroiditis,ensp2entrezid(list_semT$FALSE.)),]
af_hashimoto@compareClusterResult<- af_hashimoto@compareClusterResult[af_hashimoto@compareClusterResult$Cluster ==clustersWithGeneSeed(af_hashimoto,ensp2entrezid(list_semH$FALSE.)), ]



##################






# funcion que devuelve un dataframe con aquellos clusters que comparten pathway biologico en ambos hpo.
same_pathway <- function(afh,aft){
  hashimoto_clusters <-c()
  tiroiditis_clusters <- c()
  hCluster_pathway <- c()
  tCluster_pathway <- c()
  hCluster_geneRatio <- c()
  tCluster_geneRatio <- c()
  hCluster_pvalue <- c()
  tCluster_pvalue <- c()
  for(i in 1:length(afh@compareClusterResult$Cluster)){
    for(j in 1:length(aft@compareClusterResult$Cluster)){
      if(afh@compareClusterResult$Description[i] == aft@compareClusterResult$Description[j] ){
        hashimoto_clusters <- append(hashimoto_clusters,afh@compareClusterResult$Cluster[i])
        tiroiditis_clusters <- append(tiroiditis_clusters,aft@compareClusterResult$Cluster[j])
        hCluster_pathway <- append(hCluster_pathway,afh@compareClusterResult$Description[i])
        tCluster_pathway <- append(tCluster_pathway,aft@compareClusterResult$Description[j])
        hCluster_geneRatio <- append(hCluster_geneRatio,afh@compareClusterResult$GeneRatio[i])
        tCluster_geneRatio <- append(tCluster_geneRatio,aft@compareClusterResult$GeneRatio[j])
        hCluster_pvalue <- append(hCluster_pvalue,aft@compareClusterResult$p.adjust[i])
        tCluster_pvalue <- append(tCluster_pvalue,aft@compareClusterResult$p.adjust[j])
        
      }
    }
  }
  df <- data.frame(hashimoto_clusters,tiroiditis_clusters,hCluster_pathway,tCluster_pathway,hCluster_geneRatio,hCluster_pvalue,tCluster_geneRatio,tCluster_pvalue)
  return(df)
}

# aplciamos restricciones 
df <- same_pathway(af_hashimoto,af_tiroiditis)
df$hCluster_geneRatio <- sapply(df$hCluster_geneRatio, function(x) eval(parse(text=x))) #paso el gene ratio de fraccion a decimal
df$tCluster_geneRatio <- sapply(df$tCluster_geneRatio, function(x) eval(parse(text=x))) # same

# df <- df[df$hCluster_geneRatio>=0.3,]
df <- df[df$tCluster_geneRatio>=0.7,]
setwd(WD_results)
write.csv(df,"comparacion_clusters.csv")
