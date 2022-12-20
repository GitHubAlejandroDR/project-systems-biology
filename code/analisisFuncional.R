### esto creo que me sobra
args=commandArgs(trailingOnly = TRUE)
WD_code <- args[1]
WD_data <- args[2]
WD_results <- args[3]

# cargamos los .R necesarios
source("./comunidades.R")

### Funciones necesarias  ###

# para pasar de ensp id a entrez id
ensp2entrezid <- function(gene_uniprot){
  require(clusterProfiler)
  require(org.Hs.eg.db)
  gnm <- tryCatch(bitr(gene_uniprot,fromType="ENSEMBLPROT",toType="ENTREZID",OrgDb="org.Hs.eg.db"),error=function(e){NA})
  return(as.character(gnm[2]))
}

# para obtener los entrez ids pasando como parametro el num de cluster
obtenerEntrezID <- function(numcluster,data){
  ENSP_genes <- substr(getNodesIn(data,clusterids=c(numcluster)),6,20)
  genes_ENTREZ <- c()
  for(i in 1:length(ENSP_genes)){
    if(!is.na(ensp2entrezid(ENSP_genes[i]))){
      genes_ENTREZ <- c(genes_ENTREZ,as.character(ensp2entrezid(ENSP_genes[i])))
    }
  }
  return(genes_ENTREZ)
}

# para obtener los clusters con los entrez IDs de los genes que los forman. Obtenemos listas de listas que se usaran como parametro para realizar el analisis funcional
clusters_hpo <- function(data){
  lista <- list()
  for (i in 1:data$numbers[3]){
    name_clus <- as.name(paste('X',i, sep = ""))
    lista[as.character(name_clus)] <- list(c(obtenerEntrezID(i,data)))
  }
  return(lista)
}

# tenemos que filtrar los clusters que tengan genes semilla
# se le pasa o data_hashimoto o data_tiroiditis
clustersWithGeneSeed <- function(af,data){
  gene_seed <- data
  clusters_seleccionados <- c()
  for(i in 1:length(af@geneClusters)){
    for(j in 1:length(af@geneClusters[i][[1]])){
      for(l in 1:length(gene_seed)){
        print(af@geneClusters[i][[1]][j])
        print(as.character(gene_seed[l]))
        print(i)
        print(j)
        print(l)
        if(af@geneClusters[i][[1]][j] == as.character(gene_seed[l])){
          print("if ok")
          clusters_seleccionados <- append(clusters_seleccionados,af@compareClusterResult$Cluster[i])
        }
        
      }
    }
  }
  return(clusters_seleccionados)
}

##################

# ANÁLISIS FUNCIONAL A LOS CLUSTERS
af_hashimoto <- compareCluster(geneCluster = clusters_hpo(data_LC_hashimoto), fun = enrichKEGG, pvalueCutoff= 0.05)
af_tiroiditis <- compareCluster(geneCluster = clusters_hpo(data_LC_tiroiditis),fun=enrichKEGG,pvalueCutoff=0.05)

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

df <- same_pathway(af_hashimoto,af_tiroiditis)
df$hCluster_geneRatio <- sapply(df$hCluster_geneRatio, function(x) eval(parse(text=x)))
df$tCluster_geneRatio <- sapply(df$tCluster_geneRatio, function(x) eval(parse(text=x)))
df <- df[df$hCluster_geneRatio>=0.5,]
df <- df[df$tCluster_geneRatio>=0.5,]
setwd(WD_results)
write.csv(df,"comparacion_clusters.csv")
