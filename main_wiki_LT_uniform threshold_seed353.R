library(igraph)
edge<-read.csv("Wiki.csv",header=TRUE)
edge2<-as.matrix(edge[1:103688,1:2])
#modularity#############
wiki<-graph_from_data_frame(edge2, directed = FALSE)


components(wiki)   # more than 1 components
co <- components(wiki)
wiki1 <- induced.subgraph(wiki, which(co$membership == which.max(co$csize)))

graph<-wiki1


initial_parallel_process <- function(ncore){
  stopifnot(require(snowfall))
  sfInit(parallel = TRUE, cpus = ncore)
  sfLibrary(Matrix)
  sfLibrary(data.table)
  sfLibrary(igraph)
  sfLibrary(assertthat)
  sfExport("graph")
  sfSource('wiki_LT_uniform_seed353.R')
}


simulation_parallel <- function(seeding_strategy, Threshold){
  return(sfLapply(Threshold, simulate_parallel, graph=graph,seeding_strategy = seeding_strategy))
}


strategies <- c('knowall','know30','know15','knowallcluster','know30cluster','know15cluster',
                'random','randomcluster','onehop')

initial_parallel_process(28) 
simulate_strategies <- list()
start.time <- Sys.time()

threshold<-seq(0.15,0.625,by=0.025)

for(i in 1:length(strategies)){
  simulate_strategies[[i]] <- simulation_parallel(strategies[i], threshold)  
}

sfStop()
end.time <- Sys.time()
end.time - start.time
saveRDS(simulate_strategies,'wiki_LT_uniform_seed353.rds')    


