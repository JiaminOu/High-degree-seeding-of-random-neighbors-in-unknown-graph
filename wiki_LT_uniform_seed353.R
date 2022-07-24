#one hop function is corrected on 03/08/2021#
#set.seed() function changed to more random numbers on 25/08/2021
#for the same strategy, seeds will be different from networks structure#

library(igraph)
library(assertthat)

stopifnot(require(data.table))
stopifnot(require(Matrix))

sample_threshold <- function(nNode, prob, value,seed){
  threshold <- rep(value, ceiling(nNode * prob))
  L<-length(threshold)
  set.seed(seed)
  reorder <- sample(x = 1:L, size = L, replace = F)
  Rthreshold<-threshold[reorder]
  return(Rthreshold[1:nNode])
}


### calculate percentage of ones among neighbors for each individual. 

calculate_value <- function(node, infection, each_neighbors){
  return(mean(infection[each_neighbors[[node]]] == 1)) ### percentage with infection
}


####know x% ########
knowsome2 <- function(seed2, node_infor, xpercent, nNode, nSeed){
  #randomseq<-runif(1000,min=0,max=1000000)
  #set.seed(floor(randomseq[seed2]))
  seed2=seed2
  sample_seed <- sample(c(1:nNode), xpercent * nNode, replace = F)
  sample_infor <- node_infor[sample_seed, ]
  sort_degree <- sort(sample_infor$degree, index.return=TRUE, decreasing=TRUE) ## high-low
  sort_order <- sort_degree$ix[1:nSeed]
  node_seed<-sample_infor$id[sort_order]
  return(node_seed)
}

knowsome <- function(xpercent, seed2, nNode, nSeed,node_infor) {
  node_seed_matrix <- do.call('rbind', lapply(1:seed2, knowsome2, node_infor, xpercent, nNode, nSeed))
  # here do.call should be check in order to produce a correct matrix
  return(node_seed_matrix)
}


####know x% clustered; choose the neighbours of highest degrees########
knowsomeCluster2 <- function(seed2, node_infor, xpercent, nNode, nSeed, each_neighbors){
  #randomseq<-runif(1000,min=0,max=1000000)
  #set.seed(floor(randomseq[seed2]))
  seed2=seed2
  sample_seed <- sample(c(1:nNode), xpercent * nNode, replace = F)
  sample_infor <- node_infor[sample_seed, ]
  maxnode<-sample_infor$id[which.max(sample_infor$degree)]
  
  node_seed1<-as.numeric(each_neighbors[[maxnode]])
  
  if ( length(node_seed1)>=(nSeed-1)){
    sample_infor2 <- node_infor[as.numeric(node_seed1), ]
    sort_degree2 <- sort(sample_infor2$degree, index.return=TRUE, decreasing=TRUE) ## choose the neighbours with highest degrees#
    sort_order2 <- sort_degree2$ix[1:(nSeed-1)]
    node_seed<-sample_infor2$id[sort_order2]
    node_seed<-c(node_seed,maxnode)
  }else{
    node_seed<-c(node_seed1,maxnode)
    remain<-nSeed-length(node_seed1)
    
    #all neighbours of the selected seeds#
    node_seed1<-unique(unlist(lapply(1:length(node_seed), NEI,node_seed,each_neighbors)))
    node_seed1=node_seed1[!node_seed1 %in% node_seed]
    remain=nSeed-length(node_seed)
    
    while (length(node_seed1)>0 && length (node_seed1)<remain){
      node_seed<-unique(c(node_seed,node_seed1))
      #all neighbours of the selected seeds#
      node_seed1<-unique(unlist(lapply(1:length(node_seed), NEI,node_seed,each_neighbors)))
      node_seed1=node_seed1[!node_seed1 %in% node_seed]
      remain=nSeed-length(node_seed)
    }
    if (length(node_seed1)>0 && length (node_seed1)>=remain){
      sample_infor2 <- node_infor[as.numeric(node_seed1), ]
      sort_degree2 <- sort(sample_infor2$degree, index.return=TRUE, decreasing=TRUE) ## choose the neighbours with highest degrees#
      sort_order2 <- sort_degree2$ix[1:remain]
      node_seed2<-sample_infor2$id[sort_order2]
      node_seed<-c(node_seed,node_seed2)}
  }
  return(node_seed)
}

NEI<-function(nn,nodeIDlist,each_neighbors){
  node<-nodeIDlist[nn]
  return(as.numeric(each_neighbors[[node]]))
}

knowsomeCluster <- function(seed2,node_infor, xpercent, nNode, nSeed, each_neighbors) {
  node_seed_matrix <- do.call('rbind', lapply(1:seed2, knowsomeCluster2, node_infor, xpercent, nNode, nSeed,each_neighbors))
  # here do.call should be check in order to produce a correct matrix
  return(node_seed_matrix)
}

####random ########
randomseeds2<-function(seed2,nNode,nSeed){
  randomseq<-runif(1000,min=0,max=1000000)
  #set.seed(floor(randomseq[seed2]))
  seed2=seed2
  node_seed<-sample(c(1:nNode), nSeed, replace = F)
  return(node_seed)
}

randomseeds<-function(seed2,nNode,nSeed){
  node_seed_matrix <- do.call('rbind', lapply(1:seed2, randomseeds2, nNode,nSeed))
  return(node_seed_matrix)
}


####random cluster_neighbours are also random########
randomCluster2 <- function(seed2, nSeed,nNode,node_infor,each_neighbors){
  #randomseq<-runif(1000,min=0,max=1000000)
  #set.seed(floor(randomseq[seed2]))
  seed2=seed2
  randomnode<- sample(c(1:nNode), 1, replace = F)
  node_seed1<-as.numeric(each_neighbors[[randomnode]])
  length(node_seed1)>=(nSeed-1)
  
  if ( length(node_seed1)>=(nSeed-1)){
    node_seed1<-sample(as.numeric(each_neighbors[[randomnode]]),nSeed-1,replace = F)
    node_seed<-c(node_seed1,randomnode)
  }else {
    node_seed<-unique(c(node_seed1,randomnode))
    
    #all neighbours of the selected seeds#
    node_seed1<-unique(unlist(lapply(1:length(node_seed), NEI,node_seed,each_neighbors)))
    node_seed1=node_seed1[!node_seed1 %in% node_seed]
    remain=nSeed-length(node_seed)
    
    while (length(node_seed1)>0 && length (node_seed1)<remain){
      node_seed<-unique(c(node_seed,node_seed1))
      #all neighbours of the selected seeds#
      node_seed1<-unique(unlist(lapply(1:length(node_seed), NEI,node_seed,each_neighbors)))
      node_seed1=node_seed1[!node_seed1 %in% node_seed]
      remain=nSeed-length(node_seed)
    }
    if (length (node_seed1)>=remain){
      node_seed2<-sample(node_seed1,remain,replace = F)
      node_seed<-c(node_seed,node_seed2)}
  }
  return(node_seed)
}

randomCluster <- function(seed2,nSeed,nNode,node_infor,each_neighbors) {
  node_seed_matrix <- do.call('rbind', lapply(1:seed2, randomCluster2, nSeed,nNode,node_infor,each_neighbors))
  # here do.call should be check in order to produce a correct matrix
  return(node_seed_matrix)
}


#################one hop########################
onehopseeds3<-function(node,each_neighbors){
  #set.seed(seed)
  return(sample(each_neighbors[[node]], 1, replace = F))
}


onehopseeds2<-function(seed2,nNode,nSeed,each_neighbors){
  #randomseq<-runif(1000,min=0,max=1000000)
  #set.seed(floor(randomseq[seed2]))
  seed2=seed2
  node_seed_nei<-sample(c(1:nNode), nSeed*2, replace = F)
  node_seed <- unique(as.numeric(do.call('cbind', lapply(node_seed_nei, onehopseeds3, each_neighbors))))[1:nSeed]
  return(node_seed)
}

onehopseeds<-function(seed2,nNode,nSeed,each_neighbors){
  node_seed_matrix <- do.call('rbind', lapply(1:seed2, onehopseeds2, nNode,nSeed,each_neighbors))
  return(node_seed_matrix)
}



LT2<-function(i,node_seed_matrix,nNode,threshold,each_neighbors){
  node_seed=node_seed_matrix[i,]
  node_value <- rep.int(0, nNode) ### node_value
  infection <- rep.int(0, nNode)  ### infection status 0 uninfected 1 infected
  node_value[as.numeric((node_seed))] <- 1 ### assign value for seed nodes
  infection[as.numeric((node_seed))] <- 1 ### assign infection status for seed nodes
  new_infected <- list()
  
  ####
  day_infected <- rep(0,15) ### infected num daily
  day_infected[1]=sum(infection == 1)
  ####
  
  day <- 1
  max_day <- 15
  Diff=1
  
  while(day <= (max_day - 1) ){
    if (Diff==0){day_infected[(day+1):15]=day_infected[day]
    break}
    not_infected <- which(infection == 0)
    old_infected <- which(node_value > threshold)
    node_value[not_infected] <- unlist(lapply(not_infected, calculate_value,
                                              infection, each_neighbors))
    new_infected[[day+1]] <- setdiff(which(node_value > threshold), old_infected)
    infection[new_infected[[day+1]]] <- 1 ### if exceed threshold, infected
    day <- day + 1
    day_infected[day] <- sum(infection == 1)
    Diff<-day_infected[day]-day_infected[day-1]
  }
  return(day_infected)
}


LT <- function(i,node_seed_matrix,nNode,threshold,each_neighbors) {
  day_infected_matrix <- do.call('rbind', lapply(1:i, LT2,node_seed_matrix,nNode,threshold,each_neighbors))
  # here do.call should be check in order to produce a correct matrix
  return(day_infected_matrix)
}



#####sfExport##################################################

social_threshold_network1 <- function(nSeed,graph,
                                      seed2,#seed2 control the rounds of seeds selection#
                                      threshold_value,
                                      seeding_strategy){
  #build the network#
  #graph<-graph_from_edgelist(edge2, directed = FALSE)
  nNode=vcount(graph)
  
  adj_matrix <- igraph::as_adjacency_matrix(graph, type = 'both')
  
  each_neighbors <- which(adj_matrix > 0, arr.ind = TRUE)
  ### for every individual, find the connections
  each_neighbors <- split(each_neighbors[, 2], each_neighbors[, 1])
  
  degree_real<-rep(0,nNode) #multiple links are calculated in the degree by igraph#
  for (i in c(1:nNode)){
    degree_real[i]<-length(each_neighbors[[i]])
  }
  
  #############
  
  #use seed to control the assigned thresholds#
  
  node_infor <- data.table(id = 1:nNode, degree = degree_real)
  
  if(seeding_strategy =='know30'){ 
    node_seed_matrix<-knowsome(xpercent=0.3, seed2, nNode, nSeed,node_infor)
  } else if(seeding_strategy =='know15'){ 
    node_seed_matrix<-knowsome(xpercent=0.15, seed2, nNode, nSeed,node_infor)
  }else if(seeding_strategy =='knowall'){ 
    sort_degree <- sort(node_infor$degree, index.return=TRUE, decreasing=TRUE) ## high-low
    node_seed <- sort_degree$ix[1:nSeed]
    Vector1<-rep(node_seed,times=seed2)
    T2<-matrix(Vector1,nrow=nSeed)
    node_seed_matrix<-t(T2)
  }else if(seeding_strategy =='know30cluster'){ 
    node_seed_matrix<-knowsomeCluster(seed2,node_infor, xpercent=0.3, nNode, nSeed, each_neighbors)
  }else if(seeding_strategy =='know15cluster'){ 
    node_seed_matrix<-knowsomeCluster(seed2,node_infor, xpercent=0.15, nNode, nSeed, each_neighbors)
  }else if(seeding_strategy =='knowallcluster'){ 
    node_seed_matrix<-knowsomeCluster(seed2,node_infor, xpercent=1, nNode, nSeed, each_neighbors)
  }else if(seeding_strategy =='random'){ 
    node_seed_matrix<-randomseeds(seed2,nNode,nSeed)
  }else if(seeding_strategy =='randomcluster'){
    node_seed_matrix<-randomCluster(seed2,nSeed,nNode,node_infor,each_neighbors)
  }else if(seeding_strategy =='onehop'){ 
    node_seed_matrix<-onehopseeds(seed2,nNode,nSeed,each_neighbors) 
  }
  
  
  #######
  threshold <- rep(threshold_value,nNode) 
  
  n_active_matrix<-LT(seed2,node_seed_matrix,nNode,threshold,each_neighbors)
  
  return(n_active_matrix)
}


simulate_parallel <- function(Threshold,graph,seeding_strategy){
  return(social_threshold_network1(nSeed=353,
                                   graph=graph,
                                   seed2=10000,#seed2 control the rounds of seeds selection#
                                   threshold_value=Threshold,
                                   seeding_strategy=seeding_strategy))
}

