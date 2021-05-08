library(igraph)
library(threejs)

# Read a matrix
metagenomes <- read.table("/Users/cdiaz/Career/ProfessionalTraining/1_HTS/UNAM-Networks-Workshop/Networks-Day4/Matrix/matrix_virus.txt", header = T, row.names =1, sep='\t')
head(metagenomes, 5)
# Creating the vertices for the network
viruses <- data.frame(id=1:nrow(metagenomes),
                      name=rownames(metagenomes),
                      shape=rep('circle',nrow(metagenomes)),
                      color=rep('blue',nrow(metagenomes)),
                      label.color=rep('black',nrow(metagenomes)))

samples <- data.frame(id=1:ncol(metagenomes)+nrow(viruses),
                      name=colnames(metagenomes),
                      shape=rep('square',ncol(metagenomes)),
                      color=rep('red',ncol(metagenomes)),
                      label.color=rep('black',ncol(metagenomes)))

vertices <- rbind(viruses,samples)
vertices$size = rep(8,nrow(vertices))
vertices$label.color = rep('black',nrow(vertices))

# Creating the edges for the network
edges <- data.frame(id1=NA, id2=NA, width=NA)[numeric(0),]
for (coln in 1:ncol(metagenomes)) {
  for (rown in 1:nrow(metagenomes)) {
    if (metagenomes[rown,coln]>0) {
      edges[nrow(edges)+1,] = list(id1=rown, id2=coln+nrow(metagenomes), width=metagenomes[rown, coln])
    }
  }
}

edges$color = rep('black',nrow(edges))

# Checking the network plot with edge attributes
metagen_net <- graph_from_data_frame(d=edges, v=vertices, directed = FALSE)
plot(metagen_net)

# Community detection
kc = fastgreedy.community(metagen_net)
sizes(kc)
membership(kc)
plot(kc, metagen_net)

# Interactive networks
inet <- graphjs(metagen_net, vertex.size = 1)
i <-  membership(kc)
sizes(kc)
g <- set_vertex_attr(metagen_net, "color", value = c("yellow", "blue", "red", 'orange', 'cyan', 'grey', 'black','green')[i])
graphjs(g)



# Likelihood of getting a random graph looking the same as our metagenomes networks
random <- erdos.renyi.game(n = gorder(metagen_net), p.or.m = ed, type = "gnp")
plot(random)
ed2 <- edge_density(random)
mean_distance(random, directed = FALSE)

# Generate 1000 random graphs
gl <- vector('list', 1000)
for(i in 1:1000){
  gl[[i]] <- erdos.renyi.game(n = gorder(metagen_net), p.or.m = ed, type = "gnp")}
gl.apls <- unlist(lapply(gl, mean_distance, directed = FALSE))
hist(gl.apls, xlim = range(c(1.5, 10)))
abline(v = metagen_net.apl, col = "red", lty = 3, lwd = 2)
mean(gl.apls > metagen_net.apl)

# Optional visualization
m <- layout_as_tree(metagen_net)
plot(metagen_net, vertex.label.color = "black", layout = m)


# Some characteristics of the network
V(metagen_net)
E(metagen_net)
gsize(metagen_net) # number of edges
gorder(metagen_net) # number of vertices
vertex_attr(metagen_net)
edge_attr(metagen_net)
ed <- edge_density(metagen_net)
farthest_vertices(metagen_net) # The vertices that are farthest apart
get_diameter(metagen_net) # Path sequence between the 2 most separated vertices
mean_distance(metagen_net, directed = FALSE)

# Degrees (importance of the nodes)
outDegree <- degree(metagen_net, mode = c('out'))
table(outDegree)
hist(outDegree, breaks = 30)
which.max(outDegree)

# Betweenness (also measures the importance of a node)
between <- betweenness(g, directed = FALSE)
hist(between, breaks = 30)
plot(metagen_net, 
     vertex.label = NA,
     edge.color = 'black',
     vertex.size = sqrt(between)+1,
     edge.arrow.size = 0.05,
     layout = layout_nicely(metagen_net))

# Centrality
central <- eigen_centrality(metagen_net)
which.max(central$vector)
plot(metagen_net,
     vertex.label.color = "black", 
     vertex.label.cex = 0.6,
     vertex.size = 25*(central$vector),
     edge.color = 'gray88',
     main = "Metagenomes Network")

# Visualizing important nodes and edges
PelagibacterPh_HTVC008M <- make_ego_graph(g, diameter(g), nodes = 'Pelagibacter phage HTVC008M', mode = c("all"))[[1]]
dists <- distances(PelagibacterPh_HTVC008M, "Pelagibacter phage HTVC008M") # geodesic distances of all vertices from vertex PelibacterPhage
colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
V(PelagibacterPh_HTVC008M)$color <- colors[dists+1]
plot(PelagibacterPh_HTVC008M, 
     vertex.label = dists, 
     vertex.label.color = "white",
     vertex.label.cex = .6,
     edge.color = 'black',
     vertex.size = 7,
     edge.arrow.size = .05,
     main = "Geodesic Distances from Pelagibacter Phage HTVC008M")

# Neighborhood
n1 <- neighbors(metagen_net, 'Guymas', mode = c('all'))
n2 <- neighbors(metagen_net, 'SRR127461', mode = c('all'))
n3 <- neighbors(metagen_net, 'SRR2046238', mode = c('all'))
n4 <- neighbors(metagen_net, 'ERR2021511', mode = c('all'))
n5 <- neighbors(metagen_net, 'SRR3577362', mode = c('all'))
intersection(n1,n2,n3,n4,n5)
intersection(n1,n2,n3,n4)
intersection(n1,n3,n5)


