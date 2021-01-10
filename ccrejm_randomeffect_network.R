
library(Rfast)
library(igraph)
library(ggplot2)
library(gridExtra)

setwd('.../crejm-code')
load(".../crejm-code/postselection.est.RData")

Sigma1.final<-postselection.est$Sigma1.est

#----------------------------------------------------
# Now create the heatmaps

library(tidyverse)
library(lattice)
library(viridisLite)
library(RColorBrewer)
library(ggpubr)

#1. Player Random Effect Correlation Matrix
coul <- colorRampPalette(brewer.pal(8, "RdBu"))
h1<-levelplot(cov2cor(Sigma1.final),
              col.regions=c("#B2182B","#C13739","#D15748",
                            "#DF755D","#EC9374","#F5AF8E",
                            "#F9C6AD","#FADBC9","#E7E0DB",
                            "#D3E4ED","#B9D9E9","#9DCBE1",
                            "#7EB8D7","#5BA2CB","#3E8DC0",
                            "#2F79B6","#2166AC"))

plot(h1)

#2. Get the Network

library(igraph)
library(readxl)

network.1 <- read_excel(".../crejm-code/crejm_network.xlsx", 
                        sheet = "randeff_cov_network")

network.1<-as.matrix(network.1)

g <- graph.adjacency(network.1,mode="undirected",weighted=TRUE,diag=FALSE)
iso <- V(g)[degree(g)==0]
g <- delete.vertices(g, iso)
E(g)$color[E(g)$weight <0] <- 'red'
E(g)$color[E(g)$weight >0] <- 'blue'


# Plot the graph
plot.igraph(g,layout=layout_in_circle,vertex.frame.color="white",asp=0,
            vertex.color = "white",vertex.label.cex=1.5,edge.width=2)




#-------------------------------------------------------

