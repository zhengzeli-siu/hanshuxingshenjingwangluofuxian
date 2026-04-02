file.path <- "~/Real_data/Other_alternatives/Result Analysis"
graph.path <- paste(file.path, "Graph", sep="")
library(igraph)
source(paste(file.path,"AAL.info.R",sep=""))

layout <- function(list, option=c("xy","xz","yz")){
  p <- length(list)
  result <- matrix(NA, nrow=p, ncol=2)
  coord.vec <- unlist(coordinate.list) # length is 3p
  x.vec <- rep(NA, p)
  y.vec <- rep(NA, p)
  z.vec <- rep(NA, p)
  
  for(j in 1:p){
    x.vec[j] <- coord.vec[3*(j-1) + 1]
    y.vec[j] <- coord.vec[3*(j-1) + 2]
    z.vec[j] <- coord.vec[3*(j-1) + 3]
  }
  
  if(option=="xy"){
    result[,1] <- x.vec
    result[,2] <- y.vec
  }
  if(option=="xz"){
    result[,1] <- x.vec
    result[,2] <- z.vec
  }
  if(option=="yz"){
    result[,1] <- y.vec
    result[,2] <- z.vec
  }
  
  return(result)
}

graph.plotter <- function(adj.mat, coordinate.list, name.tag.vec, 
                          option=c("xy","xz","yz"), file.path, title, color){
  p <- ncol(adj.mat)
  
  # compute number of neighbors of each node
  nbr.count <- colSums(adj.mat)
  
  net <- graph_from_adjacency_matrix(adj.mat, mode="undirected") %>%
    set_vertex_attr("label", value = name.tag.vec)
  
  plot(net, vertex.label=name.tag.vec, vertex.shape = "circle",
       vertex.color = "red", vertex.frame.color = "red",
       vertex.label.font = 1, vertex.label.family = "Helvetica",
       vertex.label.cex = 0.3, vertex.label.dist = 0.4, vertex.label.degree = pi/4,
       edge.color=color,  edge.width = 1.2, vertex.size = 0.5*(nbr.count+1),
       layout=layout(coordinate.list, "yz"), edge.curved = 0.8)
  return(0)
}

save.path <- paste(file.path, "Graphs", sep="/")

source(paste(file.path, "AdjMat.ABIDE.R", sep=""))
colnames(G.autism.and) <- name.tag.vec
row.names(G.autism.and) <- name.tag.vec
graph.plotter(G.autism.and, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.SCV.Autism", "black")
#dev.off()
graph.plotter(G.control.and, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.SCV.Control", "black")
#dev.off()

source(paste(file.path, "AdjMat.ADHD.R", sep=""))
graph.plotter(G.ADHD.or, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.SCV.ADHD", "black")
# dev.off()
graph.plotter(G.control.or, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.SCV.Control", "black")
# dev.off()
# 
source(paste(file.path, "AdjMat.ABIDE.spar.R", sep=""))
graph.plotter(G.autism.sym, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.two.pct.Autism", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.two.pct.Control", "darkblue")
# dev.off()
# 
source(paste(file.path, "AdjMat.ADHD.spar.R", sep=""))
graph.plotter(G.ADHD.sym, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.two.pct.ADHD", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.two.pct.Control", "darkblue")
# dev.off()

source(paste(file.path, "AdjMat.ADHD.spar.FGLasso.R", sep=""))
graph.plotter(G.ADHD.sym, coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.ADHD", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.Control", "darkblue")
# dev.off()


#Solea


solea_ADHD.austim=get(load("~/Real_data/npFGM/results/ADHD1_npFGM.RData"))
solea_ADHD.control=get(load("~/Real_data/npFGM/results/ADHD2_npFGM.RData"))
solea_ABIDE.austim=get(load("~/Real_data/npFGM/results/ABIDE1_npFGM.RData"))
solea_ABIDE.control=get(load("~/Real_data/npFGM/results/ABIDE2_npFGM.RData"))


one_percent_check <- function(mat) {
  diag(mat)=0
  one_count <- sum(mat == 1)
  total_elements <- length(mat)-dim(mat)[1]
  one_percentage <- one_count / total_elements
  return(one_percentage)
}


sapply(solea_ADHD.austim, one_percent_check)
sapply(solea_ADHD.control, one_percent_check)
sapply(solea_ABIDE.austim, one_percent_check)
sapply(solea_ABIDE.control, one_percent_check)

make_symmetric_and <- function(mat) {
  # Ensure the matrix is square
  if (nrow(mat) != ncol(mat)) {
    stop("Input matrix must be square")
  }
  
  # Get the matrix dimensions
  n <- nrow(mat)
  
  # Iterate through the upper triangular part of the matrix
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Apply the AND operation: both (i, j) and (j, i) must be 1, otherwise both are set to 0
      mat[i, j] <- mat[i, j] & mat[j, i]
      mat[j, i] <- mat[i, j]  # Ensure symmetry
    }
  }
  
  return(mat)
}

diag(solea_ADHD.austim[[3]])=0
diag(solea_ADHD.control[[1]])=0
diag(solea_ABIDE.austim[[4]])=0
diag(solea_ABIDE.control[[1]])=0


graph.plotter(solea_ADHD.austim[[3]], coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.ADHD", "darkgreen")
# dev.off()
graph.plotter(solea_ADHD.control[[1]], coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.Control", "darkgreen")
# dev.off()
graph.plotter(solea_ABIDE.austim[[4]], coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.ADHD", "darkgreen")
# dev.off()
graph.plotter(solea_ABIDE.control[[1]], coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.Control", "darkgreen")
# dev.off()








####################FGDNN######################

upper_triangle <- function(mat) {
  mat[upper.tri(mat, diag = TRUE)]
}


####################################
##############ADHD 1################
####################################
lambda.all=seq(1000,4000,10)

sparse=rep(NA,length(lambda.all))

d=read.csv("~/Real_data/fgDNN/results/ADHD_1.csv")
est=est.final=array(NA,c(length(lambda.all),p,p))

for(i in 1:p){
  for(j in 1:p){
    input_str <- d$selected[which((d$i==i)&(d$j==j))]
    
    
    clean_text <- gsub("tensor\\(|\\)", "", input_str)
    clean_text <- gsub("\\[|\\]", "", clean_text)
    elements <- strsplit(clean_text, ", ")[[1]]
    selected_num <- as.numeric(elements)
    selected_num=selected_num[-1]
    
    
    input_str <- d$lambda[which((d$i==i)&(d$j==j))]
    clean_text <- gsub("\\[|\\]", "", input_str)
    elements <- strsplit(clean_text, ", ")[[1]]
    lambda.value <- as.numeric(elements)
    lambda.value=lambda.value[-1]
    
    add=setdiff(lambda.all,lambda.value)
    selected_num=c(selected_num,rep(0,length(add)))
    
    
    
    est[,i,j]=selected_num>0
  }
}
for(r in 1:length(lambda.all)){  
  for(i in 1:p){
    for(j in 1:p){
      est.final[r,i,j]=all(est[r,i,j],est[r,j,i])
    }
  }
  est.final[r,,]=as.numeric(est.final[r,,])
  diag(est.final[r,,])=0
  sparse[r]=sum(est.final[r,,])/(length(est.final[r,,])-p)
}


graph.plotter(make_symmetric_and(est.final[300,,]), coordinate.list, name.tag.vec, "yz", save.path, "fgDNN.ADHD.two.pct.ADHD", "pink")













####################################
##############ADHD 2################
####################################
lambda.all=seq(1000,4000,10)

sparse=rep(NA,length(lambda.all))

d=read.csv("~/Real_data/fgDNN/results/ADHD_2.csv")
est=est.final=array(NA,c(length(lambda.all),p,p))

for(i in 1:p){
  for(j in 1:p){
    input_str <- d$selected[which((d$i==i)&(d$j==j))]
    
    
    clean_text <- gsub("tensor\\(|\\)", "", input_str)
    clean_text <- gsub("\\[|\\]", "", clean_text)
    elements <- strsplit(clean_text, ", ")[[1]]
    selected_num <- as.numeric(elements)
    selected_num=selected_num[-1]
    
    
    input_str <- d$lambda[which((d$i==i)&(d$j==j))]
    clean_text <- gsub("\\[|\\]", "", input_str)
    elements <- strsplit(clean_text, ", ")[[1]]
    lambda.value <- as.numeric(elements)
    lambda.value=lambda.value[-1]
    
    add=setdiff(lambda.all,lambda.value)
    selected_num=c(selected_num,rep(0,length(add)))
    
    
    
    est[,i,j]=selected_num>0
  }
}
for(r in 1:length(lambda.all)){  
  for(i in 1:p){
    for(j in 1:p){
      est.final[r,i,j]=all(est[r,i,j],est[r,j,i])
    }
  }
  est.final[r,,]=as.numeric(est.final[r,,])
  diag(est.final[r,,])=0
  sparse[r]=sum(est.final[r,,])/(length(est.final[r,,])-p)
}


graph.plotter(make_symmetric_and(est.final[300,,]), coordinate.list, name.tag.vec, "yz", save.path, "fgDNN.ADHD.two.pct.ADHD", "pink")







####################################
##############ABIDE 1################
####################################
lambda.all=seq(1000,6000,20)

sparse=rep(NA,length(lambda.all))

d=read.csv("~/Real_data/fgDNN/results/ABIDE_1.csv")
est=est.final=array(NA,c(length(lambda.all),p,p))

for(i in 1:p){
  for(j in 1:p){
    input_str <- d$selected[which((d$i==i)&(d$j==j))]
    
    
    clean_text <- gsub("tensor\\(|\\)", "", input_str)
    clean_text <- gsub("\\[|\\]", "", clean_text)
    elements <- strsplit(clean_text, ", ")[[1]]
    selected_num <- as.numeric(elements)
    selected_num=selected_num[-1]
    
    
    input_str <- d$lambda[which((d$i==i)&(d$j==j))]
    clean_text <- gsub("\\[|\\]", "", input_str)
    elements <- strsplit(clean_text, ", ")[[1]]
    lambda.value <- as.numeric(elements)
    lambda.value=lambda.value[-1]
    
    add=setdiff(lambda.all,lambda.value)
    selected_num=c(selected_num,rep(0,length(add)))
    
    
    
    est[,i,j]=selected_num>0
  }
}
for(r in 1:length(lambda.all)){  
  for(i in 1:p){
    for(j in 1:p){
      est.final[r,i,j]=all(est[r,i,j],est[r,j,i])
    }
  }
  est.final[r,,]=as.numeric(est.final[r,,])
  diag(est.final[r,,])=0
  sparse[r]=sum(est.final[r,,])/(length(est.final[r,,])-p)
}


graph.plotter(est.final[250,,], coordinate.list, name.tag.vec, "yz", save.path, "fgDNN.ADHD.two.pct.ADHD", "pink")






####################################
##############ABIDE 2################
####################################

lambda.all=seq(1000,8000,20)
sparse=rep(NA,length(lambda.all))

d=read.csv("~/Real_data/fgDNN/results/ABIDE_2.csv")
est=est.final=array(NA,c(length(lambda.all),p,p))

for(i in 1:p){
  for(j in 1:p){
    input_str <- d$selected[which((d$i==i)&(d$j==j))]
    
    
    clean_text <- gsub("tensor\\(|\\)", "", input_str)
    clean_text <- gsub("\\[|\\]", "", clean_text)
    elements <- strsplit(clean_text, ", ")[[1]]
    selected_num <- as.numeric(elements)
    selected_num=selected_num[-1]
    
    
    input_str <- d$lambda[which((d$i==i)&(d$j==j))]
    clean_text <- gsub("\\[|\\]", "", input_str)
    elements <- strsplit(clean_text, ", ")[[1]]
    lambda.value <- as.numeric(elements)
    lambda.value=lambda.value[-1]
    
    add=setdiff(lambda.all,lambda.value)
    selected_num=c(selected_num,rep(0,length(add)))
    
    
    
    est[,i,j]=selected_num>0
  }
}
for(r in 1:length(lambda.all)){  
  for(i in 1:p){
    for(j in 1:p){
      est.final[r,i,j]=all(est[r,i,j],est[r,j,i])
    }
  }
  est.final[r,,]=as.numeric(est.final[r,,])
  diag(est.final[r,,])=0
  sparse[r]=sum(est.final[r,,])/(length(est.final[r,,])-p)
}


graph.plotter(est.final[350,,], coordinate.list, name.tag.vec, "yz", save.path, "fgDNN.ADHD.two.pct.ADHD", "pink")
