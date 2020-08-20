#Reconstruct a network from a random matrix, not taking time series into account.

#'@export
inputThreshold <- threshold(inputArray, min = -Inf, max = Inf)
construct_random <- function(inputArray, inputThreshold,...) UseMethod("construct_random")

#'@export
construct_random.igraph <- function(inputArray, inputThreshold,...){
rowsinMatrix <- nrow(inputArray)
colsinMatrix <- ncol(inputArray)
randomMatrices <-randn(rowsinMatrix, rowsinMatrix)
A<- threshold(randomMatrices,inputThreshold,...)
}
G<-graph_from_literal(A, simplify=FALSE)#allow for nodes