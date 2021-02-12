#' Creates a minimum spanning tree connecting the sensors.
#' @export
# distance <-distances(inverse.sqrt)
#'
#' @param inputArray inputArray
#' @param distance distance
#' @param ... ...
#'
#' @export
construct_correlation_spanning <- function(inputArray, distance, ...) {
  UseMethod("construct_correlation_spanning")
}

#' @export
construct_correlation_spanning.igraph <- function(inputArray, distance, ...) {
  # rowsinMatrix <- nrow(inputArray)
  # colsinMatrix <- ncol(inputArray)
  # correlation_matrix <- cor(inputArray, method = "pearson") #create correlation matrix
  #
  # if (distance == inverse.sqrt){ #Distance matrix
  #   D <- sqrt(2*(1-correlation_matrix))
  # } else{
  #   D <- 1 - square(correlation_matrix)
  # }
  #
  # MST <- mst(D)
  #
  # asSparse <- sparseMatrix(MST) #convert to sparse matrix
  # asAdjacency <- as_adjacency_matrix(asSparse)#convert to adjacency matrix
  #
  # G <- graph_from_adjacency_matrix(asAdjacency)#generate graph
}
