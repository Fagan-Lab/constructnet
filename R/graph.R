#graph.R
#--------
#status: Finished draft and some simple tests
#' Utilities for creating and interacting with graph objects.
#author: Stefan McCabe (stefanmccabe at gmail dot com)
#converted by: Zhaoyi Zhuang

#' @param A input MATRIX
#'
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix.
#' @param remove_self_loops If True, remove the diagonal of the matrix before creating the
#'                          graph object.
#'
#' @export
create_graph <- function(A, mode = NULL, remove_self_loops=T) {
    #Flexibly creating a networkx graph from a adjacency matrix

    #Returns
    #-------
    #G
    #   A graph
  if(remove_self_loops) {
    diag(A) = 0
  }
  if (is.null(mode)) {
    if(suppressWarnings(all(A == t(A), tolerance = (1e-05 * abs(t(A)) + 1e-08)))){
      G <- graph_from_adjacency_matrix(A, mode = "undirected")
      G <- simplify(G, remove.multiple = T, remove.loops = T)
    } else {
      G <- graph_from_adjacency_matrix(A, mode = "directed")
      G <- simplify(G, remove.multiple = T, remove.loops = T)
    }
  } else {
    G <- graph_from_adjacency_matrix(A, mode = mode)
    G <- simplify(G, remove.multiple = T, remove.loops = T)
  }

  G
}

#test
# q= matrix(c(1, 4, 3, 2, 5, 1, 3, 4, 7), nrow=3, ncol=3, byrow = TRUE)
# net <- create_graph(q)
# plot(net)

