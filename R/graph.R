#' Utilities for creating and interacting with graph objects.
#'
#' @param A input MATRIX
#'
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix.
#' @param remove_self_loops If True, remove the diagonal of the matrix before creating the
#'                          graph object.
#'
#' @export
create_graph <- function(A, mode = NULL, remove_self_loops = T) {
  if (remove_self_loops) {
    diag(A) <- 0
  }

  if (is.null(mode)) {
<<<<<<< HEAD
    if (suppressWarnings(all(A == t(A), tolerance = (1e-05 * abs(t(A)) + 1e-08)))) {
=======
    if(suppressWarnings(all(A == t(A), tolerance = (1e-05 * abs(t(A)) + 1e-08)))){
>>>>>>> f7c42ad1971a9f7a310bb57b82247c8c100addbc
      G <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
      G <- igraph::simplify(G, remove.multiple = T, remove.loops = T)
    } else {
      G <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
      G <- igraph::simplify(G, remove.multiple = T, remove.loops = T)
    }
  } else {
    G <- igraph::graph_from_adjacency_matrix(A, mode = mode)
    G <- igraph::simplify(G, remove.multiple = T, remove.loops = T)
  }

  G
}
