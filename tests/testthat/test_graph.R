library(igraph)
library(testthat)


q= matrix(c(1, 4, 3, 2, 5, 1, 3, 4, 7), nrow=3, ncol=3, byrow = TRUE)
G <- graph_from_adjacency_matrix(q, mode = "directed")
G <- simplify(G, remove.multiple = T, remove.loops = T)

test_that("creating a networkx graph from an adjacency matrix", {
  expect_equal(create_graph(q)[1], G[1])
  expect_equal(create_graph(q)[2], G[2])
  expect_equal(create_graph(q)[3], G[3])


})