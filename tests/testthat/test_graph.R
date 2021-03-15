context("Testing Network Distance")

test_that("creating a networkx graph from an adjacency matrix", {
  q <- matrix(c(1, 4, 3, 2, 5, 1, 3, 4, 7), nrow = 3, ncol = 3, byrow = TRUE)
  G <- igraph::graph_from_adjacency_matrix(q, mode = "directed")
  G <- igraph::simplify(G, remove.multiple = T, remove.loops = T)

  expect_equal(create_graph(q)[1], G[1])
  expect_equal(create_graph(q)[2], G[2])
  expect_equal(create_graph(q)[3], G[3])
})
