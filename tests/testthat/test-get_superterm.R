test_that("get_superterm does not zero-weight a term shared by every gene set in a community", {
  skip_if_not_installed("tm")

  # Reviewer's own example: a community whose gene sets all contain
  # "INTERFERON" -- exactly the word that should end up naming it. Under the
  # pre-fix standard TF-IDF (IDF = log2(N/d_t)), a term present in every
  # document of the corpus gets d_t == N, so IDF = log2(1) = 0 and the term
  # is dropped from the label entirely. A second, unrelated community is
  # included for a realistic multi-community call.
  geneset_names <- c(
    "INTERFERON_ALPHA_RESPONSE",
    "INTERFERON_GAMMA_RESPONSE",
    "TYPE_I_INTERFERON_SIGNALING",
    "INTERFERON_STIMULATED_GENES",
    "ANTIVIRAL_INTERFERON_DEFENSE",
    "GLYCOLYSIS_PATHWAY",
    "PENTOSE_PHOSPHATE_PATHWAY",
    "FRUCTOSE_METABOLISM"
  )
  community_membership <- c(1, 1, 1, 1, 1, 2, 2, 2)

  st <- get_superterm(geneset_names, community_membership, n_terms = 3)

  label_comm1 <- st$summary$superterm[st$summary$community == 1]
  terms_comm1 <- strsplit(label_comm1, "/")[[1]]

  expect_true("Interferon" %in% terms_comm1)
  expect_equal(terms_comm1[1], "Interferon")
})

test_that("get_superterm still returns a valid mapping/summary shape", {
  skip_if_not_installed("tm")

  geneset_names <- c(
    "INTERFERON_ALPHA_RESPONSE",
    "INTERFERON_GAMMA_RESPONSE",
    "GLYCOLYSIS_PATHWAY",
    "PENTOSE_PHOSPHATE_PATHWAY"
  )
  community_membership <- c(1, 1, 2, 2)

  st <- get_superterm(geneset_names, community_membership)

  expect_named(st, c("mapping", "summary"))
  expect_equal(nrow(st$mapping), length(geneset_names))
  expect_named(st$mapping, c("geneset", "community", "superterm"))
  expect_named(st$summary, c("community", "superterm", "n_genesets"))
  expect_setequal(st$mapping$community, c(1, 2))
})
