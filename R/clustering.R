# Choosing and validating the number of subtypes.
#
# The original pipeline fixed k = 3 for k-means with no justification. Here k
# is chosen by consensus clustering (Monti et al., 2003) using the proportion
# of ambiguous clustering (PAC; Senbabaoglu et al., 2014), and reported with
# mean silhouette width and bootstrap stability (adjusted Rand index). Choices of
# k that produce a cluster smaller than 5% of the cohort are excluded.

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) sum(choose(x, 2))
  n <- sum(tab)
  idx <- comb2(tab); ai <- comb2(rowSums(tab)); bj <- comb2(colSums(tab))
  expected <- ai * bj / choose(n, 2)
  (idx - expected) / ((ai + bj) / 2 - expected)
}

pac_score <- function(consensus, lower = 0.1, upper = 0.9) {
  v <- consensus[lower.tri(consensus)]
  mean(v > lower & v < upper)
}

choose_k <- function(x, k_range = 2:8, reps = 200, p_item = 0.8, seed = 42, out_dir = NULL,
                     min_frac = 0.05) {
  x <- scale(x)
  d <- t(x)  # ConsensusClusterPlus clusters columns
  cc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    d, maxK = max(k_range), reps = reps, pItem = p_item, pFeature = 1, clusterAlg = "km",
    distance = "euclidean", seed = seed, plot = NULL, verbose = FALSE))
  dist_x <- stats::dist(x)
  res <- do.call(rbind, lapply(k_range, function(k) {
    cl <- cc[[k]]$consensusClass
    sil <- mean(cluster::silhouette(cl, dist_x)[, "sil_width"])
    data.frame(k = k, pac = pac_score(cc[[k]]$consensusMatrix), silhouette = sil,
               min_cluster = min(table(cl)))
  }))
  # A "subtype" of a handful of patients is not clinically usable: only k whose
  # smallest cluster holds >= max(10, min_frac * n) samples are eligible.
  res$eligible <- res$min_cluster >= max(10, ceiling(min_frac * nrow(x)))
  if (!any(res$eligible)) stop("no k gives clusters of a usable size")
  best <- res$k[res$eligible][which.min(res$pac[res$eligible])]
  res$selected <- res$k == best
  if (!is.null(out_dir)) utils::write.csv(res, file.path(out_dir, "k_selection.csv"), row.names = FALSE)
  list(table = res, k = res$k[res$selected], consensus = cc)
}

cluster_kmeans <- function(x, k, seed = 42) {
  set.seed(seed)
  stats::kmeans(scale(x), centers = k, nstart = 50, iter.max = 100)$cluster
}

bootstrap_stability <- function(x, k, n_boot = 100, seed = 42) {
  ref <- cluster_kmeans(x, k, seed)
  set.seed(seed)
  aris <- replicate(n_boot, {
    idx <- sample(nrow(x), replace = TRUE)
    km <- stats::kmeans(scale(x)[idx, , drop = FALSE], centers = k, nstart = 10)
    # assign every sample to its nearest bootstrap centre, compare with reference
    cen <- km$centers
    xs <- scale(x)
    lab <- apply(xs, 1, function(r) which.min(colSums((t(cen) - r)^2)))
    adjusted_rand(ref, lab)
  })
  c(mean = mean(aris), lo = unname(stats::quantile(aris, 0.025)), hi = unname(stats::quantile(aris, 0.975)))
}

# Agreement between clusters and histology: if ARI with LUAD/LUSC is high, the
# "subtypes" are mostly histology.
histology_agreement <- function(clusters, histology) {
  ok <- !is.na(histology)
  tab <- table(cluster = clusters[ok], histology = histology[ok])
  list(table = tab, ari = adjusted_rand(clusters[ok], histology[ok]),
       chisq_p = suppressWarnings(stats::chisq.test(tab)$p.value),
       cramers_v = sqrt(suppressWarnings(stats::chisq.test(tab)$statistic) / (sum(tab) * (min(dim(tab)) - 1))))
}
