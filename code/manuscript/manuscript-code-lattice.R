library(STexampleData)
library(ggspavis)
library(viridis)
library(RColorBrewer)
library(spdep)
library(scater)
library(tmap)
library(patchwork)
library(dplyr)
library(tidyr)

source("code/utils.R")

spe <- Visium_mouseCoronal()
rownames(spe) <- rowData(spe)$gene_name

####
#### Source: https://lmweber.org/BestPracticesST/chapters/workflow-Visium-mouseCoronal.html
# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

hist(colData(spe)$sum, breaks = 25)
qc_lib_size <- colData(spe)$sum < 600

# check spatial pattern of discarded spots
plotQC(spe,
  type = "spots",
  discard = "qc_lib_size"
)

# plot library size vs. number of cells per spot
hist(colData(spe)$detected, breaks = 20)
qc_detected <- colData(spe)$detected < 2000

# check spatial pattern of discarded spots
plotQC(spe,
  type = "spots",
  discard = "qc_detected"
)

# histogram of mitochondrial read proportions
hist(colData(spe)$subsets_mito_percent, breaks = 20)

# select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 28
table(qc_mito)

# check spatial pattern of discarded spots
plotQC(spe,
  type = "spots",
  discard = "qc_mito"
)

discard <- qc_lib_size | qc_detected | qc_mito

colData(spe)$discard <- discard

# remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)

spe <- logNormCounts(spe)

# plot
plotVisium(spe, fill = "Nrgn", assay = "logcounts", image = FALSE)


# add some values in 'colData' to annotate spots
colData(spe)$sum <- colSums(counts(spe))
colData(spe)$logsum <- log(colSums(counts(spe)))


plotSpots(spe, annotate = "detected") +
  scale_colour_gradient2() +
  geom_point(size = 1)

plotSpots(spe, annotate = "logsum") +
  scale_colour_viridis() +
  geom_point(size = 1.5)

plotSpots(spe, annotate = "detected") +
  scale_color_distiller(palette = "RdYlBu") +
  geom_point(size = 1.5)

plotSpots(spe, annotate = "logsum") +
  scale_color_distiller(palette = "RdYlBu") +
  geom_point(size = 1.5)

pLibSize <- plotSpots(spe, annotate = "sum") +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  geom_point(size = 1) +
  labs(color = "Library\nsize", title = "") +
  theme_void()

pDetected <- plotSpots(spe, annotate = "detected") +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  geom_point(size = 1) +
  labs(color = "Genes\ndetected", title = "") +
  theme_void()


plotSpots(spe, annotate = "sum") +
  scale_color_distiller(palette = "YlGn", direction = 1) +
  geom_point(size = 1)


plotVisium(spe, fill = "sum", trans = "log")

saveRDS(spe, "data/Visium_mouseCorona.rds")

# create sf object

coords <- as.data.frame(spatialCoords(spe))

coords$pxl_row_in_fullres <- -coords$pxl_row_in_fullres

spsf <- st_as_sf(coords,
  coords = c("pxl_col_in_fullres", "pxl_row_in_fullres")
)

Nrgn <- logcounts(spe)["Nrgn", ] |>
  as.matrix() |>
  as.data.frame()

colnames(Nrgn) <- "Nrgn"

spsf <- cbind(spsf, colData(spe), Nrgn)


plot(spsf)

# find 5 nearest neighbours
st_buffer(st_geometry(spsf), dist = 80) |>
  poly2nb(snap = 10) |>
  nb2listw() -> direct_neigbours

plot(st_buffer(st_geometry(spsf), dist = 70), border = "grey")
plot(direct_neigbours, st_geometry(spsf), add = TRUE)



glance_htest <- function(ht) {
  c(ht$estimate,
    "Std deviate" = unname(ht$statistic),
    "p.value" = unname(ht$p.value)
  )
}
# global moran's I
moran.test(spsf[["detected"]], listw = direct_neigbours, randomisation = FALSE) %>% glance_htest()

# local moran's I
loc <- localmoran(spsf[["detected"]], listw = direct_neigbours)
# extract the effect size
locEffect <- loc[, 1]
# extract the p-value and adjust for multiple testing
p.val.adj <- loc[, 5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

# plot(spsf)


pLocEff_direct <- ggplot() +
  geom_sf(data = spsf, aes(color = locEffect), size = 1) +
  scale_color_distiller(palette = "RdYlBu", limits = c(-2, 8)) +
  theme_void()

pLocEff_direct


# ggsave("misc/visium_overview_poster.pdf", width = 8, height = 6)


plot(st_buffer(st_geometry(spsf), dist = 70), border = "grey")
plot(direct_neigbours, st_geometry(spsf), add = TRUE)


#
spsf %>%
  knearneigh(k = 2) %>%
  knn2nb() %>%
  nb2listw(style = "W") -> knn5

glance_htest <- function(ht) {
  c(ht$estimate,
    "Std deviate" = unname(ht$statistic),
    "p.value" = unname(ht$p.value)
  )
}
# global moran's I
moran.test(spsf[["detected"]], listw = knn5, randomisation = FALSE) %>% glance_htest()

# local moran's I
loc <- localmoran(spsf[["detected"]], listw = knn5)
# extract the effect size
locEffect <- loc[, 1]
# extract the p-value and adjust for multiple testing
p.val.adj <- loc[, 5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

fill_range <- seq(-3, 8, by = 0.1)

pLocEff_knn <- ggplot() +
  geom_sf(data = spsf, aes(color = locEffect), size = 1) +
  scale_color_distiller(palette = "RdYlBu", limits = c(-2, 8)) +
  theme_void()

pLocEff_knn

####
# Find out neighbour size
####

st_coordinates(spsf) |> head()

# sort by X
st_coordinates(spsf) |>
  as.data.frame() |>
  arrange(X, Y) |>
  filter(X == 3137) |>
  mutate(Y.difference = Y - lag(Y))

st_coordinates(spsf) |>
  as.data.frame() |>
  arrange(Y, X) |>
  filter(Y == -8075) |>
  mutate(X.difference = X - lag(X))

st_buffer(st_geometry(spsf), dist = 70) |> plot()

st_buffer(st_geometry(spsf), dist = 70) |> poly2nb()


## gene
p <- pDetected + pLocEff_direct + pLocEff_knn
p


# Nrgn

pDetected <- plotVisium(spe,
  fill = "Nrgn", assay = "logcounts",
  image = FALSE, facets = NULL
) +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  geom_point(size = 5) +
  labs(color = "Nrgn", title = "") +
  theme_void() +
  labs(title = "", subtitle = "")

pDetected


spsf %>%
  knearneigh(k = 10) %>%
  knn2nb() %>%
  nb2listw(style = "W") -> knn5

# local moran's I
loc <- localmoran(spsf[["Nrgn"]], listw = knn5)
# extract the effect size
locEffect <- loc[, 1]
# extract the p-value and adjust for multiple testing
p.val.adj <- loc[, 5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

fill_range <- seq(-3, 8, by = 0.1)

pLocEff_knn <- ggplot() +
  geom_sf(data = spsf, aes(color = locEffect), size = 1) +
  scale_color_distiller(palette = "RdYlBu") +
  theme_void()

pLocEff_knn



st_buffer(st_geometry(spsf), dist = 80) |>
  poly2nb(snap = 10) |>
  nb2listw() -> direct_neigbours

# Moran's I
moran.test(spsf[["Nrgn"]], listw = direct_neigbours)


# local moran's I
loc <- localmoran(spsf[["Nrgn"]], listw = direct_neigbours)
# extract the effect size
locEffect <- loc[, 1]
# extract the p-value and adjust for multiple testing
p.val.adj <- loc[, 5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

fill_range <- seq(-3, 8, by = 0.1)

pLocEff_direct <- ggplot() +
  geom_sf(data = spsf, aes(color = locEffect), size = 0.3) +
  scale_color_distiller(
    palette = "RdYlBu", direction = -1,
    limits = c(-4, 4)
  ) +
  labs(color = "locI(Nrgn)") +
  theme_void()

pLocEff_direct

pLocEff_pVal <- ggplot() +
  geom_sf(data = spsf, aes(color = -log10(p.val.adj)), size = 0.3) +
  scale_color_distiller(palette = "BuGn", direction = 1) +
  labs(color = "-log10(adj. p-val.)") +
  theme_void()

pLocEff_pVal

pLocEff_direct


spsf |>
  ggplot(aes(x = Nrgn, y = locEffect, fill = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  theme_light()


pGeneLoc <- spsf |>
  ggplot(aes(x = Nrgn, y = locEffect, fill = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  theme_light() +
  geom_density_2d(color = "cornflowerblue")


pLocPval <- spsf |>
  ggplot(aes(fill = Nrgn, x = locEffect, y = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_light() +
  geom_density_2d(color = "cornflowerblue")

psupploc <- pGeneLoc + pLocPval

## gene
p <- pDetected + pLocEff_direct + pLocEff_knn
p

p2 <- pDetected + pLocEff_direct + pLocEff_pVal

p2


## Example neighbourhood
##
##
(sfe <- HeNSCLCData())

# Empty cells
colData(sfe)$is_empty <- colData(sfe)$nCounts < 1
# Select, sum negative control probes
(neg_inds <- str_detect(rownames(sfe), "^NegPrb")) %>% sum()


colData(sfe)$prop_neg <- colSums(counts(sfe)[neg_inds, ]) / colData(sfe)$nCounts
# Remove low quality cells
sfe <- sfe[, !sfe$is_empty & sfe$prop_neg < 0.1]
# Re-calculate stats
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$is_neg <- neg_inds
# log Counts
sfe <- logNormCounts(sfe)

sfePoly <- sfe

colGraph(sfe, "knn10") <- findSpatialNeighbors(sfe,
  method = "knearneigh",
  dist_type = "idw", k = 10,
  style = "W"
)


colGraph(sfePoly, "poly2nb") <- colGeometries(sfePoly)[["cellSeg"]] |>
  poly2nb(snap = 2) |> # snap = 10
  nb2listw(zero.policy = TRUE)


bbox_use <- st_bbox(c(xmin = 8000, xmax = 15000, ymin = 161000, ymax = 166000))

sfe <- runUnivariate(sfe,
  features = c("KRT17"),
  type = "localmoran"
)

pKnn <- plotLocalResult(sfe, "localmoran",
  features = c("KRT17"), ncol = 1,
  colGeometryName = "cellSeg",
  bbox = bbox_use
)
pKnn <- pKnn + scale_fill_distiller(palette = "RdYlBu", limits = c(-13, 13)) +
  labs(title = "Local Moran's I\n(Contiguos Neighbours)", fill = "locI(KRT17)")

pKnn

# poly
sfePoly <- runUnivariate(sfePoly,
  features = c("KRT17"),
  type = "localmoran",
  zero.policy = TRUE
)


pPoly <- plotLocalResult(sfePoly, "localmoran",
  features = c("KRT17"), ncol = 1,
  colGeometryName = "cellSeg",
  bbox = bbox_use
)


pPoly <- pPoly + scale_fill_distiller(palette = "RdYlBu", limits = c(-13, 13)) +
  labs(title = "Local Moran's I\n(Contiguos Neighbours)", fill = "locI(KRT17)")

pPoly


pLocalDiff <- (pKnn + pPoly) +
  plot_layout(guides = "collect")

pLocalDiff


fig2 <- (pLocalDiff / p2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 10, hjust = 0, vjust = 0)
  )

fig2

ggsave("outs/localMoranOverview.pdf", fig2, width = 8, height = 6)


###

data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  pivot_longer(everything(), names_to = "neighbourhood", values_to = "Ii") |>
  ggplot(aes(x = Ii)) +
  geom_histogram() +
  facet_grid(neighbourhood ~ .)



pknncont <- data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  ggplot(aes(x = knn, y = cont)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "cornflowerblue") +
  theme_light() +
  labs(
    x = "Local Moran's I \n(10 Nearest Neighbours)",
    y = "Local Moran's I \n(Contiguos Neighbours)"
  ) +
  coord_fixed(ratio = 1) +
  geom_density_2d(color = "cornflowerblue")


data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  ggplot(aes(x = knn, y = cont)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "cornflowerblue") +
  theme_light() +
  labs(
    x = "Local Moran's I \n(10 Nearest Neighbours)",
    y = "Local Moran's I \n(Contiguos Neighbours)"
  ) +
  coord_fixed(ratio = 1) +
  geom_smooth(method = "lm", col = "red")


<<<<<<< Updated upstream
layout <- "
AAAA
AAAA
AAAA
BBCC
BBCC
BBCC
"

psupp <- pknncont + pGeneLoc + pLocPval +
  plot_layout(design = layout) +
=======
psupp <- pGeneLoc + pLocPval + pknncont + pDistcont +
>>>>>>> Stashed changes
  plot_annotation(tag_levels = "A")

psupp

ggsave("outs/localMoranOverviewSupp.pdf", psupp, width = 9, height = 6)
