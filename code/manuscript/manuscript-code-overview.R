library(ggspavis)
library(STexampleData)
library(patchwork)
library(scater)

### code adapted from the ggspavis vignette ###
# https://bioconductor.org/packages/release/bioc/vignettes/ggspavis/inst/doc/ggspavis_overview.html

spe <- Visium_mouseCoronal()
rownames(spe) <- rowData(spe)$gene_name
colData(spe)$sum <- colSums(counts(spe))


plotVisium(spe, highlight = "in_tissue")


### code adapted from the banksy vignette ###
# https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/domain-segment.html

#### starmap sample ####

library(Banksy)

library(data.table)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scran)

#' Change paths accordingly
gcm_path <- "raw_data/well11processed_expression_pd.csv"
mdata_path <- "raw_data/well11_spatial.csv"

#' Gene cell matrix
gcm <- fread(gcm_path)
genes <- gcm$GENE
gcm <- as.matrix(gcm[, -1])
rownames(gcm) <- genes

#' Spatial coordinates and metadata
mdata <- fread(mdata_path, skip = 1)
headers <- names(fread(mdata_path, nrows = 0))
colnames(mdata) <- headers
#' Orient spatial coordinates
xx <- mdata$X
yy <- mdata$Y
mdata$X <- max(yy) - yy
mdata$Y <- max(xx) - xx
mdata <- data.frame(mdata)
rownames(mdata) <- colnames(gcm)

locs <- as.matrix(mdata[, c("X", "Y", "Z")])

#' Create SpatialExperiment
se <- SpatialExperiment(
  assay = list(processedExp = gcm),
  spatialCoords = locs,
  colData = mdata
)

lambda <- 0.8
k_geom <- 30
npcs <- 50
aname <- "processedExp"
se <- Banksy::computeBanksy(se, assay_name = aname, k_geom = k_geom)

set.seed(1000)
se <- Banksy::runBanksyPCA(se, lambda = lambda, npcs = npcs)

set.seed(1000)
se <- Banksy::clusterBanksy(se, lambda = lambda, npcs = npcs, resolution = 0.8)

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]

p <- plotColData(se, x = "X", y = "Y", point_size = 0.01, colour_by = cnames[1]) +
  scale_color_manual(values = pals::glasbey()) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none") 

#### visium sample ####

spe <- Visium_mouseCoronal()

# subset to only points in tissue
spe <- subset(spe, ,in_tissue == TRUE)

# normalise the sample 
spe <- logNormCounts(spe)

lambda <- 0.5
k_geom <- 6
npcs <- 50
aname <- "logcounts"
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs)

set.seed(1000)
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.8)

cnames <- colnames(colData(spe))
cnames <- cnames[grep("^clust", cnames)]

colData(spe) <- cbind(colData(spe),(spatialCoords(spe)))
colData(spe)

q <- plotColData(spe, x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", point_size = 1.5, colour_by = cnames[1]) +
  scale_y_reverse() +
  coord_equal() +
  scale_color_manual(values = pals::glasbey()) +
  theme_void() +
  theme(legend.position = "none") +
  xlab('X') +
  ylab('Y')

plotVisium(spe, annotate = cnames[1], image = TRUE) +
  scale_color_manual(values = pals::glasbey())

