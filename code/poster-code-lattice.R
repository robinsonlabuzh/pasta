library(STexampleData)
library(ggspavis)
library(viridis)
library(RColorBrewer)
library(spdep)
library(scater)
library(tmap)
library(patchwork)

spe <- Visium_mouseCoronal()

#### 
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
plotQC(spe, type = "spots", 
       discard = "qc_lib_size")

# plot library size vs. number of cells per spot
hist(colData(spe)$detected, breaks = 20)
qc_detected <- colData(spe)$detected < 2000

# check spatial pattern of discarded spots
plotQC(spe, type = "spots", 
       discard = "qc_detected")

# histogram of mitochondrial read proportions
hist(colData(spe)$subsets_mito_percent, breaks = 20)

# select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 28
table(qc_mito)

# check spatial pattern of discarded spots
plotQC(spe, type = "spots", 
       discard = "qc_mito")

discard <- qc_lib_size | qc_detected | qc_mito

colData(spe)$discard <- discard

# remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)


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
  scale_color_distiller(palette = "Spectral") +
  geom_point(size = 1.5)
  
plotSpots(spe, annotate = "logsum") +
  scale_color_distiller(palette = "Spectral") +
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

# create sf object

coords <- as.data.frame(spatialCoords(spe))

coords$pxl_row_in_fullres <- -coords$pxl_row_in_fullres 

spsf <- st_as_sf(coords,
         coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

spsf <- cbind(spsf, colData(spe))


plot(spsf)

# find 5 nearest neighbours
spsf %>% knearneigh(k = 6) %>% knn2nb() %>% 
  nb2listw(style = "W") -> knn6

glance_htest <- function(ht) c(ht$estimate, 
                               "Std deviate" = unname(ht$statistic), 
                               "p.value" = unname(ht$p.value))
# global moran's I
moran.test(spsf[['detected']], listw = knn6, randomisation = FALSE) %>% glance_htest()
 
# local moran's I
loc <- localmoran(spsf[['detected']], listw = knn6)
#extract the effect size
locEffect <- loc[,1]
#extract the p-value and adjust for multiple testing
p.val.adj <- loc[,5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

#plot(spsf)

pLocEff <- ggplot() + 
  geom_sf(data = spsf, aes(color = locEffect), size = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_void()

ggplot() + 
  geom_sf(data = spsf, aes(color = locEffect), size = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_void()

ggplot() + 
  geom_sf(data = spsf, aes(color = p.val.adj), size = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_void()

p <- pDetected + pLocEff

ggsave("misc/visium_overview_poster.pdf", width = 8, height = 6)



