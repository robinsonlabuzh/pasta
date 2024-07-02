suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatstat)
  library(mixR)
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(seg)
  library(stats)
  library(RANN)
  library(STexampleData)
  library(patchwork)
  library(reshape2)
  library(sf)
  library(ncf)
  library(rgeoda)
  library(Voyager)
  library(SpatialFeatureExperiment)
  library(SFEData)
  library(scran)
  library(scater)
  library(tmap)
  library(spdep)
  library(dixon)
  library(stringr)
  library(magrittr)
})

.ppp <- \(spe, marks = NULL) {
  xy <- spatialCoords(spe)
  w <- owin(range(xy[, 1]), range(xy[, 2]))
  m <- if (is.character(marks)) {
    stopifnot(
      length(marks) == 1,
      marks %in% c(
        gs <- rownames(spe),
        cd <- names(colData(spe))))
    if (marks %in% gs) {
      assay(spe)[marks, ]
    } else if (marks %in% cd) {
      spe[[marks]]
    } else stop("'marks' should be in ",
      "'rownames(.)' or 'names(colData(.))'")
  }
  ppp(xy[, 1], xy[, 2], window = w, marks = factor(m))
}


pValuesHotspotMarks <- function(pp, alpha = 0.05){
  # Code source: https://idblr.rbind.io/post/pvalues-spatial-segregation/
  
  # Significant p-values assumming normality of the Poisson process
  ## relrisk() computes standard errors based on asymptotic theory, assuming a Poisson process
  # call relative risk function
  f1 <- relrisk(pp_sel,se=TRUE)
  
  z <- qnorm(alpha/2, lower.tail = F)     # z-statistic
  f1$u <- f1$estimate + z*f1$SE           # Upper CIs
  f1$l <- f1$estimate - z*f1$SE           # Lower CIs
  mu_0 <- as.vector(table(spatstat.geom::marks(pp))/pp$n) # null expectations by type
  f1$p <- f1$estimate # copy structure of pixels, replace values
  for (i in 1:length(f1$p)) {
    f1$p[[i]]$v <- factor(ifelse(mu_0[i] > f1$u[[i]]$v, "lower",
                                ifelse( mu_0[i] < f1$l[[i]]$v, "higher", "none")),
                          levels = c("lower", "none", "higher"))
  }
  return(f1)
}

pValuesHotspot <- function(pp, alpha = 0.05){
  # Code source: https://idblr.rbind.io/post/pvalues-spatial-segregation/
  # density estimate for all marks
  f1 <- density(unmark(pp), se = TRUE)
  # Significant p-values assumming normality of the Poisson process
  z <- qnorm(alpha/2, lower.tail = F)     # z-statistic
  f1$u <- f1$estimate + z*f1$SE           # Upper CIs
  f1$l <- f1$estimate - z*f1$SE           # Lower CIs
  mu_0 <- intensity(unmark(pp)) # null expectations by type
  f1$p <- f1$estimate # copy structure of pixels, replace values
  f1$p$v <- factor(ifelse(mu_0 > f1$u$v, "lower",
                          ifelse( mu_0 < f1$l$v, "higher", "none")),
                   levels = c("lower", "none", "higher"))
  
  return(f1)
}

#' Functional Principal Component Analysis
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param r the functional domain
#' @param knots the number of knots
#' @param pve the proportion of variance explained
#'
#' @return a list with components of fpca.face
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' # calculate the Gcross metric for alpha and beta cells
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", r_seq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id"), ncores = 2
#' )
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_id
#' )
#' # extract the functional response matrix
#' mat <- metric_res %>%
#'     select(ID, r, rs) %>%
#'     spread(ID, rs) %>%
#'     select(!r)
#' # create a dataframe as required by pffr
#' dat <- data.frame(ID = colnames(mat))
#' dat$Y <- t(mat)
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, `[`, 1))
#' dat$patient_id <- factor(sapply(split_data, `[`, 2))
#' dat$image_id <- factor(sapply(split_data, `[`, 3))
#' # calculate fPCA
#' mdl <- functionalPCA(dat = dat, r = metric_res$r |> unique(), knots = 30, pve = 0.99)
#' @import dplyr
functionalPCA <- function(dat, r, knots, pve = 0.95) {
  # calculate the fPCA - this is a bit a pointless wrapper until now
  res <- refund::fpca.face(Y = dat$Y, center = TRUE, argvals = r, knots = knots, pve = pve)
  return(res)
}

#' Plot a biplot from an fPCA analysis
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param res the output from the fPCA calculation
#'
#' @return a list with components of fpca.face
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' # calculate the Gcross metric for alpha and beta cells
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", r_seq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id"), ncores = 2
#' )
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_id
#' )
#' # extract the functional response matrix
#' mat <- metric_res %>%
#'     select(ID, r, rs) %>%
#'     spread(ID, rs) %>%
#'     select(!r)
#' # create a dataframe as required by pffr
#' dat <- data.frame(ID = colnames(mat))
#' dat$Y <- t(mat)
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, `[`, 1))
#' dat$patient_id <- factor(sapply(split_data, `[`, 2))
#' dat$image_id <- factor(sapply(split_data, `[`, 3))
#' # calculate fPCA
#' mdl <- functionalPCA(dat = dat, r = metric_res$r |> unique(), knots = 30, pve = 0.99)
#' p <- plotFpca(dat = dat, res = mdl)
#' print(p)
#' @import dplyr
plotFpca <- function(dat, res, colourby) {
  scores_df <- res$scores %>% as.data.frame()
  # plot fCPA results - assumes same order of fPCA results and input data
  p <- ggplot(scores_df, aes(scores_df[, 1], scores_df[, 2], colour = (dat[[colourby]]))) +
    scale_color_continuous(type = "viridis", name = colourby) +
    geom_point(size=0.75) +
    coord_equal() +
    theme_light()
  return(p)
}

#PRE: list of point pattern, corresponding celltypes of interest, functions to evaluate
#POST: result of the metric
metricRes <- function(plot_by, pp, celltype, fun, bootstrap, continuous, f){
  if(continuous){
    pp <- subset(pp, select = plot_by)
  }
  else{
    pp_sel <- pp[[plot_by]]
    pp <- subset(pp_sel, marks == celltype)
  }
  if(bootstrap){
    metric.res <- lohboot(pp, fun = fun, f = f)
  }
  else{
    metric.res <- do.call(fun, args = list(X=pp, f=f))
  }
  metric.res$type <- celltype
  metric.res$plot_by <- plot_by
  return(metric.res)
}

#PRE: celltypes, function to calculation and edge correction method
#POST: dataframe of 
metricResToDF <- function(plot_by, celltype, pp, fun, edgecorr, bootstrap, continuous, f){
  lapply(plot_by, function(u) {
    metricRes(u, fun = fun, pp = pp, celltype = celltype, bootstrap, continuous, f)  %>%
      as.data.frame()
  }) %>% bind_rows
}

# [MR: try to write the above a little more compactly]

#PRE: Celltypes of interest, function to analyse, edge correction to perform
#POST: plot of the metric
plotMetric <- function(plot_by, pp, celltype = NULL, x, fun, edgecorr, bootstrap = FALSE, continuous = FALSE, f = NULL){
  #calculate the metric and store as dataframe
  res_df <- metricResToDF(plot_by, celltype, pp, fun, edgecorr, bootstrap, continuous, f)
  #plot the curves
  p <- ggplot(res_df, aes(x=.data[[x]], y=.data[[edgecorr]], colour= (plot_by),  show.legend = FALSE))+
    geom_line(size=1) +
    {if(bootstrap)geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, show.legend = FALSE)}+
    ggtitle(paste0(fun, '-function'))+
    geom_line(aes(x=.data[[x]],y=theo),linetype = "dashed", size=1)+
    ylab(edgecorr) +
    #guides(colour = guide_legend(override.aes = list(shape = 23, size = 2))) +
    theme_light()
  
  return(p)
}
