source("code/utils.R")
library('spatialFDA')
spe <- readRDS("data/spe.rds")

#subset the data to only look at sample ID 0.01, 0.06 and 0.26
# list("-0.29", "0.01", "0.06")
#zstack_list <- list("-0.04", '-0.09', '-0.14', '-0.19', '-0.24', '-0.29', '0.01', '0.06', '0.11', '0.16', '0.21', "0.26")

#define the Z-stacks that you want to compare
zstack_list <- list("-0.09", "0.01", "0.21")

#define the celltype that you want to compare across the stacks - hereby we assume independence across the z-stacks which is an assumption that can be challenged
celltype_ls <- "OD Mature"

selectZstacks <- function(zstack, spe){
  sub <- spe[, spe$sample_id == zstack]
  pp <- .ppp(sub, marks = "cluster_id")
  return(pp)
}

pp_ls <- lapply(zstack_list, selectZstacks, spe)
names(pp_ls) <- zstack_list

setRiprasWindows <- function(pp){
  Window(pp) <- ripras(pp)
  return(pp)
}
#the entire point patterns with the ripras windows
pp <- lapply(pp_ls, setRiprasWindows)

separateMarks <- function(pp){
  #split the multitype point process into several single type processes
  ppls <- split(pp)
  return (ppls)
}
#the point patterns separated by their marks
pp_ls <- lapply(pp, separateMarks)

# local L function
L_odmature_lisa <- localL(pp_ls[['0.01']]$`OD Mature`)

df <- as.data.frame(L_odmature_lisa)
dfm <- reshape2::melt(df, "r")

get_sel <- dfm %>% filter(r > 200.5630 & r < 201.4388, variable != "theo") %>%
  mutate(sel = value) %>% select(variable, sel)

dfm <- dfm %>% left_join(get_sel)

p <- ggplot(dfm, aes(x=r, y=value, group=variable, colour=sel)) +
  geom_line() + 
  scale_color_continuous(type = "viridis") +
  geom_vline(xintercept = 200) +
  theme(legend.position = "none") +
  theme_light() +
  ggtitle("LISA curves of slice 0.01") +
  labs(colour="Value at r~200")

ppdf <- as.data.frame(pp[['0.01']]) %>% filter(marks=="OD Mature")
ppdf$sel <- get_sel$sel # assume they are in same order

q <- ggplot(ppdf, aes(x=x, y=y, colour=sel)) + 
  geom_point(size=0.75) +
  scale_color_continuous(type = "viridis") +
  theme(legend.position = "none") +
  theme_light()+
  ggtitle("Points coloured by LISA value at r ~ 200")

# normalise the values
df_fdob <- asinh(df %>% as.matrix / 50) %>% as.data.frame()

# extract the functional response matrix
mat <- df_fdob %>%
  select(!c(r,theo))
# create a dataframe as required by pffr
dat <- data.frame(ID = colnames(mat))
dat$Y <- t(mat)
dat$sel <- get_sel$sel

# perform functional PCA
res <- functionalPCA(dat = dat, r = df_fdob$r |> unique())
scores_df <- res$scores %>% as.data.frame()
scores_df$sel <- get_sel$sel
# plot a biplot
p_biplot <- ggplot(scores_df, aes(scores_df[, 1], scores_df[, 2], colour = sel)) +
  geom_point() +
  coord_equal() +
  theme_light() +
  xlab("functional PC1") +
  ylab("functional PC2") +
  scale_color_continuous(type = "viridis") +
  ggtitle("Biplot of LISA curves") +
  xlab('PC1') +
  ylab('PC2') +
  labs(colour="Value at r~200")

res <- calcMetricPerFov(spe, 'OD Mature', subsetby = 'sample_id', fun = 'Lest', marks = 'cluster_id', r_seq=NULL, by = c('Animal_ID','sample_id'))
res <- subset(res, sample_id %in% c('-0.09', '0.01', '0.21'))

#assemble plots for paper
p_L <- plotMetricPerFov(res, theo = TRUE, correction = "iso", x = "r", imageId = 'sample_id') + theme(legend.position = 'right') + labs(colour="Slice")
p_local <- wrap_plots(list(p,p_biplot), nrow = 2, guides = 'collect')
L_plot <- wrap_plots(list(p_L, p_local), nrow = 2, heights = c(1, 2))


pls <- lapply(zstack_list, function(zstack){
  pp_sel <- pp_ls[[zstack]][[celltype_ls]]
  p <- pp_sel |> as.data.frame() |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(color = "black", alpha = 0.3, size=0.75) +
    theme_light() +
    coord_fixed() +
    ggtitle(paste0(zstack, ', (N = ',npoints(pp_sel), ')')) + 
    theme(legend.position = 'none')
  if(zstack == '0.01'){
    p <- pp_sel |> as.data.frame() |> 
      ggplot(aes(x = x, y = y, col = get_sel$sel)) +
      scale_color_continuous(type = "viridis") +
      geom_point(size=0.75) +
      theme_light() +
      coord_fixed() +
      ggtitle(paste0(zstack, ', (N = ',npoints(pp_sel), ')')) + 
      theme(legend.position = 'none')
  }
  return(p)
})

p1 <- wrap_plots(pls, guides = 'collect', nrow = 3)
p_total <- wrap_plots(list(p1, L_plot), ncol = 2) + plot_annotation(tag_levels = 'A')

ggsave('outs/pp_example.pdf', plot = p_total, width = 8, height = 9)

### code for figure of poster

spe <- readRDS("data/spe.rds")

#subset the data to only look at sample ID 0.01
sub <- spe[, spe$sample_id == 0.01]
(pp <- .ppp(sub, marks = "cluster_id"))

#split the multitype point process into several single type processes
#first, set the marks of the point process to be factors
marks(pp) <- factor(marks(pp))
ppls <- split(pp)

selection <- c('OD Mature', 'Ependymal', 'Microglia')
pp_sel <-  subset(pp, marks %in% selection, drop = TRUE)

#PRE: list of point pattern, corresponding celltypes of interest, functions to evaluate
#POST: result of the metric
metricResBoot <- function(ppls, celltype, fun){
  metric.res <- lohboot(ppls[[celltype]], fun = fun)
  metric.res$type <- celltype
  return(metric.res)
}
#PRE: celltypes, function to calculation and edge correction method
#POST: dataframe of 
metricResBootToDF <- function(celltype_ls, ppls, fun, edgecorr){
  res_ls <- lapply(celltype_ls, metricResBoot, fun = fun, ppls = ppls)
  #stick all values into a dataframe
  res_df <- c()
  for(i in 1:length(celltype_ls))res_df <- rbind(res_df, res_ls[[i]])
  return(res_df)
}

#PRE: Celltypes of interest, function to analyse, edge correction to perform
#POST: plot of the metric
plotMetric <- function(celltype_ls, ppls, fun, edgecorr){
  res_df <- metricResBootToDF(celltype_ls, ppls, fun, edgecorr)
  #plot the curve
  p <- ggplot(res_df, aes(x=r, y=res_df[[edgecorr]], color= type),size = 1, shape = 1)+
    geom_line(size = 2)+
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.25, show.legend = FALSE)+
    ggtitle(paste0(fun, '-function'))+
    geom_line(aes(x=r,y=theo, color = 'Poisson'), size = 2, linetype = "dashed")+
    ylab(edgecorr) +
    scale_color_manual(name='Point Processes',
                       breaks=c('Ependymal', 'Microglia', 'OD Mature', 'Poisson'),
                       values=c('Ependymal'='red', 'Microglia'='dark green', 'OD Mature'='blue', 'Poisson'='black'))+
    guides(fill = "none", color = guide_legend(override.aes = list(shape = 21)) )  +
    theme_light()+
    ylab('K(r): isotropic correction')
  return(p)
}

celltype_ls <- c("Ependymal", "OD Mature", "Microglia")
p_K <- plotMetric(celltype_ls, ppls, 'Kest', 'iso')

ggsave('outs/Kest.pdf', plot = p_K, width = 5, height = 3)

### SI figure pp

# homogeneous functions
res_ls <- lapply(list('Kest', 'Lest', 'pcf'), function(fun){
  res <- calcMetricPerFov(spe, 'OD Mature', subsetby = 'sample_id', fun = fun, marks = 'cluster_id', r_seq=NULL, by = c('Animal_ID','sample_id'))
  res <- subset(res, sample_id %in% c('-0.09', '0.01', '0.21'))
  return(res)
})

p_ls_homo <- lapply(res_ls, function(res){plotMetricPerFov(res, theo = TRUE, correction = "iso", x = "r", image_id = 'sample_id') + theme(legend.position = 'right') + labs(colour="Slice")})
p_homo <- wrap_plots(p_ls_homo)

# inhomogeneous functions
res_ls <- lapply(list('Kinhom', 'Linhom'), function(fun){
  res <- calcMetricPerFov(spe, 'OD Mature', subsetby = 'sample_id', fun = fun, marks = 'cluster_id', r_seq=NULL, by = c('Animal_ID','sample_id'), correction = 'isotropic') 
  res <- subset(res, sample_id %in% c('-0.09', '0.01', '0.21'))
  return(res)
})

res_pcf <- calcMetricPerFov(spe, 'OD Mature', subsetby = 'sample_id', fun = 'pcfinhom', marks = 'cluster_id', r_seq=NULL, by = c('Animal_ID','sample_id')) 
res_pcf <- subset(res_pcf, sample_id %in% c('-0.09', '0.01', '0.21'))

p_ls_inhomo <- lapply(res_ls, function(res){plotMetricPerFov(res, correction = "iso", theo = TRUE, x = "r", image_id = 'sample_id') + theme(legend.position = 'right') + labs(colour="Slice")})
p <- plotMetricPerFov(res_pcf, correction = "iso", theo = TRUE, x = "r", image_id = 'sample_id') + theme(legend.position = 'right') + labs(colour="Slice")

p_ls_inhomo <- wrap_plots(p_ls_inhomo)
p_inhomo <- wrap_plots(p_ls_inhomo, p, widths=c(2,1))

# locally scaled functions 

res_ls <- lapply(list('Kscaled', 'Lscaled'), function(fun){
  res <- calcMetricPerFov(spe, 'OD Mature', subsetby = 'sample_id', fun = fun, marks = 'cluster_id', r_seq=NULL, by = c('Animal_ID','sample_id'))
  res <- subset(res, sample_id %in% c('-0.09', '0.01', '0.21'))
  return(res)
})

p_ls_scaled <- lapply(res_ls, function(res){plotMetricPerFov(res, correction = "iso", theo = TRUE, x = "r", image_id = 'sample_id') + theme(legend.position = 'right') + labs(colour="Slice")})
p_scaled <- wrap_plots(p_ls_scaled)


p <- wrap_plots(list(p_homo, p_inhomo, p_scaled), nrow = 3, guides = 'collect') + plot_annotation(tag_levels = 'A')

ggsave('outs/pp_function_comparison.pdf', plot = p, width = 10, height = 10)
