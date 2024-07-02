source("code/utils.R")

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
  ggtitle("LISA curves of slice 0.01")

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
res <- functionalPCA(dat = dat, r = df_fdob$r |> unique(), knots = 30, pve = 0.99)
# plot a biplot
p_biplot <- plotFpca(dat = dat, res = res, colourby = 'sel')
p_biplot <- p_biplot + ggtitle("Biplot of LISA curves") +
  xlab('PC1') +
  ylab('PC2')
  


#assemble plots for paper
p_L <- plotMetric(plot_by = zstack_list, pp, celltype_ls, 'r', 'Lest', 'iso', bootstrap = TRUE)
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
      ggplot(aes(x = x, y = y, col = res$scores[,1])) +
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
