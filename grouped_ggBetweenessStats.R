grouped_ggBetweenessStats <- function(data, title_Text = NULL){
  require(ggstatsplot)
  require(viridis)
  ggstatsplot::grouped_ggbetweenstats(
    data = data,
    x = algorithm,
    y = scores,
    grouping.var = Metric,
    pairwise.comparisons = F, # display significant pairwise comparisons
    # pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
    # p.adjust.method = "bonferroni", # method for adjusting p-values for multiple comparisons
    # conf.level = 0.99, # changing confidence level to 99%
    ggplot.component = list( # adding new components to `ggstatsplot` default
      ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())
      ),
    k = 3,
    # title.prefix = "Comparison",
    caption = substitute(paste(
      italic("MCMC methods"),
      " vs non-MCMC methods"
    )),
    package = "RColorBrewer",
    palette = "Spectral",
    messages = FALSE,
    nrow = 2,
    # title.text = "MCMC vs non-MCMC methods for Cancer nepolitan network"
    title.text = title_Text
  )
}

library(readxl)
library(dplyr)

res <- read_excel("data/gold standard data analysis.xlsx", sheet = 2) %>% as.data.frame()

# sn.data.1 <- data.frame(score = res[2:16,2], algorithm = "mh")
# sn.data.2 <- data.frame(score = res[2:16,5], algorithm = "ns")
# sn.data.3 <- data.frame(score = res[2:16,8], algorithm = "har")
# sn.data.4 <- data.frame(score = res[2:16,11], algorithm = "gs")
# sn.data.5 <- data.frame(score = res[2:16,14], algorithm = "iamb")
# sn.data.6 <- data.frame(score = res[2:16,17], algorithm = "inter.iamb")
# sn.data.7 <- data.frame(score = res[2:16,20], algorithm = "fast.iamb")
# sn.data.8 <- data.frame(score = res[2:16,23], algorithm = "rsmax2")
# sn.data.9 <- data.frame(score = res[2:16,26], algorithm = "mmpc")
# sn.data.10 <- data.frame(score = res[2:16,29], algorithm = "hc")
# sn.data.11 <- data.frame(score = res[2:16,32], algorithm = "tabu")
# sn.data.12 <- data.frame(score = res[2:16,35], algorithm = "mmhc")
# sn.data.13 <- data.frame(score = res[2:16,38], algorithm = "h2pc")
sn.data <- base::rbind(
  data.frame(V1 = res[2:16,2], V2 ="mh", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,5], V2 = "ns", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,8], V2 = "har", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,11], V2 = "gs", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,14], V2 = "iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,17], V2 = "inter.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,20], V2 = "fast.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,23], V2 = "rsmax2", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,26], V2 = "hc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,29], V2 = "tabu", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,32], V2 = "mmhc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,35], V2 = "h2pc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2)
)
sn.data = cbind(sn.data, Metric = "Sensitivity")
sn.data$scores = as.numeric(sn.data$scores)
grouped_ggBetweenessStats(sn.data, 
                          title_Text = "MCMC vs non-MCMC methods for Cancer Neapolitan network"
                          )
# ggplot(sn.data, aes(x=algorithm, y = scores)) + 
#   geom_violin(trim = F) +
#   stat_summary(fun = median, geom="point", size=2, color="red") +
#   geom_boxplot(width=0.1)
# 
# # ggplot(ToothGrowth, aes(x=dose, y=len)) + 
#   geom_violin(trim=T, fill="gray")+
#   labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")+
#   geom_boxplot(width=0.1)+
#   theme_classic()
# Change color by groups
# dp <- ggplot(sn.data, aes(x=algorithm, y=scores, fill=algorithm)) + 
#   geom_violin(trim=FALSE) +
#   geom_boxplot(width=0.1, fill="white")+
#   labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")
# dp + theme_classic()
# # Continusous colors
# dp + scale_fill_brewer(palette="Blues") + theme_classic()
# # Discrete colors
# dp + scale_fill_brewer(palette="Dark2") + theme_minimal()
# # Gradient colors
# dp + scale_fill_brewer(palette="RdBu") + theme_minimal()
grouped_ggBetweenessStats(sn.data)


# sp.data.1 <- data.frame(score = res[2:16,3], algorithm = "mh")
# sp.data.2 <- data.frame(score = res[2:16,6], algorithm = "ns")
# sp.data.3 <- data.frame(score = res[2:16,9], algorithm = "har")
# sp.data.4 <- data.frame(score = res[2:16,12], algorithm = "gs")
# sp.data.5 <- data.frame(score = res[2:16,15], algorithm = "iamb")
# sp.data.6 <- data.frame(score = res[2:16,18], algorithm = "inter.iamb")
# sp.data.7 <- data.frame(score = res[2:16,21], algorithm = "fast.iamb")
# sp.data.8 <- data.frame(score = res[2:16,24], algorithm = "rsmax2")
# sp.data.9 <- data.frame(score = res[2:16,27], algorithm = "mmpc")
# sp.data.10 <- data.frame(score = res[2:16,30], algorithm = "hc")
# sp.data.11 <- data.frame(score = res[2:16,33], algorithm = "tabu")
# sp.data.12 <- data.frame(score = res[2:16,36], algorithm = "mmhc")
# sp.data.13 <- data.frame(score = res[2:16,39], algorithm = "h2pc")
sp.data = bind_rows(
  data.frame(V1 = res[2:16,3], V2 = "mh", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,6], V2 = "ns", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,9], V2 = "har", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,12], V2 = "gs", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,15], V2 = "iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,18], V2 = "inter.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,21], V2 = "fast.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,24], V2 = "rsmax2", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,27], V2 = "hc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,30], V2 = "tabu", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,33], V2 = "mmhc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,36], V2 = "h2pc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2)
  )
sp.data = cbind(sp.data, Metric = "Specificity")
sp.data$scores = as.numeric(sp.data$scores)
grouped_ggBetweenessStats(sp.data)
# f1.data.1 <- data.frame(score = res[2:16,4], algorithm = "mh")
# f1.data.2 <- data.frame(score = res[2:16,7], algorithm = "ns")
# f1.data.3 <- data.frame(score = res[2:16,10], algorithm = "har")
# f1.data.4 <- data.frame(score = res[2:16,13], algorithm = "gs")
# f1.data.5 <- data.frame(score = res[2:16,16], algorithm = "iamb")
# f1.data.6 <- data.frame(score = res[2:16,19], algorithm = "inter.iamb")
# f1.data.7 <- data.frame(score = res[2:16,22], algorithm = "fast.iamb")
# f1.data.8 <- data.frame(score = res[2:16,25], algorithm = "rsmax2")
# f1.data.9 <- data.frame(score = res[2:16,28], algorithm = "mmpc")
# f1.data.10 <- data.frame(score = res[2:16,31], algorithm = "hc")
# f1.data.11 <- data.frame(score = res[2:16,34], algorithm = "tabu")
# f1.data.12 <- data.frame(score = res[2:16,37], algorithm = "mmhc")
# f1.data.13 <- data.frame(score = res[2:16,40], algorithm = "h2pc")
f1.data <- bind_rows(
  data.frame(V1 = res[2:16,4], V2 = "mh", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,7], V2 = "ns", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,10], V2 = "har", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,13], V2 = "gs", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,16], V2 = "iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,19], V2 = "inter.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,22], V2 = "fast.iamb", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,25], V2 = "rsmax2", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,28], V2 = "hc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,31], V2 = "tabu", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,34], V2 = "mmhc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2),
  data.frame(V1 = res[2:16,37], V2 = "h2pc", stringsAsFactors = F) %>% dplyr::rename(scores = V1, algorithm = V2)
)
f1.data = cbind(f1.data, Metric = "F1_score")
f1.data$scores = as.numeric(f1.data$scores)
grouped_ggBetweenessStats(f1.data)

# plt.data = base::rbind(sn.data, sp.data, f1.data)
grouped_ggBetweenessStats(plt.data)
