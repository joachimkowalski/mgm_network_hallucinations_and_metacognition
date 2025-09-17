#Which estimator to choose?
#from: Isvoranu, A. M., & Epskamp, S. (2023). Which estimation method to choose in network psychometrics? Deriving guidelines for applied researchers. Psychological methods, 28(4), 925.
library("shiny")
runGitHub("AdelaIsvoranu/simulation_exploration_app")
#simulation data that would reflect my dataset the closest -> MAGNA and DASS datasets, skewed data, n = 150, with and without transformations
#parameters -> specificity, sensitivity, precision, bridge sensitivity, precision and specificity (+top 25% of precision and sensitivity parameters), correlation and fade maximum false edge
#MGM with EBIC and 0.25 gamma is quite conservative - very specific and precise, however the sensitivity/discovery parameters are quite low. However, given my sample size, I can live with that.
#MGM CV with 10 folds is suboptimal than model selection with gamma 0 in almost every parameter. ebic glasso is quite sensitive but too un-specific and im-precise.
#evidently ranked transformation is necessary for skewed data

#load data
data <- read.csv("SONATINA_pre_networks_onlyclinical210125.csv")

#load packages -> at this point i'm not sure if all of those are necessary, but too anxious to remove them
library("rcompanion")
library("mice")
library("bootnet")
library("qgraph")
library("dplyr") 
library("mgm")

# List of variables to include
include <- c(
  "PANSS_positive",
  "PANSS_negative",
  "PANSS_cognitivedisorganization",
  "PANSS_hosility",
  "RHS1_total",
  "SCL1_socialph",
  "SCL1_vegetative",
  "SCL1_agoraph",
  "SCL1_depress",
  "IVIR_negative",
  "IVIR_positive",
  "MCQ1_positive",
  "MCQ1_uncontrollability",
  "MCQ1_cogconfidence",
  "MCQ1_lossofcontrol",
  "MCQ1_selfconscious",
  "I_CAS_strategies",
  "ACS1_focusing",
  "ACS1_shifting",
  "ACS1_divided",
  "DL_non_forced",
  "DL_forced_right",
  "DL_forced_left")

#Create the subset
subset <- data[, include]

#filtering based on DL_filter -> unreliable Dichotic Listening data is out; however, this is just one patient
subset[c("DL_non_forced", "DL_forced_right", "DL_forced_left")][data$DL_filter == 0, ] <- NA

#shorter names to be displayed on nodes
names(subset) <- c(
  "PANSS1",
  "PANSS2",
  "PANSS3",
  "PANSS4",
  "RHS",
  "SCL1",
  "SCL2",
  "SCL3",
  "SCL4",
  "IVI1",
  "IVI2",
  "MCQ1",
  "MCQ2",
  "MCQ3",
  "MCQ4",
  "MCQ5",
  "CAS",
  "ACS1",
  "ACS2",
  "ACS3",
  "DL_N",
  "DL_R",
  "DL_L")

#longer names for the plot legend
Etykiety <- scan(#"path to my hard drive",
                 what = "character", sep = "\n")


#imputation of missing data -> 11 sets of pmm averaged
layout(t(1:1))
md.pattern(subset)

imputation11 <- mice(subset, m = 11, defaultMethod = "pmm", seed = 1)
imputation11_subsets <- complete(imputation11, "all")
imputed_subset <- imputation11_subsets[[1]]
for (i in seq_along(imputation11_subsets)[-1]) {
  imputed_subset <- imputed_subset + imputation11_subsets[[i]]
}
imputed_subset <- imputed_subset / length(imputation11_subsets)


#ranking to account to skewed structure for MGM (see Isvoranu & Epskam, 2023)
imputed_ranked_subset <- imputed_subset
imputed_ranked_subset <- as.data.frame(lapply(imputed_subset, rank))


#model MGM
#Node prediction: https://jonashaslbeck.com/Predictability-in-network-models/
data_mgm <- imputed_ranked_subset
matrix_mgm <- as.matrix(data_mgm)
p <- ncol(matrix_mgm)

set.seed(1) 

#model specification
fit_obj <- mgm(data = matrix_mgm,
               type = rep('g', p),
               level = rep(1, p),
               lambdaSel = 'EBIC',
               lambdaGam = 0.25,
               ruleReg = 'OR', 
               verbatim = TRUE)

pred_obj <- predict(object = fit_obj, 
                    data = matrix_mgm,
                    errorCon = "R2")
#R^2
pred_obj$errors 

#network plot
layout(t(1:1))
qgraph(fit_obj$pairwise$wadj, # weighted adjacency matrix as input
       layout = 'spring', 
       pie = pred_obj$error[,2], # provide errors as input
       pieColor = rep('#377EB8',p),
       edge.color = fit_obj$pairwise$edgecolor, #why is it green, though? 
       labels = colnames(matrix_mgm),
       groups = list(symptoms = 1:9,
                     metacognitive = 10:17,
                     attention = 18:23),
       nodeNames = Etykiety,
       legend.cex = 0.22,
       theme = "colorblind")


#PANSS2 and DL_L are isolated -> re-analysis of the model with these two variables pruned
imputed_ranked_pruned_subset <- imputed_ranked_subset %>% select(-DL_L, -PANSS2)

#shorter list for the legend
Etykiety2 <- scan(#"location on hard drive",
                 what = "character", sep = "\n")

#MGM model estimation, second attempt 
data_mgm2 <- imputed_ranked_pruned_subset
matrix_mgm2 <- as.matrix(data_mgm2)
p2 <- ncol(matrix_mgm2)

set.seed(1) 

#model specification
fit_obj2 <- mgm(data = matrix_mgm2,
               type = rep('g', p2),
               level = rep(1, p2),
               lambdaSel = 'EBIC',
               lambdaGam = 0.25,
               ruleReg = 'OR', 
               verbatim = TRUE)

pred_obj2 <- predict(object = fit_obj2, 
                    data = matrix_mgm2,
                    errorCon = "R2")

#R^2 explained
pred_obj2$errors 

#colorblind friendly -> the mgm default is green :c pls, change it Jonas Haslbeck! 
fit_obj2$pairwise$edgecolor <- gsub("darkgreen", "blue", fit_obj2$pairwise$edgecolor)

#network plot
layout(t(1:1))
png("network_plot.png", width = 1700, height = 1100, res = 300)
mgm_graph <- qgraph(fit_obj2$pairwise$wadj, # weighted adjacency matrix as input
       layout = 'spring', 
       pie = pred_obj2$error[,2], # provide errors as input
       pieColor = rep('#377EB8',p2),
       posCol = "blue",
       edge.color = fit_obj2$pairwise$edgecolor,
       labels = colnames(matrix_mgm2),
       nodeNames = Etykiety2,
       groups = list(symptoms = 1:8,
                     "metacognitive factors" = 9:16,
                     "attention functioning" = 17:21),
       cut = 0,
       legend.cex = 0.22,
       theme = "colorblind")
dev.off()


#model stability via boostrap
res_obj2 <- resample(object = fit_obj2, data = matrix_mgm2, nB = 5000, nCores = 4)

#stability and centrality plots
png("res_plot2.png", width = 3000, height = 18000, res = 300)
stability_plot <- plotRes(object = res_obj2, axis.ticks = c(-.2, -.1, 0, .1, .2, .3, .4, .5), quantiles = c(0.05, 0.95), labels = colnames(matrix_mgm2))
dev.off()

png("centrality_plot.png",width = 1240, height = 840, res = 200)
centrality_plot <- centralityPlot(fit_obj2$pairwise$wadj, include = c("Strength", "Closeness"), scale = "z-scores", labels = colnames(matrix_mgm2), orderBy = "Strength")
dev.off()

#case-drop boostraping to estimate stability of network centrality estimates
bootstrap_model_mgm <- estimateNetwork(matrix_mgm2, default = "mgm", level = rep(1, p2), type = rep('g', p2), tuning = 0.25, verbose = TRUE, criterion = "EBIC", rule = "OR")
bootstrap_model_mgm2 <- bootnet(bootstrap_model_mgm, type = "case", nBoots = 5000, nCores = 4, statistics = c("strength", "closeness"))
png("centrality_stability_plot.png",width = 1240, height = 840, res = 200)
plot(bootstrap_model_mgm2, c("strength", "closeness"))
dev.off()

#extracting errors/R^2 to order them by strength and paste into centrality plots
centrality_plot_strength <- centrality_plot[["data"]] %>%
  filter(measure == "Strength") %>%
  select(node, measure, value)
centrality_plot_strength_ordered <- centrality_plot_strength %>%
  select(node, value) %>%
  arrange(desc(value))
rsquared_ordered <- merge(pred_obj2$errors, centrality_plot_strength_ordered, 
                          by.x = "Variable", by.y = "node", 
                          all.x = TRUE)
rsquared_ordered <- rsquared_ordered %>%
  arrange(desc(value))
rsquared_ordered <- rsquared_ordered[match(colnames(matrix_mgm2), rsquared_ordered$Variable), ] #have to reorder them again... 

#centrality plots with R^2
labels_combined <- paste(rsquared_ordered$Variable, sprintf("%.2f", rsquared_ordered$R2), sep = " - ")
labels_combined <- as.matrix(labels_combined)
png("centrality_plot_r2.png",width = 1240, height = 840, res = 200)
centrality_plot_r2 <- centralityPlot(fit_obj2$pairwise$wadj, include = c("Strength", "Closeness"), labels = labels_combined, scale = "z-scores", orderBy = "Strength", theme_bw = TRUE) 
dev.off()

#extracting values of edges to describe several strongest ones in the paper
wadj_matrix <- fit_obj2$pairwise$wadj
wadj_values <- as.vector(wadj_matrix)
top_indices <- order(wadj_values, decreasing = TRUE)[1:80] #they are pairwise, so 40x2
row_col_indices <- arrayInd(top_indices, dim(wadj_matrix))
variable_names <- colnames(matrix_mgm2)
top_values <- data.frame(
  Variable_A = variable_names[row_col_indices[, 1]],
  Variable_B = variable_names[row_col_indices[, 2]],
  Value = wadj_values[top_indices])
print(top_values)

#saving the wadj matrix for supplement
wadj_matrix <- fit_obj2$pairwise$wadj
colnames(wadj_matrix) <- rownames(wadj_matrix) <- colnames(matrix_mgm2)
write.csv(wadj_matrix, file = "wadj_matrix.csv", row.names = TRUE)

#mean weight from wadj but excluding zeroes
mean_nonzero <- mean(fit_obj2$pairwise$wadj[fit_obj2$pairwise$wadj != 0])
print(mean_nonzero)

#average predictability
mean_r2 <- mean(pred_obj2$errors$R2)
sd_r2 <- sd(pred_obj2$errors$R2)
print(mean_r2)
print(sd_r2)
