
# Load necessary libraries
library(dplyr)
library(readxl)
library(tidyr)
library(dplyr)
library(mgcv); 
library(VGAM)
library(MASS)
library(rmutil)
library(fql)
library(actuar)
library(boot)
library(sandwich)
library(lmtest)
library(robust)
library(robustbase)
library(openxlsx)
library('RSimNorm')
library('MicrobiomeStat')
library(phyloseq)
library(ANCOMBC)
library(GUniFrac)
library(eulerr)

# Read the first file
data1 <- read.delim("v35.genus.matrix.upload.yz_trimesters.txt", header = TRUE, sep = "\t")

# Process GA_at_collection by adding weeks and days
data1 <- data1 %>%
  mutate(GA_at_collection = GA_at_collection__weeks + GA_at_collection__days / 7)


# Read the second file
data2 <- read.delim("Race spreadsheet.txt", header = TRUE, sep = "\t")

head(data1)
colnames(data1)
dim(data1)

head(data2)
dim(data2)

merged_data <- merge(data1, data2, by = "Individual_Name", all.x = TRUE)
merged_data$MGI.name
dim(merged_data)
colnames(merged_data)
final_data <- merged_data[,-c(7,253)]
final_data
colnames(final_data)
table(final_data$Race)
final_data$raceblack= as.numeric(final_data$Race=='Black')
final_data$racewhite = as.numeric(final_data$Race=='White')
final_data$total_taxa_abundance = rowSums(final_data[,7:250])

final_data_first_observation_2nd_trimester <- final_data %>%
  filter(GA_at_collection >= 13, GA_at_collection < 28) %>%
  group_by(Individual_Name) %>%
  filter(row_number() == 1)
dim(final_data_first_observation_2nd_trimester)
hist(final_data_first_observation_2nd_trimester$GA_at_collection)

glm_coef <- function(data, indices) {
  fit <- glm(formula, family = poisson(link="log"), data = data[indices,])
  return(coef(fit))
}
winsor<-function(vector,quantile){
  q=quantile(vector,quantile)
  vector[vector>q]=round(q)
  return(vector)
}
whichsignificant = function(est,varname='p_beta1'){
  est=data.frame(est)
  est[,varname]
  return(which(est[,varname]<0.05))
}
whichsignificant.adj = function(est,varname='p_beta1'){
  est=data.frame(est)
  est[,varname]
  return(which(p.adjust(est[,varname],method = 'fdr')<0.05))
}

##########################################################################
analysis_data=final_data_first_observation_2nd_trimester

countdata = t(analysis_data[,7:250])
dim(countdata)
#count data: row taxa; col sample
sum(rowMeans(countdata != 0)>0.1)
colSums(countdata > 0) / nrow(countdata)
colMeans(countdata!= 0)

countdata = countdata[rowMeans(countdata != 0)>0.1,]

dim(countdata)
dim(analysis_data)
formula = y~ raceblack + preterm_37 +    GA_at_collection + offset(logtc)


iters = dim(countdata)[1]
est =matrix(999, iters, 9);
colnames(est)<-c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')

est2=matrix(999, iters, 9);
colnames(est2)<-c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')

est3=matrix(999, iters, 9);
colnames(est3)<-c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')

est4=matrix(999, iters, 9);
colnames(est4)<-c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')

est5=matrix(999, iters, 9);
colnames(est5) = c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')

est6=matrix(999, iters, 9);
colnames(est6)<-c('beta1','beta2','beta3','se_beta1','se_beta2','se_beta3','p_beta1','p_beta2','p_beta3')
suppressWarnings(
  for (III in 1:iters) {
    print(III)
    set.seed(III)
    y <- winsor(countdata[III, ], 1)
    datam <- cbind(y = y, analysis_data, logtc = log(colSums(countdata)))
    colnames(datam)[1] <- 'y'
    datas <- data.frame(datam)
    
    # Negative Binomial
    tryCatch({
      fitg <- glm.nb(formula, link = "log", data = datas)
      est[III, 1:9] <- as.numeric(summary(fitg)$coefficients[2:4, -3])
    }, error = function(e) {
      message("Error in Negative Binomial model at iteration ", III, ": ", e$message)
    })
    
    # Poisson
    tryCatch({
      fitg2 <- glm(formula, poisson(link = "log"), data = datas)
      est2[III, 1:9] <- as.numeric(summary(fitg2)$coefficients[2:4, -3])
    }, error = function(e) {
      message("Error in Poisson model at iteration ", III, ": ", e$message)
    })
    
    # Poisson Bootstrap
    tryCatch({
      fitg3coef <- glm_coef(data = datas)
      fitg3results <- boot(data = datas, statistic = glm_coef, R = 1000)
      fitg3boot_se <- apply(fitg3results$t[complete.cases(fitg3results$t),], 2, sd)
      Zvalue <- fitg3coef / fitg3boot_se
      Pvalue <- 2 * (1 - pnorm(abs(Zvalue)))
      est3[III, 1:9] <- c(fitg3coef[2:4], fitg3boot_se[2:4], Pvalue[2:4])
    }, error = function(e) {
      message("Error in Poisson Bootstrap model at iteration ", III, ": ", e$message)
    })
    
    # Poisson Sandwich
    tryCatch({
      fitg5 <- glm(formula, family = poisson(link = 'log'), data = datas)
      est5[III, 1:9] <- as.numeric(coeftest(fitg5, vcov. = vcovHC(fitg5, type = "HC3"))[2:4, c(1, 2, 4)])
    }, error = function(e) {
      message("Error in Poisson Sandwich model at iteration ", III, ": ", e$message)
    })
    

  }
)


feature.dat = countdata
meta.dat = analysis_data
meta.dat = data.frame(meta.dat)
rownames(meta.dat) = colnames(feature.dat) = paste0('X',1:dim(meta.dat)[1])
# nb poi boot fql sandwich qp
linda.obj  <- linda(feature.dat = feature.dat,meta.dat = meta.dat,formula = '~ preterm_37 + GA_at_collection + raceblack',
                    feature.dat.type = 'count',zero.handling = c('pseudo-count'),alpha=0.05,is.winsor = T,outlier.pct = 0.00)



zico.obj.preterm_37 <- ZicoSeq(meta.dat = data.frame(meta.dat), feature.dat = as.matrix(feature.dat), grp.name = 'preterm_37',  adj.name = c('GA_at_collection','raceblack'),
                    feature.dat.type = "count", 
                    # Filter to remove rare taxa 
                    prev.filter = 0.05, mean.abund.filter = 0, max.abund.filter = 0, min.prop = 0, 
                    # Winsorization to replace outliers 
                    is.winsor = FALSE,# outlier.pct = 0.03, winsor.end = 'top',
                    # Posterior sampling to impute zeros 
                    is.post.sample = TRUE, post.sample.no = 25, 
                    # Multiple link functions to capture diverse taxon-covariate relation 
                    link.func = list( function(x) log(x+10e-10)),
                    stats.combine.func = max, 
                    # Permutation-based multiple testing correction 
                    perm.no = 100, strata = NULL, 
                    # Reference-based multiple stage normalization 
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2, 
                    # Family-wise error rate control 
                    is.fwer = FALSE, verbose = TRUE, return.feature.dat = FALSE) 



# Creating a phyloseq object

assays = S4Vectors::SimpleList(counts = as.matrix(feature.dat))
smd = S4Vectors::DataFrame(meta.dat)
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
#head(meta.dat)
# Running ANCOMBC2
ancombc2.obj <- ancombc2(
  data = tse,#tax_level = 'Species',
  #assay.type = "counts",
  fix_formula = "preterm_37 + GA_at_collection__weeks + raceblack",  p_adj_method = "fdr",
  alpha = 0.05,
  neg_lb = FALSE,
  prv_cut = 0
)



 # pretermbirth
sets <- list(
  NegativeBinomial = whichsignificant.adj(est,varname = 'p_beta2'),
  Poisson = whichsignificant.adj(est2,varname = 'p_beta2'),
  Bootstrap = whichsignificant.adj(est3,varname = 'p_beta2'),
  Sandwich = whichsignificant.adj(est5,varname = 'p_beta2'),
  linda = which(linda.obj$output$preterm_37$reject),
  ZicoSeq = which(zico.obj.preterm_37$p.adj.fdr<0.05),
  ancombc2 = which(ancombc2.obj$res$q_preterm_37<0.05)
)
sets
# Fit the Euler diagram
fit <- euler(sets,shape = "ellipse")
plot(fit, quantities = TRUE)
# Plot
plot_obj <- plot(fit, quantities = TRUE)
# Save the plot using ggsave
ggsave("euler_diagramPTBseed2023.0605.png", plot = plot_obj, width = 10, height = 8, dpi = 300)


#preterm birth
heatmapvalues = data.frame(
  NB = p.adjust(est[,'p_beta2'],method = 'BH'),
  POI = p.adjust(est2[,'p_beta2'],method = 'BH'),
  Bootstrap = p.adjust(est3[,'p_beta2'],method = 'BH'),
  Sandwich = p.adjust(est5[,'p_beta2'],method = 'BH'),
  Linda = linda.obj$output$preterm_37$padj,
  Zicoseq = zico.obj.preterm_37$p.adj.fdr,
  Ancombc2 = ancombc2.obj$res$q_preterm_37
)

# Replace NaN values with 1
heatmapvalues[is.na(heatmapvalues)] <- 1

heatmapvalues
rownames(heatmapvalues)
library(reshape2)
# Reshape data for ggplot2
melted_heatmapvalues <- melt(as.matrix(heatmapvalues))

# Rename columns for clarity
colnames(melted_heatmapvalues) <- c("Virus", "Measure", "Value")

library(ggplot2)
# Create the heatmap
cutoff=0.05
# Create the heatmap with custom colors and border
# Create the heatmap with custom colors and border
heatmap_plot = ggplot(melted_heatmapvalues, aes(x=Measure, y=Virus, fill=Value)) +
  geom_tile() +
  scale_fill_gradient2(low =  "#5bb963", mid = "#d9efdb", high ="white", midpoint = 0.5, 
                       limits = c(min(melted_heatmapvalues$Value), max(melted_heatmapvalues$Value))) +
  geom_tile(data = subset(melted_heatmapvalues, Value < cutoff), color = "black", fill = NA, linewidth = 0.5) +
  labs(x = "Methods", y = "taxa", fill = "adjust P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_blank())
heatmap_plot

ggsave(filename = "preterm_2ndtrimester.0614.jpeg", plot = heatmap_plot, width = 6.5, height = 8.5, dpi = 300)
