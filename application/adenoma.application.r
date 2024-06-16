library("readxl")
library(plyr)
library(mgcv); library(MASS);
library(fql);
library(betareg)
library(RColorBrewer)
library(GUniFrac)
library(boot)
library(sandwich)
library(lmtest)
library(robust)
library(robustbase)
library(lmtest)
library('RSimNorm')
library('MicrobiomeStat')
library(phyloseq)
library(ANCOMBC)

load("adenomas.RData")

class(data.obj$meta.dat)
dim(data.obj$meta.dat)
head(data.obj$meta.dat)
data.obj$meta.dat$Extraction_treatment_Group
data.obj$meta.dat$Diagnosis
data.obj$meta.dat$Diagnosis2
table(data.obj$meta.dat$Diagnosis)
table(data.obj$meta.dat$Diagnosis2)
table(data.obj$meta.dat$gender)

class(data.obj$otu.tab)
head(data.obj$otu.tab)
dim(data.obj$otu.name)

lapply(data.obj$abund.list, dim)
lapply(data.obj$abund.list, rownames)
class(data.obj$abund.list$Species)
head(data.obj$abund.list$Species)
identical(as.numeric(data.obj$otu.tab),as.numeric(data.obj$abund.list$Species))

data.obj$tree

data.obj$otu.name.full

data.obj$size.factor


sum(data.obj$abund.list$Genus['Proteobacteria;Eikenella',]==0)
length(data.obj$abund.list$Genus['Proteobacteria;Eikenella',])

#genus with prevelance over 10%
otu_prevalence <- apply(data.obj$abund.list$Genus, 1, function(x) sum(x > 0) / length(x))
otu_prevalence
otu_filtered <- data.obj$abund.list$Genus[otu_prevalence >= 0.05, ]
dim(otu_filtered)
countdata=otu_filtered
dim(countdata)


#WINSORLIZATION
winsor<-function(vector,quantile){
  q=quantile(vector,quantile)
  vector[vector>q]=round(q)
  return(vector)
}

glm_coef <- function(data, indices) {
  fit <- glm(formula, family = poisson(link="log"), data = data[indices,])
  return(coef(fit))
}


tc=log(colSums(countdata))
length(tc)


#noraml=0,abnormal=1
diag=as.numeric(data.obj$meta.dat$Diagnosis2=='abnormal')


batch1=as.numeric(data.obj$meta.dat$SEQ_pool==1)
batch2=as.numeric(data.obj$meta.dat$SEQ_pool==2)



smoke=data.obj$meta.dat$smoke
gender=as.numeric(data.obj$meta.dat$gender)-1


polyps=data.obj$meta.dat$POLYPS-1
NUMPLYP_risk=data.obj$meta.dat$NUMPLYP_risk


ageg=data.obj$meta.dat$ageg
agegrp=data.obj$meta.dat$agegrp-1

final_data = cbind(diag,smoke,gender,batch1,batch2,tc)
meta.dat = final_data
dim(meta.dat)
countdata = countdata[,complete.cases(meta.dat)]
meta.dat = meta.dat[complete.cases(meta.dat),]
final_data = final_data[complete.cases(final_data),]


formula=y~diag+gender+smoke+batch1+batch2+offset(tc)
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
    
    y <- winsor(countdata[III, ], 1)
    datam <- cbind(y = y, final_data, logtc = log(colSums(countdata)))
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
    
    #fql
    tryCatch({
      fitg4<- fql(formula,datas)
      est4[III,1:9]=as.numeric(unlist(fitg4[[1]][2:4,]))
    }, error = function(e) {
      message("Error in FQL model at iteration ", III, ": ", e$message)
    })

    # Poisson Sandwich
    tryCatch({
      fitg5 <- glm(formula, family = poisson(link = 'log'), data = datas)
      est5[III, 1:9] <- as.numeric(coeftest(fitg5, vcov. = vcovHC(fitg5, type = "HC3"))[2:4, c(1, 2, 4)])
    }, error = function(e) {
      message("Error in Poisson Sandwich model at iteration ", III, ": ", e$message)
    })
    
    # Quasi-Poisson
    tryCatch({
      fitg6 <- glm(formula, family = quasipoisson(link = 'log'), data = datas)
      est6[III, 1:9] <- as.numeric(summary(fitg6)$coefficients[2:4, -3])
    }, error = function(e) {
      message("Error in Quasi-Poisson model at iteration ", III, ": ", e$message)
    })
  }
)

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


whichsignificant(est,varname='p_beta1')
whichsignificant(est2)
whichsignificant(est3)
whichsignificant(est4)
whichsignificant(est5)
whichsignificant(est6)



whichsignificant.adj(est)
whichsignificant.adj(est2)
whichsignificant.adj(est3)
whichsignificant.adj(est4)
whichsignificant.adj(est5)
whichsignificant.adj(est6)

whichsignificant(est,varname='p_beta3')
whichsignificant(est2,varname='p_beta3')
whichsignificant(est3,varname='p_beta3')
whichsignificant(est4,varname='p_beta3')
whichsignificant(est5,varname='p_beta3')
whichsignificant(est6,varname='p_beta3')


whichsignificant.adj(est,varname='p_beta3')
whichsignificant.adj(est2,varname='p_beta3')
whichsignificant.adj(est3,varname='p_beta3')
whichsignificant.adj(est4,varname='p_beta3')
whichsignificant.adj(est5,varname='p_beta3')
whichsignificant.adj(est6,varname='p_beta3')

feature.dat = countdata
meta.dat = final_data
meta.dat = data.frame(meta.dat)
rownames(meta.dat) = colnames(feature.dat) = paste0('X',1:dim(meta.dat)[1])
# nb poi boot fql sandwich qp
linda.obj  <- linda(feature.dat = feature.dat,meta.dat = meta.dat,formula = '~ diag + gender + smoke + batch1 + batch2',
                    feature.dat.type = 'count',zero.handling = c('pseudo-count'),alpha=0.05,is.winsor = T,outlier.pct = 0.00)

which(linda.obj$output$diag$pvalue<0.05)
which(linda.obj$output$diag$padj<0.05)

which(linda.obj$output$smoke$pvalue<0.05)
which(linda.obj$output$smoke$reject)


zico.obj.diag <- ZicoSeq(meta.dat = data.frame(meta.dat), feature.dat = as.matrix(feature.dat), grp.name = 'diag',  adj.name = c('gender' ,'smoke' ,'batch1' ,'batch2'),
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

which(zico.obj$p.raw <0.05)
which(zico.obj$p.adj.fdr <0.05)

zico.obj <- ZicoSeq(meta.dat = data.frame(meta.dat), feature.dat = as.matrix(feature.dat), grp.name = 'smoke',  adj.name = c('gender' ,'diag' ,'batch1' ,'batch2'),
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

which(zico.obj$p.raw <0.05)
which(zico.obj$p.adj.fdr <0.05)
# Creating a phyloseq object

assays = S4Vectors::SimpleList(counts = as.matrix(feature.dat))
smd = S4Vectors::DataFrame(meta.dat)
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
#head(meta.dat)
# Running ANCOMBC2
ancombc2.obj <- ancombc2(
  data = tse,tax_level = 'Species',
  #assay.type = "counts",
  fix_formula = "diag + gender + smoke + batch1 + batch2",  p_adj_method = "fdr",
  alpha = 0.05,
  neg_lb = FALSE,
  prv_cut = 0
)

which(ancombc2.obj$res$p_diag<0.05)
which(ancombc2.obj$res$q_diag<0.05)
which(ancombc2.obj$res$q_smoke<0.05)


heatmapvalues = data.frame(
  NegativeBinomial = p.adjust(est[,'p_beta1'],method = 'BH'),
  Poisson = p.adjust(est2[,'p_beta1'],method = 'BH'),
  #FQL = p.adjust(est4[,'p_beta1'],method = 'BH'),
  Bootstrap = p.adjust(est3[,'p_beta1'],method = 'BH'),
  Sandwich = p.adjust(est5[,'p_beta1'],method = 'BH'),
  #QuasiPoisson =p.adjust(est6[,'p_beta1'],method = 'BH'),
  Linda = linda.obj$output$diag$padj,
  Zicoseq = zico.obj.diag$p.adj.fdr,
  Ancombc2 = ancombc2.obj$res$q_diag
)
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
  labs(x = "Methods", y = "Virus", fill = "adjust P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_blank())
heatmap_plot

ggsave(filename = "adenomas_diagnosis0614.jpeg", plot = heatmap_plot, width = 6.5, height = 8.5, dpi = 300)

# diag
sets <- list(
  NegativeBinomial = whichsignificant.adj(est,varname = 'p_beta1'),
  Poisson = whichsignificant.adj(est2,varname = 'p_beta1'),
  Bootstrap = whichsignificant.adj(est3,varname = 'p_beta1'),
  Sandwich = whichsignificant.adj(est5,varname = 'p_beta1'),
  Linda = which(linda.obj$output$diag$reject),
  Zicoseq = which(zico.obj.diag$p.adj.fdr<0.05),
  Ancombc2 = which(ancombc2.obj$res$q_diag<0.05)
)
sets
# Fit the Euler diagram
fit <- euler(sets,shape = "ellipse")
fit <- euler(sets)
plot(fit, quantities = TRUE)
# Plot
plot_obj <- plot(fit, quantities = TRUE)
plot_obj
# Save the plot using ggsave
ggsave("adenomasdiageuler.png", plot = plot_obj, width = 10, height = 8, dpi = 300)

library(UpSetR)
# Find the union of all sets
all_elements <- unique(unlist(sets))

# Create an empty binary matrix
binary_matrix <- matrix(0, nrow = length(all_elements), ncol = length(sets))
rownames(binary_matrix) <- all_elements
colnames(binary_matrix) <- names(sets)

# Fill the binary matrix
for (i in 1:length(sets)) {
  binary_matrix[rownames(binary_matrix) %in% sets[[i]], i] <- 1
}

# Convert to data frame
binary_df <- as.data.frame(binary_matrix)

# Generate UpSet plot
upset(binary_df, sets = colnames(binary_df), keep.order = TRUE)
# Customize the colors
upset(
  binary_df, 
  sets = colnames(binary_df), 
  keep.order = TRUE,
  main.bar.color = "darkblue",  # Color of the main bar plot
  matrix.color = "darkgreen",   # Color of the matrix elements
  sets.bar.color = "darkred",   # Color of the set size bars
  text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.5),  # Scale text sizes
  mainbar.y.label = "Intersection Size",    # Label for the main bar plot
  sets.x.label = "Set Size"                  # Label for the set size bars
)

