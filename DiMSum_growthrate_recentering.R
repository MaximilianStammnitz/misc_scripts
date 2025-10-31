## 31.10.2025
## maximilian.stammnitz@crg.eu

###############################################
## Pre-processing of DiMSum growth rate data ##
###############################################


## 0. Environment ##
####################

## Libraries
packages <- c("viridis")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)

## Paths
setwd("") ## set your directory above the DiMSum folders


## 1. Input DiMSum tables ##
############################

## Input function
load_conc <- function(conc){
  load(paste0("conc", conc, "/conc_", conc, "_fitness_replicates.RData")) ## may need to adjust this path encoding to correctly load your files  
  data <- rbind(all_variants, synonymous)
  data <- data[-which(data$WT == T)[-1],] ## remove duplicated WT
  return(data)}

sPCA_data <- list("1" = load_conc("01"), "2" = load_conc("02"),
                  "3" = load_conc("03"), "4" = load_conc("04"),
                  "5" = load_conc("05"), "6" = load_conc("06"),
                  "7" = load_conc("07"), "8" = load_conc("08"),
                  "9" = load_conc("09"), "10" = load_conc("10"),
                  "11" = load_conc("11"), "12" = load_conc("12"),
                  "13" = load_conc("13"), "14" = load_conc("14"),
                  "15" = load_conc("15"), "16" = load_conc("16"),
                  "17" = load_conc("17"), "18" = load_conc("18"),
                  "19" = load_conc("19"), "20" = load_conc("20"),
                  "21" = load_conc("21"))

## Clean up environment
rm(packages, install_if_missing, load_blocks)


## 2. Normalise growth rates, using stops & synonymous WT variants ##
#####################################################################

## Error-weighted mean of growth rate estimates
sPCA_data <- lapply(sPCA_data, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
sPCA_data <- lapply(sPCA_data, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})

## Error-weighted mean of stop mutations
stops.sPCA <- lapply(sPCA_data, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## Error-weighted mean of WT/synonymous mutations
synWT.sPCA <- lapply(sPCA_data, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## Calculate coefficients for linear transformation
scaling <- cbind(data.frame(rbind(do.call(c, stops.sPCA), do.call(c, synWT.sPCA))))
rownames(scaling) <- c("stops", "syn")
colnames(scaling) <- 1:21

## Linear transformation (stop mode consistently to 0)
coefficients <- matrix(NA, ncol = 21, nrow = 2)
colnames(coefficients) <-  colnames(scaling)
for (i in 1:21){
  
  tmp.lm <- lm(formula = c(0, scaling[2,i]) ~ scaling[,i]) ## here you linearly rescale the growth rates between zero and the conc.-specific synonymous WT mode 
  tmp.lm <- summary(tmp.lm)
  coefficients[1,i] <- tmp.lm$coefficients[[2]]
  coefficients[2,i] <- tmp.lm$coefficients[[1]]  
  
}
for (i in 1:length(sPCA_data)){
  
  sPCA_data[[i]]$gr_normalised <- sPCA_data[[i]]$growthrate*coefficients[1,i] + coefficients[2,i] ## apply the coefficients to
  sPCA_data[[i]]$gr_sigma_normalised <- sPCA_data[[i]]$growthrate_sigma*coefficients[1,i]
  
}

## Clean up environment
rm(i, tmp.lm)


## 3. Visualise growth rates before vs. after normalisation ##
##############################################################

pdf("growth_rates_rescaled.pdf", height = 10, width = 25)

par(mfcol = c(1,2))
conc <- length(sPCA_data)

### original fitness
set.max <- c()
for (i in 1:conc){
  set.max <- c(set.max, max(density(sPCA_data[[i]]$growthrate)$y))
}
set.max <- max(set.max)

plot(density(sPCA_data[[1]]$growthrate), 
     bty = "n", 
     col = viridis(n = conc, alpha = 0.5)[1],
     lwd = 3,
     main = "", 
     cex.main = 3, 
     cex.axis = 1.8,
     cex.lab = 2.5, 
     xlab = "Growth rate (raw)", 
     xlim = c(-0.1, 0.5),
     ylim = c(0, set.max))

for(i in 2:conc){
  lines(density(sPCA_data[[i]]$growthrate), 
        col = viridis(n = conc, alpha = 0.5)[i],
        lwd = 3)
}

### legend
legend('topleft',
       legend = 1:conc,
       col = viridis(n = conc, alpha = 0.9),
       lwd = 3,
       cex = 1.5,
       bty = "n")

### fitness re-scaled
set.max <- c()
for (i in 1:conc){
  set.max <- c(set.max, max(density(sPCA_data[[i]]$gr_normalised)$y))
}
set.max <- max(set.max)

plot(density(sPCA_data[[1]]$gr_normalised), 
     bty = "n", 
     col = viridis(n = conc, alpha = 0.5)[1],
     lwd = 3,
     main = "", 
     cex.main = 3, 
     cex.axis = 1.8,
     cex.lab = 2.5, 
     xlab = "Growth rate (re-scaled)", 
     xlim = c(-0.1, 0.5),
     ylim = c(0, set.max))

for(i in 2:conc){
  lines(density(sPCA_data[[i]]$gr_normalised), 
        col = viridis(n = conc, alpha = 0.5)[i],
        lwd = 3)
}

dev.off()


## 4. Version ##
################

# sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Madrid
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1 tools_4.5.1
