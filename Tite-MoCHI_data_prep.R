## Standardised data preparation of Tite-MoCHI runs
## 13.10.2025
## maximilian.stammnitz@crg.eu

## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "viridis")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input data: DiMSum files from a dose-response DMS experiment ##
#####################################################################

## DiMSum input function
load_blocks <- function(conc){
  load(paste0(conc, "/", conc, "_fitness_replicates.RData"))
  data <- rbind(all_variants, synonymous)
  data <- data[-which(data$WT == T)[-1],] ## remove duplicated WT
  return(data)}

## Adjust: input list of the number of concentrations used
input.data <- list("1" = load_blocks("01"), "2" = load_blocks("02"), "3" = load_blocks("03"), 
                   "4" = load_blocks("04"), "5" = load_blocks("05"), "6" = load_blocks("06"),
                   "7" = load_blocks("07"), "8" = load_blocks("08"), "9" = load_blocks("09"), 
                   "10" = load_blocks("10"), "11" = load_blocks("11"), "12" = load_blocks("12"), 
                   "13" = load_blocks("13"), "14" = load_blocks("14"), "15" = load_blocks("15"))

## Adjust: label data with the actual concentrations used
dosages <- rep(NA, length(input.data))
dosages[1] <- 100000
dosages[length(input.data)] <- 0
for(i in 2:c(length(input.data) - 1)){
  dosages[i] <- 100000 * (1/2)^c(i-1)
}
names(input.data) <- dosages


## 2. Re-scale fitness to -1 (stops) and 0 (WT/syn) ##
######################################################

## Calculate the error-weighted mean of stops and WTs
input.data <- lapply(input.data, function(x){x$fitness_over_sigmasquared <- x$fitness/(x$sigma)**2; return(x)})
input.data <- lapply(input.data, function(x){x$fitness_one_over_sigmasquared <- 1/(x$sigma)**2; return(x)})

## STOPs
input.stops <- lapply(input.data, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$fitness_over_sigmasquared, na.rm = TRUE) / sum(x$fitness_one_over_sigmasquared, na.rm = TRUE); return(x)})

## WTs/synonymous
input.WTs <- lapply(input.data, function(x){x <- x[which(x$Nham_aa == 0),,drop = F]; x <- sum(x$fitness_over_sigmasquared, na.rm = TRUE)/sum(x$fitness_one_over_sigmasquared, na.rm = TRUE); return(x)})

## Calculate coefficients for re-scaling to -1 / 0
scaling_data <- cbind(data.frame(rbind(do.call(c, input.stops), do.call(c, input.WTs))))
colnames(scaling_data) <- 1:length(input.data)
rownames(scaling_data) <- c("STOPs", "WTs")
coefficients <- matrix(NA, ncol = ncol(scaling_data), nrow = 2)
colnames(coefficients) <- colnames(scaling_data)
rownames(coefficients) <- c("Slope", "Intercept")
for (i in 1:ncol(scaling_data)){
  
  tmp.lm <- lm(formula = c(-1,0) ~ scaling_data[,i])
  tmp.lm <- summary(tmp.lm)
  coefficients["Slope",i] <- tmp.lm$coefficients[[2]]
  coefficients["Intercept",i] <- tmp.lm$coefficients[[1]]  
  
}

## Apply to calculate re-scaled fitnesses and fitness errors
for (i in 1:length(input.data)){
  
  input.data[[i]]$fitness_rescaled <- input.data[[i]]$fitness*coefficients["Slope",i] + coefficients["Intercept",i]
  input.data[[i]]$sigma_rescaled <- input.data[[i]]$sigma*coefficients["Slope",i]
  
}

## sanity check
pdf("re-scaled_dose-response_fitness.pdf", height = 10, width = 25)

par(mfcol = c(1,2))

### original fitness
set.max <- c()
for (i in 1:length(input.data)){
  set.max <- c(set.max, max(density(input.data[[i]]$fitness)$y))
}
set.max <- max(set.max)

plot(density(input.data[[1]]$fitness), 
     bty = "n", 
     col = viridis(n = length(input.data), alpha = 0.5)[1],
     lwd = 3,
     main = "", 
     cex.main = 3, 
     cex.axis = 1.8,
     cex.lab = 2.5, 
     xlab = "Fitness", 
     xlim = c(-1.5, 0.5),
     ylim = c(0, set.max))

for(i in 2:length(input.data)){
  lines(density(input.data[[i]]$fitness), 
        col = viridis(n = length(input.data), alpha = 0.5)[i],
        lwd = 3)
}

### legend
legend('topleft',
       legend = c(paste0(round(as.numeric(names(input.data)[1:length(input.data)]),3), " nM")),
       col = viridis(n = length(input.data), alpha = 0.9),
       lwd = 3,
       cex = 1.5,
       bty = "n")

### fitness re-scaled
set.max <- c()
for (i in 1:length(input.data)){
  set.max <- c(set.max, max(density(input.data[[i]]$fitness_rescaled)$y))
}
set.max <- max(set.max)

plot(density(input.data[[1]]$fitness_rescaled), 
     bty = "n", 
     col = viridis(n = length(input.data), alpha = 0.5)[1],
     lwd = 3,
     main = "", 
     cex.main = 3, 
     cex.axis = 1.8,
     cex.lab = 2.5, 
     xlab = "Fitness (re-scaled)", 
     xlim = c(-1.5, 0.5),
     ylim = c(0, set.max))

for(i in 2:length(input.data)){
  lines(density(input.data[[i]]$fitness_rescaled), 
        col = viridis(n = 11, alpha = 0.5)[i],
        lwd = 3)
}

### legend
legend('topleft',
       legend = c(paste0(round(as.numeric(names(input.data)[1:length(input.data)]),3), " nM")),
       col = viridis(n = length(input.data), alpha = 0.9),
       lwd = 3,
       cex = 1.5,
       bty = "n")

dev.off()


## clean-up environment
rm(dosages, input.stops, input.WTs, 
   scaling_data, coefficients, tmp.lm,
   set.max, i)


## 3. Produce MoCHI input file for the Tite-MoCHI run ##
########################################################

## Only export the key columns of interest
key.columns <- c("aa_seq", "Nham_aa", "WT", "fitness_rescaled", "sigma_rescaled")
input.data.tMoCHI <- lapply(input.data, function(x){x <- x[,match(key.columns, colnames(x))]; colnames(x)[4:5] <- c("fitness", "sigma"); return(x)})

## export for Tite-MoCHI
for(i in 1:length(input.data.tMoCHI)){
  
  write.table(input.data.tMoCHI[[i]], 
              file = paste0("Tite-MoCHI_input_conc_", names(input.data.tMoCHI)[i], "_rescaled.txt"), 
              sep = "\t", quote = F, na = "", row.names = F)
  
}


## 4. Produce MoCHI input file for a "dummy double mutant" run ##
#################################################################

## Only export the key columns of interest
input.data.ddMoCHI <- lapply(input.data, function(x){x <- x[,match(key.columns, colnames(x))]; colnames(x)[4:5] <- c("fitness", "sigma"); return(x)})

## Choose one of the wildtypes (here: highest concentration)
WT.max <- input.data.ddMoCHI[[1]][which(input.data.ddMoCHI[[1]]$WT == T),,drop = F]

## Remove all main WT variants
input.data.ddMoCHI <- lapply(input.data.ddMoCHI, function(x){x <- x[-which(x$WT == T),]; return(x)})

## Encode each concentration as a mutation (should work until 19 concentrations)
alphabet <- c("G", "A", "V", "L", "M", "I", "F", 
              "Y", "W", "K", "R", "H", "D", "E", 
              "S", "T", "C", "N", "Q", "P")

for(i in 1:length(input.data.ddMoCHI)){
  
  ## add dummy AA 
  input.data.ddMoCHI[[i]]$aa_seq <- paste0(input.data.ddMoCHI[[i]]$aa_seq, alphabet[i])
  
  ## add mutant count
  input.data.ddMoCHI[[i]]$Nham_aa <- input.data.ddMoCHI[[i]]$Nham_aa + 1
  
}

## Combine
input.data.ddMoCHI <- do.call(rbind, input.data.ddMoCHI)

## Re-add the chosen WT (add proline as its unique extension)
WT.max$aa_seq <- paste0(WT.max$aa_seq, alphabet[20])
input.data.ddMoCHI$WT <- NA
input.data.ddMoCHI <- rbind(WT.max, input.data.ddMoCHI)

write.table(input.data.ddMoCHI, 
            file = "ddMoCHI_input_rescaled.txt", 
            sep = "\t", quote = F, na = "", row.names = F)


## 5. Version ##
################

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
# other attached packages:
# [1] viridis_0.6.5     viridisLite_0.4.2 scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.4     gtable_0.3.6       glue_1.8.0         gridExtra_2.3      tibble_3.3.0       pkgconfig_2.0.3   
# [11] generics_0.1.4     dplyr_1.1.4        lifecycle_1.0.4    ggplot2_4.0.0      cli_3.6.5          S7_0.2.0           grid_4.5.1         vctrs_0.6.5        compiler_4.5.1     tools_4.5.1       
# [21] pillar_1.11.1      rlang_1.1.6        stringi_1.8.7