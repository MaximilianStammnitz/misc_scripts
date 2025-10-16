## 1nt sequencing error correction in nf-core/deepmutscan
## based on WT sequencing
## 16.10.2025
## maximilian.stammnitz@crg.eu

## 0. Environment ##
####################

## libs
library(Biostrings)

## paths
setwd("/Users/mstammnitz/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/Desktop/DeepDrug/2_mature_projects/GID1A-GAI-gibberelic_acid/2_full-length_GID1A_SUNi_library/2025-06_13_dmscore_bPCA/")


## 1. Error correction at the amino acid level ##
#################################################

# ## test inputs
# GID1A.aa.input.r1 <- read.csv("GID1A_S2_rep1/intermediate_files/processed_gatk_files/variantCounts_for_heatmaps.csv")
# GID1A.aa.WT <- read.csv("GID1A_WT/intermediate_files/processed_gatk_files/GID1A_input_1_pe/variantCounts_for_heatmaps.csv")
# 
# ## append key columns
# GID1A.aa.input.r1$total_counts_per_cov_corrected <- GID1A.aa.input.r1$total_counts_per_cov
# GID1A.aa.input.r1$total_counts_corrected <- GID1A.aa.input.r1$total_counts
# 
# for (i in 1:nrow(GID1A.aa.input.r1)){
#   
#   ## only look at variant observed in the WT sequencing
#   if(is.na(GID1A.aa.WT[i,"total_counts_per_cov"]) == T){
#     next
#   }
#   
#   ## subtract the observed per-base coverage in the WT sequencing
#   GID1A.aa.input.r1[i,"total_counts_per_cov_corrected"] <- GID1A.aa.input.r1[i,"total_counts_per_cov"] - GID1A.aa.WT[i,"total_counts_per_cov"]
#   
#   ## based on this, also adjust the total_counts (cross-multiplication)
#   ## total_counts_corrected ~ total_counts * c(total_counts_per_cov_corrected / total_counts_per_cov)
#   GID1A.aa.input.r1[i,"total_counts_corrected"] <- GID1A.aa.input.r1[i,"total_counts"] * c(GID1A.aa.input.r1[i,"total_counts_per_cov_corrected"] / GID1A.aa.input.r1[i,"total_counts_per_cov"])
#   
#   ## round to nearest integer, do not allow for negative counts
#   GID1A.aa.input.r1[i,"total_counts_corrected"] <- round(GID1A.aa.input.r1[i,"total_counts_corrected"])
#   if(GID1A.aa.input.r1[i,"total_counts_corrected"] < 0){
#     GID1A.aa.input.r1[i,"total_counts_corrected"] <- 0
#     GID1A.aa.input.r1[i,"total_counts_per_cov_corrected"] <- 0
#   }
#   
# }
# 
# ## sanity check
# plot(x = GID1A.aa.input.r1$total_counts + 1, 
#      y = GID1A.aa.input.r1$total_counts_corrected + 1, 
#      xlim = c(1,10000), ylim = c(1,10000),
#      xlab = "GID1A input rep 1 (original)", ylab = "GID1A input rep1 (error-corrected)",
#      pch = 16, bty = "n", cex = 0.5, cex.lab = 1.3, log = "xy")


## 2. Error correction at the nucleotide level ##
#################################################

## function to run the WT error correction
seq_error_correct_by_WT <- function(wt_seq_count_path, input_count_path, output_file_path){
  
  ## load data (nucleotide-level counts from GATK)
  WT.counts <- read.csv(wt_seq_count_path)
  input.counts <- read.csv(input_count_path)

  ## append key columns
  input.counts$counts_corrected <- input.counts$counts
  input.counts$counts_per_cov_corrected <- input.counts$counts_per_cov
  
  # Process the GATK file
  cat("Sequencing error correction of GATK counts...\n")
  for (i in 1:nrow(WT.counts)){
    
    ## only look at variants observed in both the input and WT sequencing
    tmp.id <- match(WT.counts[i,"base_mut"], input.counts[,"base_mut"])
    if(is.na(tmp.id) == T){
      next
    }
    
    ## subtract the observed per-base coverage in the WT sequencing
    input.counts[tmp.id,"counts_per_cov_corrected"] <- input.counts[tmp.id,"counts_per_cov"] - WT.counts[i,"counts_per_cov"]
    
    ## based on this, also adjust the total_counts (cross-multiplication)
    ## total_counts_corrected ~ total_counts * c(total_counts_per_cov_corrected / total_counts_per_cov)
    input.counts[tmp.id,"counts_corrected"] <- input.counts[tmp.id,"counts"] * c(input.counts[tmp.id,"counts_per_cov_corrected"] / input.counts[tmp.id,"counts_per_cov"])
    
    ## round to nearest integer, do not allow for negative counts
    input.counts[tmp.id,"counts_corrected"] <- round(input.counts[tmp.id,"counts_corrected"])
    if(input.counts[tmp.id,"counts_corrected"] < 0){
      input.counts[tmp.id,"counts_corrected"] <- 0
      input.counts[tmp.id,"counts_per_cov_corrected"] <- 0
    }
    
  }
  
  ## visual sanity check, how does the correction affect the counts (?) 
  plot(x = input.counts$counts + 1, 
       y = input.counts$counts_corrected + 1, 
       xlim = c(1,10000), ylim = c(1,10000),
       xlab = "Input counts (original)", ylab = "Input counts (error-corrected)",
       pch = 16, bty = "n", cex = 0.3, cex.lab = 1.3, log = "xy")
  
  # Write the processed data
  write.csv(input.counts, file = output_file_path)
  
}

## render
seq_error_correct_by_WT(wt_seq_count_path = "GID1A_WT/intermediate_files/processed_gatk_files/GID1A_input_1_pe/variantCounts_filtered_by_library.csv",
                        input_count_path = "GID1A_S2_rep1/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library.csv",
                        output_file_path = "GID1A_S2_rep1/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library_err_corrected.csv")


## 3. Export the error adjusted counts into DiMSum table format ##
##################################################################

## function to render DiMSum input
generate_dimsum_input <- function(wt_seq_path, gatk_file, pos_range, output_file_path){
  
  # Parse the position range
  positions <- unlist(strsplit(pos_range, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])
  
  # Load the wild-type sequence
  seq_data <- readDNAStringSet(filepath = wt_seq_path)
  wt_seq <- seq_data[[1]]  # Extract the sequence
  wt_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)
  
  # Convert wt_seq to a character string
  wt_seq <- as.character(wt_seq)
  
  # Split the wild-type sequence into codons (groups of 3 bases)
  wt_codons <- substring(wt_seq, seq(1, nchar(wt_seq), 3), seq(3, nchar(wt_seq), 3))
  
  # Helper function to process GATK CSVs into count data
  process_gatk_file <- function(gatk_csv) {
    
    # Load the input GATK CSV file
    gatk_data <- read.csv(gatk_csv, stringsAsFactors = FALSE)
    
    # Initialize a data frame for results
    results <- data.frame(
      nt_seq = character(),
      count = numeric(),
      stringsAsFactors = FALSE)
    
    # Iterate over each row in the input data
    for (i in 1:nrow(gatk_data)) {
      
      # Extract the mutation info
      codon_mut <- gatk_data$codon_mut[i]
      counts <- gatk_data$counts_corrected[i] ## here now selecting the error-corrected counts (!)
      
      # Create a mutable copy of the wild-type codons
      mutated_codons <- wt_codons
      
      # Apply the mutation
      mutations <- strsplit(codon_mut, ", ")[[1]]
      for (mutation in mutations) {
        codon_position <- as.numeric(sub(":.*", "", mutation))
        new_codon <- sub(".*>", "", mutation)
        
        # Replace the codon at the specified position
        mutated_codons[codon_position] <- new_codon
      }
      
      # Convert the mutated codons back to a sequence string
      mutated_seq_string <- paste(mutated_codons, collapse = "")
      
      # Add the result to the data frame
      results <- rbind(results, data.frame(nt_seq = mutated_seq_string, count = counts))
    }
    
    return(results)
  }
  
  # Process the GATK file
  cat("Processing GATK file...\n")
  processed_data <- process_gatk_file(gatk_file)
  
  # Write the processed data to a file without column names
  write.table(processed_data, file = output_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}

## render
generate_dimsum_input(wt_seq_path = "GID1A_SUNi_ref_small.fasta",
                      pos_range = "352-1383",
                      gatk_file = "GID1A_S2_rep1/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library_err_corrected.csv",
                      output_file_path = "GID1A_S2_rep1/intermediate_files/dimsum_input_err_corrected.tsv")


## 4. Version ##
################

# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] Biostrings_2.76.0   GenomeInfoDb_1.44.2 XVector_0.48.0      IRanges_2.42.0      S4Vectors_0.46.0   
# [6] BiocGenerics_0.54.0 generics_0.1.4     
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.7              compiler_4.5.1          R6_2.6.1                tools_4.5.1            
# [5] GenomeInfoDbData_1.2.14 rstudioapi_0.17.1       crayon_1.5.3            UCSC.utils_1.4.0       
# [9] jsonlite_2.0.0   