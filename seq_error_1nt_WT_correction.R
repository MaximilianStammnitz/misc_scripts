## 1nt sequencing error correction in nf-core/deepmutscan
## incl. count heatmap re-making
## based on WT sequencing
## 16.10.2025
## maximilian.stammnitz@crg.eu

## 0. Environment ##
####################

## libs
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)
library(Biostrings)

## paths
setwd("")

## define input folder
input <- "" ## adjust


## 1. Error correction at the nucleotide level ##
#################################################

## function to run the WT sequencing error correction (nt level)
seq_error_correct_by_WT_nt <- function(wt_seq_count_path, input_count_path, output_file_path){
  
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
  # plot(x = input.counts$counts + 1, 
  #      y = input.counts$counts_corrected + 1, 
  #      xlim = c(1,10000), ylim = c(1,10000),
  #      xlab = "Input counts (original)", ylab = "Input counts (error-corrected)",
  #      pch = 16, bty = "n", cex = 0.3, cex.lab = 1.3, log = "xy")
  
  # Write the processed data
  write.csv(input.counts, file = output_file_path)
  
}

## render
seq_error_correct_by_WT_nt(wt_seq_count_path = "GID1A_WT/intermediate_files/processed_gatk_files/GID1A_input_1_pe/variantCounts_filtered_by_library.csv",
                           input_count_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library.csv"),
                           output_file_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library_err_corrected.csv"))


## 2. Propagate to the amino acid level & make error-corrected count heatmaps ##
################################################################################

## to re-make the AA level input table after WT sequencing error correction
prepare_gatk_data_for_counts_heatmaps_error_corrected <- function(gatk_file_path, aa_seq_file_path, output_csv_path, threshold = 3) {
  
  # Load the raw GATK data
  raw_gatk <- read.table(gatk_file_path, sep = ",", header = TRUE)
  
  # Read the wild-type amino acid sequence from the text file
  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))  # Split the sequence into individual amino acids
  
  # Summarize counts-per-cov for each unique aa mutation in pos_mut
  aggregated_data <- raw_gatk %>%
    group_by(pos_mut) %>%
    summarize(total_counts_per_cov_corrected = sum(counts_per_cov_corrected, na.rm = TRUE),
              total_counts_corrected = sum(counts_corrected, na.rm = TRUE))  # Also sum the counts
  
  # Extract the wild-type position and mutations from 'pos_mut'
  aggregated_data <- aggregated_data %>%
    mutate(
      wt_aa = sub("(\\D)(\\d+)(\\D)", "\\1", pos_mut),  # Wild-type amino acid (e.g., S)
      position = as.numeric(sub("(\\D)(\\d+)(\\D)", "\\2", pos_mut)),  # Position (e.g., 3)
      mut_aa = sub("(\\D)(\\d+)(\\D)", "\\3", pos_mut)   # Mutant amino acid (e.g., R)
    )
  
  # Replace 'X' with '*', indicating the stop codon
  aggregated_data <- aggregated_data %>%
    mutate(mut_aa = ifelse(mut_aa == "X", "*", mut_aa))
  
  # Replace 'X' with '*' in the wild-type amino acid sequence as well
  wt_seq <- ifelse(wt_seq == "X", "*", wt_seq)
  
  # Define all 20 standard amino acids and the stop codon "*"
  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
  
  # Create a list of all positions in the wild-type sequence
  all_positions <- 1:length(wt_seq)
  
  # Create a complete grid of all possible combinations of positions and amino acids
  complete_data <- expand.grid(mut_aa = all_amino_acids, position = all_positions)
  
  # Merge the summarized data with the complete grid (filling missing entries with 0)
  heatmap_data <- complete_data %>%
    left_join(aggregated_data, by = c("mut_aa", "position")) %>%
    mutate(total_counts_per_cov_corrected = ifelse(is.na(total_counts_per_cov_corrected), 0, total_counts_per_cov_corrected),
           wt_aa = wt_seq[position])  # Assign the wild-type amino acid
  
  # Set variants with counts < threshold to NA
  heatmap_data <- heatmap_data %>%
    mutate(
      total_counts_per_cov_corrected = ifelse(total_counts_corrected < threshold, NA, total_counts_per_cov_corrected),
      total_counts_corrected = ifelse(total_counts_corrected < threshold, NA, total_counts_corrected)
    )
  
  # Fill pos_mut column
  heatmap_data <- heatmap_data %>%
    mutate(
      pos_mut = ifelse(is.na(pos_mut),
                       paste0(wt_aa, position, mut_aa),
                       pos_mut)
    )
  
  # Save the aggregated data to a CSV file
  write.csv(heatmap_data, file = output_csv_path, row.names = FALSE)
  print(paste("Aggregated data saved to:", output_csv_path))
}

## render
prepare_gatk_data_for_counts_heatmaps_error_corrected(gatk_file_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library_err_corrected.csv"),
                                                      aa_seq_file_path = paste0(input, "/intermediate_files/aa_seq.txt"),
                                                      output_csv_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_for_heatmaps_err_corrected.csv"),
                                                      threshold = 3)

## functions to plot the WT sequencing error corrected count heatmaps
counts_heatmap_error_corrected <- function(input_csv_path, threshold = 3, output_pdf_path, img_format = "pdf") {
  
  # Inner function to add padding to the last row, adding 21 amino acids per position
  pad_heatmap_data_long <- function(heatmap_data_long, min_non_na_value, num_positions_per_row = 75) {
    
    all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
    
    max_position <- max(heatmap_data_long$position)
    num_missing_positions <- num_positions_per_row - (max_position %% num_positions_per_row)
    
    if (num_missing_positions < num_positions_per_row) {
      new_positions <- (max_position + 1):(max_position + num_missing_positions)
      
      # Add all 21 amino acid variants for each new position
      padding_data <- expand.grid(
        mut_aa = all_amino_acids,  # All possible amino acids
        position = new_positions  # New positions to be padded
      )
      
      # Set placeholder values for the added positions to the exact smallest non-NA value
      padding_data$total_counts_corrected <- min_non_na_value  # Set to the smallest non-NA value
      padding_data$wt_aa <- "Y"  # Set wild-type amino acid to 'Y'
      padding_data$wt_aa_pos <- paste0("Y", padding_data$position)  # Create wt_aa_pos with correct positions
      padding_data$row_group <- max(heatmap_data_long$row_group)  # Set row group to the current last group
      
      # Add the new padding rows to heatmap_data_long
      heatmap_data_long <- dplyr::bind_rows(heatmap_data_long, padding_data)
    }
    
    return(heatmap_data_long)
  }
  
  # Load the CSV data
  heatmap_data <- read.csv(input_csv_path)
  
  # Check if the necessary column exists in the data
  if (!"total_counts_corrected" %in% colnames(heatmap_data)) {
    stop("The column 'total_counts_corrected' is not found in the data.")
  }
  
  # Create heatmap_data_long by selecting necessary columns
  heatmap_data_long <- heatmap_data %>%
    select(mut_aa, position, total_counts_corrected, wt_aa)  # Use 'total_counts_corrected'
  
  # Find the smallest non-NA value in total_counts_corrected
  min_non_na_value <- min(heatmap_data_long$total_counts_corrected, na.rm = TRUE) ## here use the error corrected counts 
  
  # Group positions by rows (75 positions per row) and calculate row_group
  heatmap_data_long <- heatmap_data_long %>%
    mutate(row_group = ((position - 1) %/% 75) + 1)  # Grouping positions into rows
  
  # Apply padding to add missing positions at the end of the last row, using the calculated min value
  heatmap_data_long <- pad_heatmap_data_long(heatmap_data_long, min_non_na_value)
  
  # Convert positions to numeric, sort them, and create wt_aa_pos for the plot
  heatmap_data_long <- heatmap_data_long %>%
    mutate(position = as.numeric(position)) %>%  # Ensure position is numeric
    arrange(position) %>%  # Sort by position
    mutate(wt_aa_pos = factor(paste0(wt_aa, position), levels = unique(paste0(wt_aa, position))))  # Create sorted factor levels for wt_aa_pos
  
  # Add a column to identify synonymous mutations (where mut_aa == wt_aa)
  heatmap_data_long <- heatmap_data_long %>%
    mutate(synonymous = mut_aa == wt_aa)
  
  # Definiere die korrekte Reihenfolge der Aminosäuren
  amino_acid_order <- c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K",
                        "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Bearbeite heatmap_data_long und erstelle syn_positions gleichzeitig
  syn_positions <- heatmap_data_long %>%
    mutate(mut_aa = factor(mut_aa, levels = amino_acid_order),
           # Berechne die x-Koordinate, die pro Gruppe immer von 1 bis 75 verläuft
           x = as.numeric(factor(wt_aa_pos, levels = unique(wt_aa_pos))) - ((row_group - 1) * 75),
           y = as.numeric(factor(mut_aa, levels = amino_acid_order))) %>%
    filter(synonymous == TRUE)
  
  # Calculate the number of row groups and adjust plot height dynamically
  num_row_groups <- max(heatmap_data_long$row_group)
  plot_height <- num_row_groups * 4
  
  # Set the limits for the color scale, ignoring NA (negative values are now NA)
  min_count <- min(heatmap_data_long$total_counts_corrected, na.rm = TRUE)
  max_count <- max(heatmap_data_long$total_counts_corrected, na.rm = TRUE)
  max_position <- max(heatmap_data$position)
  
  # Create the heatmap plot with explicit handling for positions > max_position
  heatmap_plot <- ggplot(heatmap_data_long, aes(x = wt_aa_pos, y = mut_aa, fill = total_counts_corrected)) +
    scale_fill_gradientn(colours = c(alpha("blue", 0), "blue"), na.value = "grey35", trans = "log",  # Apply log transformation to the scale
                         limits = c(min_count, max_count),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),  # Logarithmic scale breaks
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(labels = function(x) {
      numeric_pos <- as.numeric(gsub("[^0-9]", "", x))
      ifelse(numeric_pos > max_position, " ", x)
    }) +
    geom_tile() +
    
    # Add diagonal lines for synonymous mutations using geom_segment
    geom_segment(data = syn_positions[syn_positions$position <= max_position, ],
                 aes(x = x - 0.485, xend = x + 0.485,
                     y = y - 0.485, yend = y + 0.485, color = synonymous),
                 size = 0.2) +
    
    # Manuelle Farbskala für die diagonalen Linien
    scale_color_manual(values = c("TRUE" = "grey10"), labels = c("TRUE" = "")) +
    
    theme_minimal() +
    labs(title = "Heatmap of Counts per Variant", x = "Wild-type Amino Acid", y = "Mutant Amino Acid", fill = "Counts", color = "Synonymous Mutation") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),  # Larger y-axis labels
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14),  # Larger legend title
          legend.text = element_text(size = 12),  # Larger legend text
          panel.spacing = unit(0.1, "lines"),  # Adjust panel spacing
          strip.text = element_blank(),  # Remove row group labels (facet numbers)
          strip.background = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +  # Remove minor grid lines
    facet_wrap(~ row_group, scales = "free_x", ncol = 1) +  # Group by 75 positions per row
    theme(panel.spacing = unit(0.2, "lines"))
  
  heatmap_plot <- heatmap_plot +
    geom_point(data = heatmap_data_long, aes(size = ""), colour = "black", alpha = 0)  # Invisible points for legend
  heatmap_plot <- heatmap_plot +
    guides(size = guide_legend(paste("Dropout (Counts <", threshold, ")"), override.aes = list(shape = 15, size = 8, colour = "grey35", alpha = 1)))  # Define Legend for Dropouts
  
  # Save the heatmap plot
  if (img_format == "pdf") {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150, device = cairo_pdf)
  } else {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150)
  }
  
  if (file.exists(output_pdf_path)) {
    print("Heatmap image successfully created!")
  } else {
    print("Error: Heatmap image was not created.")
  }
}
counts_per_cov_heatmap_error_corrected <- function(input_csv_path, threshold = 3, output_pdf_path, img_format = "pdf") {
  
  # Inner function to add padding to the last row, adding 21 amino acids per position
  pad_heatmap_data_long <- function(heatmap_data_long, min_non_na_value, num_positions_per_row = 75) {
    all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
    
    max_position <- max(heatmap_data_long$position)
    num_missing_positions <- num_positions_per_row - (max_position %% num_positions_per_row)
    
    if (num_missing_positions < num_positions_per_row) {
      new_positions <- (max_position + 1):(max_position + num_missing_positions)
      
      # Add all 21 amino acid variants for each new position
      padding_data <- expand.grid(
        mut_aa = all_amino_acids,  # All possible amino acids
        position = new_positions  # New positions to be padded
      )
      
      # Set placeholder values for the added positions to the exact smallest non-NA value
      padding_data$total_counts_per_cov_corrected <- min_non_na_value  # Set to the smallest non-NA value
      padding_data$wt_aa <- "Y"  # Set wild-type amino acid to 'Y'
      padding_data$wt_aa_pos <- paste0("Y", padding_data$position)  # Create wt_aa_pos with correct positions
      padding_data$row_group <- max(heatmap_data_long$row_group)  # Set row group to the current last group
      
      # Add the new padding rows to heatmap_data_long
      heatmap_data_long <- dplyr::bind_rows(heatmap_data_long, padding_data)
    }
    
    return(heatmap_data_long)
  }
  
  # Load the CSV data
  heatmap_data <- read.csv(input_csv_path)
  
  # Check if the necessary column exists in the data
  if (!"total_counts_per_cov_corrected" %in% colnames(heatmap_data)) {
    stop("The column 'total_counts_per_cov_corrected' is not found in the data.")
  }
  
  # Create heatmap_data_long by selecting necessary columns
  heatmap_data_long <- heatmap_data %>%
    select(mut_aa, position, total_counts_per_cov_corrected, wt_aa)  # Use 'total_counts_per_cov_corrected'
  
  # Find the smallest non-NA value in total_counts_per_cov_corrected
  min_non_na_value <- min(heatmap_data_long$total_counts_per_cov_corrected, na.rm = TRUE)
  
  # Group positions by rows (75 positions per row) and calculate row_group
  heatmap_data_long <- heatmap_data_long %>%
    mutate(row_group = ((position - 1) %/% 75) + 1)  # Grouping positions into rows
  
  # Apply padding to add missing positions at the end of the last row, using the calculated min value
  heatmap_data_long <- pad_heatmap_data_long(heatmap_data_long, min_non_na_value)
  
  # Convert positions to numeric, sort them, and create wt_aa_pos for the plot
  heatmap_data_long <- heatmap_data_long %>%
    mutate(position = as.numeric(position)) %>%  # Ensure position is numeric
    arrange(position) %>%  # Sort by position
    mutate(wt_aa_pos = factor(paste0(wt_aa, position), levels = unique(paste0(wt_aa, position))))  # Create sorted factor levels for wt_aa_pos
  
  # Add a column to identify synonymous mutations (where mut_aa == wt_aa)
  heatmap_data_long <- heatmap_data_long %>%
    mutate(synonymous = mut_aa == wt_aa)
  
  # Definiere die korrekte Reihenfolge der Aminosäuren
  amino_acid_order <- c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K",
                        "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Bearbeite heatmap_data_long und erstelle syn_positions gleichzeitig
  syn_positions <- heatmap_data_long %>%
    mutate(mut_aa = factor(mut_aa, levels = amino_acid_order),
           # Berechne die x-Koordinate, die pro Gruppe immer von 1 bis 75 verläuft
           x = as.numeric(factor(wt_aa_pos, levels = unique(wt_aa_pos))) - ((row_group - 1) * 75),
           y = as.numeric(factor(mut_aa, levels = amino_acid_order))) %>%
    filter(synonymous == TRUE)
  
  # Calculate the number of row groups and adjust plot height dynamically
  num_row_groups <- max(heatmap_data_long$row_group)
  plot_height <- num_row_groups * 4
  
  # Set the limits for the color scale, ignoring NA (negative values are now NA)
  min_count <- min(heatmap_data_long$total_counts_per_cov_corrected, na.rm = TRUE)
  max_count <- max(heatmap_data_long$total_counts_per_cov_corrected, na.rm = TRUE)
  max_position <- max(heatmap_data$position)
  
  # Create the heatmap plot with explicit handling for positions > max_position
  heatmap_plot <- ggplot(heatmap_data_long, aes(x = wt_aa_pos, y = mut_aa, fill = total_counts_per_cov_corrected)) +
    scale_fill_gradientn(colours = c(alpha("blue", 0), "blue"), na.value = "grey35", trans = "log",  # Apply log transformation to the scale
                         limits = c(min_count, max_count),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),  # Logarithmic scale breaks
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(labels = function(x) {
      numeric_pos <- as.numeric(gsub("[^0-9]", "", x))
      ifelse(numeric_pos > max_position, " ", x)
    }) +
    geom_tile() +
    
    # Add diagonal lines for synonymous mutations using geom_segment
    geom_segment(data = syn_positions[syn_positions$position <= max_position, ],
                 aes(x = x - 0.485, xend = x + 0.485,
                     y = y - 0.485, yend = y + 0.485, color = synonymous),
                 size = 0.2) +
    
    # Manuelle Farbskala für die diagonalen Linien
    scale_color_manual(values = c("TRUE" = "grey10"), labels = c("TRUE" = "")) +
    
    theme_minimal() +
    labs(title = "Heatmap of Counts per Coverage for Mutations", x = "Wild-type Amino Acid", y = "Mutant Amino Acid", fill = "Counts per \n Coverage", color = "Synonymous Mutation") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),  # Larger y-axis labels
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14),  # Larger legend title
          legend.text = element_text(size = 12),  # Larger legend text
          panel.spacing = unit(0.1, "lines"),  # Adjust panel spacing
          strip.text = element_blank(),  # Remove row group labels (facet numbers)
          strip.background = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +  # Remove minor grid lines
    facet_wrap(~ row_group, scales = "free_x", ncol = 1) +  # Group by 75 positions per row
    theme(panel.spacing = unit(0.2, "lines"))
  
  heatmap_plot <- heatmap_plot +
    geom_point(data = heatmap_data_long, aes(size = ""), colour = "black", alpha = 0)  # Invisible points for legend
  heatmap_plot <- heatmap_plot +
    guides(size = guide_legend(paste("Dropout (Counts <", threshold, ")"), override.aes = list(shape = 15, size = 8, colour = "grey35", alpha = 1)))  # Define Legend for Dropouts
  
  # Save the heatmap plot
  if (img_format == "pdf") {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150, device = cairo_pdf)
  } else {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150)
  }
  
  if (file.exists(output_pdf_path)) {
    print("Heatmap image successfully created!")
  } else {
    print("Error: Heatmap image was not created.")
  }
}

## render
counts_heatmap_error_corrected(input_csv_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_for_heatmaps_err_corrected.csv"),
                               threshold = 3,
                               output_pdf_path = paste0(input, "/final_plots/counts_heatmap_err_corrected.pdf"), 
                               img_format = "")
counts_per_cov_heatmap_error_corrected(input_csv_path = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_for_heatmaps_err_corrected.csv"),
                                       threshold = 3,
                                       output_pdf_path = paste0(input, "/final_plots/counts_per_cov_heatmap_err_corrected.pdf"), 
                                       img_format = "")


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
                      gatk_file = paste0(input, "/intermediate_files/processed_gatk_files/variantCounts_filtered_by_library_err_corrected.csv"),
                      output_file_path = paste0(input, "/intermediate_files/dimsum_input_err_corrected.tsv"))


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
