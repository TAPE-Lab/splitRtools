#' @rdname get_read_count
#'
#' @title extract total raw reads
#'
#' @param fastq_path A string representing the raw FastQ path
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return Total reads from raw fastq file
#'
#' @import ShortRead
#'
#'
get_read_count <- function(fastq_path
){

  message("Counting total reads!")
  fl <- fastq_path

  fq_stats <- ShortRead::countFastq(fl)

  total_reads <- fq_stats[,1]

  return(total_reads)

}

#' @rdname get_seq_run_info
#'
#' @title extract total raw reads
#'
#' @param total_reads The total reads from the sequnecing library
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return Total reads from raw fastq file
#'
#' @import readr
#'
get_seq_run_info <- function(total_reads,
                             sub_lib_fp,
                             output_folder,
                             exp_name
                             ){

  message("Assembling read filtering counts!")

  # Load the raw counts
  bc_stats_df <- readr::read_delim(paste0(sub_lib_fp, paste0("/",exp_name,".BCstats.txt")),
                                     col_names = FALSE,
                                     delim = "\t")
  # Rename the cols
  colnames(bc_stats_df) <- c("BC", "n_reads")
  # Extract the sum of reads passing the fracs
  good_qual_reads <- sum(bc_stats_df$n_reads)

  # Load the kept_binned
  binned_reads <- readr::read_delim(paste0(sub_lib_fp, paste0("/zUMIs_output/",exp_name,"kept_barcodes_binned.txt")),
                    col_names = TRUE,
                    delim = ",")

  kept_binned_reads <- sum(binned_reads$n)

  seq_run_metrics_df <- data.frame(total_reads = total_reads,
                                   quality_100_frac = (good_qual_reads/total_reads),
                                   kept_binned = (kept_binned_reads/total_reads))

  # Write to file
  readr::write_csv(seq_run_metrics_df, file = paste0(output_folder, '/',exp_name ,"/reports/",
                                                     "read_count_library_info.csv"))


  return(seq_run_metrics_df)

}

#' @rdname split_stats_output
#'
#' @title create sequencing run statistics for each sublibrary
#'
#' @param sce_split An input SCE experiment object
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE object with seq statistics embedded in the SCE metadata
#'
#' @import readr
#' @import ggplot2
#' @import RColorBrewer
#' @import Matrix
#'
split_stats_output <- function(sce_split,
                                     output_folder,
                                     exp_name){

    # Count the total mapped reads assigned to cells
    rds_map_total <- sum(assay(sce_split, "reads"))
    sce_split$mapped_reads <- Matrix::colSums(assay(sce_split, "reads"))
    total_reads = metadata(sce_split)[[1]][1,1]
    estimated_cell_n = ncol(sce_split)

    # Generate the library level info
    split_run_stats <- data.frame(estimated_cell_n = ncol(sce_split),
                                  total_reads = metadata(sce_split)[[1]][1,1],
                                  estimated_reads_cells = total_reads/estimated_cell_n,
                                  median_umi = median(sce_split$sum),
                                  median_genes = median(sce_split$detected),
                                  seq_sat = 1-(sum(sce_split$sum)/rds_map_total),
                                  intact_cell_umi_thresh = min(sce_split$sum)
    )

    # Stash the values for later
    stats_rownames <- colnames(split_run_stats)

    # transpose
    split_run_stats <- t(split_run_stats)

    colnames(split_run_stats) <- "all_wells"

    if(length(unique(sce_split$sample_id) > 1)){
    # Perform for each sub-sample
    for(i in 1:length(unique(sce_split$sample_id))){

      # subset the sce
      # Extract indices for sub-setting
      idx <- which(sce_split$sample_id == unique(sce_split$sample_id)[i])

      # Subset based on each sample
      sce_sub <- sce_split[,idx]

      # Create data frame
      split_sample_stats <- data.frame(estimated_cell_n = ncol(sce_sub),
                                         total_reads = metadata(sce_split)[[1]][1,1],
                                         estimated_reads_cells = NA,
                                         median_umi = median(sce_sub$sum),
                                         median_genes = median(sce_sub$detected),
                                         seq_sat = 1-(sum(sce_split$sum)/rds_map_total),
                                         intact_cell_umi_thresh = min(sce_split$sum))

      # transpose and cbind onto the global list
      split_sample_stats <- t(split_sample_stats)

      split_run_stats <- cbind(split_run_stats, split_sample_stats)

    }

    # Rename the column names
    colnames(split_run_stats) <- c("all_wells", unique(sce_split$sample_id))

    }

    split_run_stats <- as.data.frame(split_run_stats)

    row_names <- data.frame(vars = stats_rownames)

    split_run_stats_lab <- cbind(row_names, split_run_stats)

    # Stash in the sce object as index 2
    metadata(sce_split)[["seq_run_info"]] <- split_run_stats

    # write to csv
    readr::write_csv(split_run_stats_lab, file = paste0(output_folder, '/',exp_name ,"/reports/",
                                             "sequencing_stats.csv"))

    cell_sample_df <- split_run_stats[1, colnames(split_run_stats)[-1]]
    cell_sample_df <- as.data.frame(cell_sample_df)
    cell_sample_vec <- as.numeric(cell_sample_df[1,])

    # Call cell_n plotting function
    cell_number_samples <- data.frame(sublibrary = rep(sce_split$sub_lib_id[1], length(colnames(split_run_stats))-1),
                                      sample = colnames(split_run_stats)[2:length(colnames(split_run_stats))],
                                      abundance = cell_sample_vec)


    bar <- ggplot2::ggplot(data=cell_number_samples, aes(x=sublibrary, y=abundance, fill=sample)) +
      geom_bar(stat="identity", color="black") +
      theme_classic() +
      theme(axis.text=element_text(size=15)) + ylab('Cell number recovered')


    ggsave(plot = bar, filename = paste0(output_folder, '/',exp_name ,"/gplots/cell_abundance_barplot.png"),
           width = 3, height = 4)

    ###################################### Read filtering breakdown #############################

    # I AM HERE - - - - need to finish the proportions
    # The extra not mapped data is embedded in the zUMIs outputs
    # Probably in the .rds read counts object
    #rds_not_mapped <- readr::read_csv(file = paste0())

    # Generate the filtering parameters
    split_reads_filtering <- data.frame(
      reads_bc_qual_30 = 1-metadata(sce_split)[[1]][1,2],
      reads_bc_whitelist = 1-metadata(sce_split)[[1]][1,3],
      reads_cdna_not_mapped = NA,
      reads_cnda_intergenic = NA,
      reads_duplicates = ((rds_map_total-sum(sce_split$sum))/metadata(sce_split)[[1]][1,1]),
      reads_umi = (sum(sce_split$sum)/metadata(sce_split)[[1]][1,1])
    )

    split_reads_filtering[1,3] <- 1 - sum(split_reads_filtering[1,], na.rm = T)

    # Transpose dataframe
    split_reads_filtering_t <- t(split_reads_filtering)

    # Rename
    colnames(split_reads_filtering_t) <- exp_name

    # stash in the sce at index 3
    metadata(sce_split)[["library_stats"]] <- split_reads_filtering_t

    split_reads_filtering_gplot <- data.frame(library_proportion = split_reads_filtering_t[,exp_name],
                                              sublibrary = rep(exp_name, nrow(split_reads_filtering_t)),
                                              quality_filter = rownames(split_reads_filtering_t))

    # Call plotting function
    bar2 <- ggplot2::ggplot(data=split_reads_filtering_gplot, aes(x=sublibrary, y=library_proportion, fill=quality_filter)) +
      geom_bar(stat="identity", color="black") +
      scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) + theme_classic() +
      theme(axis.text=element_text(size=15)) + ylab('Cell number recovered')

    # write data to csv
    ggsave(plot = bar2, filename = paste0(output_folder, '/',exp_name ,"/gplots/seq_run_filtering.png"),
           width = 4, height = 4)



    readr::write_csv(split_reads_filtering, file = paste0(output_folder, '/',exp_name ,"/reports/",
                                                    "seq_run_filtering.csv"))

    return(sce_split)

}




