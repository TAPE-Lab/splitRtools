#' @rdname merge_sce_sublibs
#'
#' @title Merge the individually processed unfiltered sublibraries
#'
#' @param sl_list a list object containing the sublibraries to be merged
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE output of the merged data object
#'
#' @import SingleCellExperiment

merge_sce_sublibs <- function(merge_sce,
                              exp_name_list,
                              output_folder){

  merge_out_dir <- "sub_lib_merged"

  # extract and append the rest of the list elements into the SCE
  for(i in 2:length(exp_name_list)){

    # load the sce sublib
    exp_name <- exp_name_list[[i]]
    unfil_rds_path <- paste0(output_folder, '/',exp_name ,"/unfiltered/sce_rds_objects/",exp_name,"_sce_unfiltered.rds")
    sce_split <- readRDS(unfil_rds_path)

    # extract the fastq info
    total_reads_sublib <- data.frame(sublib_id = sce_split$sub_lib_id[1],
                              total_reads = metadata(sce_split)$library_info[1,1])

    # wipe the metadata
    metadata(sce_split) <- list()

    # combine
    merge_sce <- SingleCellExperiment::cbind(merge_sce, sce_split)

    # Add new fastq data to metadata
    total_reads <- metadata(merge_sce)[[1]]
    total_reads <- rbind(total_reads, total_reads_sublib)
    metadata(merge_sce) <- list(total_reads)

  }

  # wipe the sce_split from active memory
  rm(sce_split)

  # Save the object
  saveRDS(merge_sce, file = paste0(output_folder, '/',merge_out_dir ,"/unfiltered/",merge_out_dir,"_sce_unfiltered.rds"))
  zellkonverter::writeH5AD(merge_sce, file = paste0(output_folder, '/',merge_out_dir ,"/unfiltered/",merge_out_dir,"_sce_unfiltered.h5ad"))

  message("Unfiltered merged SCE objects written to file!!")

  return(merge_sce)
}

#' @rdname sce_merge_stats
#'
#' @title create sequencing run statistics for merged sublibraries
#'
#' @param sce_split An input merged SCE experiment object
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
sce_merge_stats <- function(sce_split,
                               output_folder,
                               exp_name){

  # Count the total mapped reads assigned to cells
  rds_map_total <- sum(assay(sce_split, "reads"))
  sce_split$mapped_reads <- Matrix::colSums(assay(sce_split, "reads"))
  total_reads = sum(metadata(sce_split)[[1]]$total_reads)
  estimated_cell_n = ncol(sce_split)

  # Generate the library level info
  split_run_stats <- data.frame(estimated_cell_n = ncol(sce_split),
                                total_reads = total_reads,
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
                                       total_reads = total_reads,
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

  # Call cell_n plotting function
  celln_df <- as.data.frame(table(sce_split$sub_lib_id, sce_split$sample_id))
  colnames(celln_df) <- c("sublibrary", "sample", "cell_number")

  # Bar xlabels could use 90 degree rotation and smaller text
  bar <- ggplot2::ggplot(data=celln_df, aes(x=sublibrary, y=cell_number, fill=sample)) +
    geom_bar(stat="identity", color="black") +
    theme_classic() +
    theme(axis.text=element_text(size=12)) + ylab('Cell number recovered')

  ggsave(plot = bar, filename = paste0(output_folder, '/',exp_name ,"/gplots/cell_abundance_barplot.png"),
         width = 3, height = 4)

  # Make the reads vs umi + gene plots
  #### To do here::

  return(sce_split)

}



