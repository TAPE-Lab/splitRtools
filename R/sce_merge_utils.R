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

merge_sce_sublibs <- function(exp_name_list,
                              output_folder){

  merge_out_dir <- "sub_lib_merged"

  # Initialize empty list for filling
  merge_counts_mtx_list <- list()
  merge_reads_mtx_list <- list()

  # ColData dataframe
  sce_coldata <- data.frame()

  # Create new the metadata
  total_reads <- data.frame(sublib_reads = c())

  # extract and append the rest of the list elements into the SCE
  for(i in 1:length(exp_name_list)){

    message(paste0("collecting sublib: ", exp_name_list[[i]]," data!"))

    # load the sce sublibs as a list
    exp_name <- exp_name_list[[i]]
    unfil_rds_path <- paste0(output_folder, '/',exp_name ,"/unfiltered/sce_rds_objects/",exp_name,"_sce_unfiltered.rds")
    sce_split <- readRDS(unfil_rds_path)

    # Create a list of counts mtx to be merged
    merge_counts_mtx_list <- append(merge_counts_mtx_list, SummarizedExperiment::assay(sce_split, 'counts'))
    merge_reads_mtx_list <- append(merge_reads_mtx_list, SummarizedExperiment::assay(sce_split, 'reads'))

    # Extract the colData
    sce_coldata <- rbind(sce_coldata, as.data.frame(colData(sce_split)))

    # extract reads
    lib_reads <- data.frame(sublib_reads = metadata(sce_split)[[1]][1,1])
    total_reads <- rbind(total_reads, lib_reads)

  }

  # wipe the sce_split from active memory
  rm(sce_split)

  message("Merging sparse matrices")
  # Merge the sparse matrices
  counts_merged_mtx <- merge_sparse(merge_counts_mtx_list)
  reads_merged_mtx <- merge_sparse(merge_reads_mtx_list)

  # Format metadata
  t_reads <- list(total_reads)

  message("Forming new merged SCE")
  # Reform the SCE
  combined_sce <- SingleCellExperiment(assays=list(counts=counts_merged_mtx,
                                                   reads=reads_merged_mtx),
                                       colData=sce_coldata,
                                       metadata=t_reads)

  # Save the object
  saveRDS(combined_sce, file = paste0(output_folder, '/',merge_out_dir ,"/unfiltered/",merge_out_dir,"_sce_unfiltered.rds"))
  zellkonverter::writeH5AD(combined_sce, file = paste0(output_folder, '/',merge_out_dir ,"/unfiltered/",merge_out_dir,"_sce_unfiltered.h5ad"))

  message("Unfiltered MERGED SCE objects written to file!! Almost there!!")

  return(combined_sce)
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


#' @rdname merge_sparse
#'
#' @title Merge sparse matrices for SCE merging
#'
#' @param  sce_list a list of sparse matrices to merge into one
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return a combined sparse matrix object
#'
#' @import Matrix
#'
merge_sparse <- function(mtx_list) {

  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in mtx_list) {

    cnold <- colnames(M)
    rnold <- rownames(M)

    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)

    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- unname(which(M != 0,arr.ind=T))
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,M@x)
  }

  new_matrix <- sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
  return(new_matrix)
}

