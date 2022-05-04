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

  # extract and append the rest of the list elements into the SCE
  for(i in 2:length(exp_name_list)){

    # load the sce sublib
    exp_name <- exp_name_list[[i]]
    unfil_rds_path <- paste0(output_folder, '/',exp_name ,"/unfiltered/sce_rds_objects/",exp_name,"_sce_unfiltered.rds")
    sce_split <- readRDS(unfil_rds_path)

    # extract the fastq info
    total_reads <- data.frame(sublib_id = sce_split$sub_lib_id[1],
                              total_reads = metadata(sce_split)$library_info[1,1])

    # wipe the metadata
    metadata(sce_split) <- list()

    # combine
    merge_sce <- cbind(merge_sce, sce_split)


  }

  # wipe the sce_split from active memory
  rm(sce_split)

  return(merge_sce)
}
