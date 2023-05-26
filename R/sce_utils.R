#' @rdname gen_split_sce
#'
#' @title generate the sce
#'
#' @param sub_lib_fp a filepath representing the sublibrary filepath
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE output and reports from the split-seq run
#'
#' @import SingleCellExperiment
#' @import readr
#' @import Matrix
#'
# To do:
# Need to auto find the relevant files

gen_split_sce <- function(
  sub_lib_fp,
  output_folder,
  inex_mode = "inex",
  ds_dge_mtx = "all",
  dge_mtx,
  gene_names,
  exp_name,
  library_stats_df
){

  message("loading UMI and Reads dge matrices")
  # Select the list
  dge_mtx_list <- paste0(sub_lib_fp, "/zUMIs_output/expression/",dge_mtx)

  # Select the specific umi and read matrices from the list
  dge_mtx_in <- readRDS(dge_mtx_list)
  dge_umi <- dge_mtx_in$umicount[[inex_mode]][[ds_dge_mtx]]
  dge_reads <- dge_mtx_in$readcount[[inex_mode]][[ds_dge_mtx]]

  # Load the gene names
  gene_names_df <- readr::read_delim(paste0(sub_lib_fp, "/zUMIs_output/expression/", gene_names), delim = "\t")

  message("Generating SCE from DGE matrix")

  col_data <- data.frame(sub_lib_id = rep(exp_name, ncol(dge_umi)))

  # Create the SCE object
  sce_split <- SingleCellExperiment::SingleCellExperiment(list(counts=dge_umi, reads=dge_reads),
                                                          colData = col_data,
                                                          metadata = list(library_info = library_stats_df))

  message("SCE created!")

  # Stash the Gene IDs and re-assign gene names
  message("Stashing gene_ids for gene names")
  row_data <- data.frame(gene_ids = rownames(sce_split))
  rowData(sce_split) <- row_data
  gene_names_df$gene_name <- tolower(gene_names_df$gene_name)
  gene_names_df$gene_name <- make.names(gene_names_df$gene_name, unique = TRUE)

  # Match the names to the vector
  match_vector <- match(rownames(sce_split), gene_names_df$gene_id)

  # Assign the names
  rownames(sce_split) <- gene_names_df$gene_name[match_vector]

  # Write the SCE to file
  message("writing unfiltered SCE to file in output folder")

  # Create output structure
  dir.create(file.path(paste0(output_folder,'/',exp_name ,  "/unfiltered/sce_rds_objects")))
  dir.create(file.path(paste0(output_folder,'/',exp_name ,  "/unfiltered/h5ad_objects")))

  # Return the SCE
  return(sce_split)
}

#' @rdname label_sce_data
#'
#' @title intake raw SCE and label with cell metadata as colData
#'
#' @param sce_split a SCE object representing the processed sublibrary
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return a SCE with labelled barcode well and sample info integrated.
#'
#' @import SingleCellExperiment
#' @import readr
#' @import readxl
#' @import scater
#' @import dplyr
#'
label_sce_data <- function(sce_split,
                            rt_bc_map,
                            lig_bc,
                            sample_map,
                            output_folder,
                            exp_name,
                            sl_index
){

  # Add per cell QC data
  sce_split <- scater::addPerCellQC(sce_split)

  # Barcode-info
  # Barcoding layouts
  rt_bc_layout <- readr::read_csv(file = rt_bc_map)
  lig_bc <- readr::read_csv(file = lig_bc)
  rt_bc_pos <- c(17,24)
  r2_bc_pos <- c(9,16)
  r3_bc_pos <- c(1,8)

  # Sample layout
  sample_rt_layout <- readxl::read_xlsx(sample_map)

  # Extract BCs from the seuruat object
  cell_barcodes <- colnames(sce_split)

  # extract the RT sections
  bc_split <- stringr::str_split(cell_barcodes, "")
  rt_bc <-  sapply(bc_split,
                   FUN = function(x) paste(x[rt_bc_pos[1]:rt_bc_pos[2]], collapse = ""))

  # Match barcodes to the rt_bc_map
  polya <- rt_bc_layout$bc[1:48]
  rhex <- rt_bc_layout$bc[49:96]

  # Match and merge the two RT codes
  mm_polya <- match(rt_bc, polya)
  mm_hex <- match(rt_bc, rhex)
  mm <- dplyr::coalesce(mm_polya, mm_hex)

  # Label the sce object with sample_id
  cell_idents <- sample_rt_layout$sample_id[mm]

  # Add new colData
  sce_split$sample_id <- cell_idents

  # Create a cocatenated well pos index
  rt_well_position <- rt_bc_layout$well_position[mm]

  # R2 (ligation)
  r2_bc <-  sapply(bc_split,
                   FUN = function(x) paste(x[r2_bc_pos[1]:r2_bc_pos[2]], collapse = ""))

  mm_r2 <- match(r2_bc, lig_bc$bc)

  r2_well_position <- lig_bc$well_position[mm_r2]

  # R3 (ligation)
  r3_bc <-  sapply(bc_split,
                   FUN = function(x) paste(x[r3_bc_pos[1]:r3_bc_pos[2]], collapse = ""))

  mm_r3 <- match(r3_bc, lig_bc$bc)

  r3_well_position <- lig_bc$well_position[mm_r3]

  # Sublibrary processing index
  sub_lib_idx_vec <- rep(exp_name, length(r3_well_position))

  # Cocat the names
  cocat_well_pos <- paste(rt_well_position, r2_well_position, r3_well_position, sep = "_")

  # Well pos with sublib_id
  pcr_bc4 <- paste0("sl", sl_index)
  full_bc <- paste(cocat_well_pos, pcr_bc4, sep = "_")

  # Add new colData
  sce_split$well_indexes <- cocat_well_pos

  # Swap the barcode with the index
  sce_split$barcode_seq <- colnames(sce_split)
  colnames(sce_split) <- full_bc

  # Save the object
  saveRDS(sce_split, file = paste0(output_folder, '/',exp_name ,"/unfiltered/sce_rds_objects/",exp_name,"_sce_unfiltered.rds"))
  zellkonverter::writeH5AD(sce_split, file = paste0(output_folder, '/',exp_name ,"/unfiltered/h5ad_objects/",exp_name,"_sce_unfiltered.h5ad"))

  message("Unfiltered SCE objects written to file!!")

  # Return the SCE object
  return(sce_split)


}

#' @rdname filter_split_sce
#'
#' @title intake labelled SCE and filter using dropletUtils
#'
#' @param sce_split
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return a SCE filtered on inctact cell barcodes
#'
#' @import SingleCellExperiment
#' @import DropletUtils
#' @import scales
#'
filter_split_sce <- function(sce_split,
                             output_folder,
                             exp_name,
                             filtering_mode = filtering_mode,
                             filter_value = filter_value){

  message("I am filtering the SCE for intact cells")

  # Extract counts
  transcript_data <- SingleCellExperiment::counts(sce_split)

  # Run the dropletUtils spline fitting function
  br.out <- DropletUtils::barcodeRanks(transcript_data)

  # Extract the data for plotting
  br.out_df <- as.data.frame(br.out)

  if(filtering_mode == "knee"){

    message("filtering by knee cutoff")

    sce_split_filter <- sce_split[, j = (sce_split$sum > metadata(br.out)$inflection)]

    # Waterfall transcripts per cell
    wfall_reads_fil <- ggplot(br.out_df, aes(x = rank, y = total)) +
      geom_point(size=0.8) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_hline(yintercept=metadata(br.out)$knee, linetype="dashed", color="dodgerblue") +
      geom_text(aes(3,metadata(br.out)$knee,label = paste0("Knee point: ", metadata(br.out)$knee), vjust = +1.5), color="dodgerblue") +
      xlab("Log Barcode Rank") + ylab("UMI count per cell")


    ggsave(plot = wfall_reads_fil, paste0(output_folder, '/',exp_name ,"/gplots/3_umi_waterfall.png"),
           device = "png", width = 7.22, height = 7.22)

  }

  else if(filtering_mode == "manual"){

    message("filtering by manual cutoff")

    sce_split_filter <- sce_split[, j = (sce_split$sum > filter_value)]

    # Waterfall transcripts per cell
    wfall_reads_fil <- ggplot(br.out_df, aes(x = rank, y = total)) +
      geom_point(size=0.8) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_hline(yintercept=filter_value, linetype="dashed", color="dodgerblue") +
      geom_text(aes(3,filter_value,label = paste0("Filter point: ", filter_value), vjust = +1.5), color="dodgerblue") +
      xlab("Log Barcode Rank") + ylab("UMI count per cell")


    ggsave(plot = wfall_reads_fil, paste0(output_folder, '/',exp_name ,"/gplots/3_umi_waterfall.png"),
           device = "png", width = 7.22, height = 7.22)

  }else{message("Filtering condition not recognised!")
         message("Please use either knee or manual (with manual tscp cutoff)")}

  # return the filtered sce
  return(sce_split_filter)

}

#' @rdname write_sce_split_lab_filt_stats
#'
#' @title intake complete SCE and write out to sce (.rds) and .H5ad
#'
#' @param sce_split a SCE object representing the processed sublibrary
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return null - writes sce data to file
#'
#' @import SingleCellExperiment
#' @import zellkonverter
#'
#'
write_sce_split_lab_filt_stats <- function(sce_split,
                                           output_folder,
                                           exp_name){


  # Create output structure
  dir.create(file.path(paste0(output_folder,'/',exp_name ,  "/filtered/sce_rds_objects")))
  dir.create(file.path(paste0(output_folder,'/',exp_name ,  "/filtered/h5ad_objects")))

  # Save the complete objects
  saveRDS(sce_split, file = paste0(output_folder, '/',exp_name ,"/filtered/sce_rds_objects/",exp_name,"_sce_filtered.rds"))
  zellkonverter::writeH5AD(sce_split, file = paste0(output_folder, '/',exp_name ,"/filtered/h5ad_objects/",exp_name,"_sce_filtered.h5ad"))

  message("Annotated SCE objects written to file!!")
  return()

}

#' @rdname split_stats_downsample_rds
#'
#' @title generate the sce downsampling statistics
#'
#' @param sce_split a SCE object representing the processed sublibrary
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE output and reports from the split-seq run
#'
#' @import SingleCellExperiment
#' @import readr
#' @import Matrix
#'
#'
split_stats_downsample_rds <- function(sce_split,
                                       output_folder,
                                       exp_name){

  # Read in the downsampling vector



  return(sce_split)
}
