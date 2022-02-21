#' @rdname run_split_pipe
#'
#' @title Run the split-seq pipeline
#'
#' @param
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE output and reports from the split-seq run
#'
#' @import SingleCellExperiment
#' @export

# Probably don't need an object for each sample?? remove this feature.
# How am I going to integrate cite-seq data?


run_split_pipe <- function(
  mode = 'single',
  n_sublibs = 1,
  data_folder,
  output_folder,
  filtering_mode = "knee",
  fastq_path, # change to extraction nextSeq_log
  rt_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  lig_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  sample_map = "../test_data_sp_5_miseq/cell_metadata.xlsx"
){

  # If run in single check for single data folder
  if(mode == 'single'){

    message("Running in single mode - single sublibrary input")
    message("Checking for directory structure")

    data_folder_abs = normalizePath(data_folder)

    output_folder_abs = normalizePath(output_folder)

    # Count the libraries try/catch if doesn't exist
    dirs = list.dirs(path = data_folder_abs, recursive = FALSE)

    # This also needs to execute a loop through the sublibraries if multiple
    if(length(dirs) != 1){

      warning("Number of possible sublibraries in data_folder is not 1")
      warning("Check data_folder directory structure")
      warning("Run in multi-mode or select a single sublibrary to proces")
      warning("Exiting")

      return()

    } else{

      message(paste0("extracting one sublibrary -- ", dirs[1]))

      # Check the directory structure of the sublibrary path
      # So that all files exits otherwise exit
      # TODO

      # Extract the experiment name from the sublib_folder
      exp_name <- rev(stringr::str_split(dirs[1], pattern = "/")[[1]])[1]

      # Check Arguments if paths exist
      dir.create(file.path(output_folder_abs))
      dir.create(file.path(output_folder_abs, exp_name))
      dir.create(file.path(output_folder_abs, exp_name, "unfiltered"))
      dir.create(file.path(output_folder_abs, exp_name, "filtered"))
      dir.create(file.path(output_folder_abs, exp_name, "gplots"))
      dir.create(file.path(output_folder_abs, exp_name, "reports"))

      # Create pipeline_metadata from a yaml log file
      # Exctract the FASTQ path
      # TODO function

      fastq_path_abs = normalizePath(fastq_path)

      # Extract the total raw read info
      total_reads <- get_read_count(fastq_path = fastq_path_abs)

      # Plot the raw reads -- this is edited out for time
      raw_bc_hist(sub_lib_fp = dirs[1],
                 output_folder = output_folder_abs,
                  exp_name = exp_name)

      # Binning stats visualisation
      # TODO

      # Get the overall run_stats
      library_stats_df <- get_seq_run_info(total_reads = total_reads,
                                       sub_lib_fp = dirs[1],
                                       output_folder = output_folder_abs,
                                       exp_name = exp_name)

      # Get the dge_mtx
      dge_mtx_fp <- paste0(exp_name, ".dgecounts.rds")

      # Get the gene_names
      gene_names_fp <- paste0(exp_name, ".gene_names.txt")

      # Generate unfiltered data
      # Put the run metadata in the SCE_metadata
      # TODO - put into the gen_split_sce
      sce_split = gen_split_sce(sub_lib_fp = dirs[1],
                                output_folder = output_folder_abs,
                                dge_mtx = dge_mtx_fp,
                                gene_names = gene_names_fp,
                                exp_name = exp_name,
                                library_stats_df = library_stats_df)

      # Label the sce
      sce_split_lab <- label_sce_cdata(sce_split = sce_split,
                                      rt_bc_map = rt_bc,
                                      lig_bc = lig_bc,
                                      sample_map = sample_map,
                                      output_folder = output_folder_abs,
                                      exp_name = exp_name)


      # Filter the intact cells in the object

      sce_split_lab_filt <- filter_split_sce(sce_split = sce_split_lab,
                                             output_folder = output_folder_abs,
                                             exp_name = exp_name)

      print("making heatmaps")
      # Make the heatmaps
      generate_barcoding_heatmaps(sce_split = sce_split_lab_filt,
                                  output_folder = output_folder_abs,
                                  exp_name = exp_name)

      print("making stats")
      # Generate read level and cell statistics summary
      # Taking the filtered sce as input
      sce_split_lab_filt_stats <- split_stats_output(sce_split = sce_split_lab_filt,
                                                     output_folder = output_folder_abs,
                                                     exp_name = exp_name)


      # Make downsampling data and embed into the sce_bject
      # Write the downsampling data to file

      # split_stats_write_function
      write_sce_split_lab_filt_stats(sce_split = sce_split_lab_filt_stats,
                                     output_folder = output_folder_abs,
                                     exp_name = exp_name)


      # create html report

    }

  } else if(mode == 'multi'){

    message("Running in multi mode - mutli sublibrary input")
    message("Checking for directory structure")

    # Check is the number of sublibraries match the input

    # Run in multi-mode

  } else {

    warning("Mode not recognised please choose from - sinlge or multi")
    warning("Single - one sublibrary")
    warning("Multi - mutliple sublibraries")
    warning("Exiting")

    return()
  }

  # Check if inputs are correct



}
