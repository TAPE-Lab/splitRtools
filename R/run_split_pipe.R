#' @rdname run_split_pipe
#'
#' @title Run the split-seq pipeline
#'
#' @param mode a string of either 'single' or 'merge'
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE output and reports from zUMIs processing of split-seq FASTQ data
#'
#' @import SingleCellExperiment
#' @export

# Bugs
# There is a bug where two lines are added to each SCE lib info metadata slot

# TODO
# Could introduce some parrallelism into the single-mode

run_split_pipe <- function(
  mode = 'single',
  n_sublibs = 1,
  data_folder,
  output_folder,
  filtering_mode = "knee",
  filter_value = 1000,
  fastq_path, # change to extraction nextSeq_log
  rt_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  lig_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  sample_map = "../test_data_sp_5_miseq/cell_metadata.xlsx"
){

    message("Running in single mode - single sublibrary input")
    message("Checking for directory structure")

    data_folder_abs = normalizePath(data_folder)

    output_folder_abs = normalizePath(output_folder)

    # Count the libraries try/catch if doesn't exist
    dirs = list.dirs(path = data_folder_abs, recursive = FALSE)
    print(paste0("There are: ",length(dirs), " sublibraries detected!"))


    if(length(dirs) != n_sublibs){

      warning("Number of possible sublibraries in data_folder is not 1")
      warning("Check data_folder directory structure")
      warning("Run in multi-mode or select a single sublibrary to proces")
      warning("Exiting")

      return()

      # This loop executes a loop through the sublibraries in data folder
    } else{

      if(mode == 'single'){

        message("Running in single sublib mode - sublibrary input analysed separately!")
      }

      # If the merging the sublibs, create a list of unfiltered SCEs
      # These will be used later to add in the SCEs processed one at a time
      if(mode == 'merge'){

        sublib_list <- list()

      }

      for(i in 1:length(dirs)){

      message(paste0("extracting raw data from sublibrary -- ", dirs[i]))

      # Check the directory structure of the sublibrary path
      # So that all files exits otherwise exit
      # TODO

      # Extract the experiment name from the sublib_folder
      exp_name <- rev(stringr::str_split(dirs[i], pattern = "/")[[1]])[1]

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

      fastq_path_dir <- normalizePath(fastq_path)

      fastq_path_abs <- paste0(fastq_path_dir, "/", exp_name)

      # Extract the total raw read info
      total_reads <- get_read_count(fastq_path = fastq_path_abs)

      # Plot the raw reads -- this is edited out for time
      # This takes a really long time and need to be modified
      # Not sustainable in current format
      #raw_bc_hist(sub_lib_fp = dirs[1],
      #           output_folder = output_folder_abs,
      #            exp_name = exp_name)

      # Binning stats visualisations
      # TODO

      # Get the overall run_stats
      library_stats_df <- get_seq_run_info(total_reads = total_reads,
                                       sub_lib_fp = dirs[i],
                                       output_folder = output_folder_abs,
                                       exp_name = exp_name)

      # Get the dge_mtx
      dge_mtx_fp <- paste0(exp_name, ".dgecounts.rds")

      # Get the gene_names
      gene_names_fp <- paste0(exp_name, ".gene_names.txt")


      # create the unfiltered SCE object
      sce_split = gen_split_sce(sub_lib_fp = dirs[i],
                                output_folder = output_folder_abs,
                                dge_mtx = dge_mtx_fp,
                                gene_names = gene_names_fp,
                                exp_name = exp_name,
                                library_stats_df = library_stats_df)

      # Label the sce with sample and sublibrary metadata
      sce_split_lab <- label_sce_cdata(sce_split = sce_split,
                                      rt_bc_map = rt_bc,
                                      lig_bc = lig_bc,
                                      sample_map = sample_map,
                                      output_folder = output_folder_abs,
                                      exp_name = exp_name)

      # add the unfilteres SCEs to the 'sublib_list' here for mergeing later
      if(mode == 'merge'){

        sublib_list_tmp <- list(sce_split_lab)
        sublib_list <- append(sublib_list, sublib_list_tmp)

      }

      # Filter the intact cells in the object with dropletUTILS
      sce_split_lab_filt <- filter_split_sce(sce_split = sce_split_lab,
                                             output_folder = output_folder_abs,
                                             exp_name = exp_name,
                                             filtering_mode = filtering_mode,
                                             filter_value = filter_value)

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
      # Change this to a UMI/gene vs reads plot
      # sce_split_lab_filt_stats_ds <- split_stats_downsample_rds(sce_split = sce_split_lab_filt,
                                     #   output_folder = output_folder_abs,
                                     #   exp_name = exp_name)

      # split_stats_write_function
      write_sce_split_lab_filt_stats(sce_split = sce_split_lab_filt_stats,
                                     output_folder = output_folder_abs,
                                     exp_name = exp_name)


      # create html report

      }
    }

    if(mode == 'merge'){

    message("Running in merge mode - merge sublibrary input")
    message("Checking for directory structure")

    # Run in merge-mode
    # create the merge output folders
    merge_out_dir <- "sub_lib_merged"
    dir.create(file.path(output_folder_abs, merge_out_dir))
    dir.create(file.path(output_folder_abs, merge_out_dir, "unfiltered"))
    dir.create(file.path(output_folder_abs, merge_out_dir, "filtered"))
    dir.create(file.path(output_folder_abs, merge_out_dir, "gplots"))
    dir.create(file.path(output_folder_abs, merge_out_dir, "reports"))

    # Check if there are less than two sublibraries
    # Cannot merge if there is only 1 so exit
    if(length(dirs) < 2){

      warning("Number of possible sublibraries in data_folder is less than 2")
      warning("Check data_folder directory structure")
      warning("Cannot merge fewer than 2 sublibaries")
      warning("Exiting")

      return()

    }

    # Take the list of SCEs that were generate earlier in single mode
    print(sublib_list)

    # wipe the metadata for re-writing

    # Adjust the metadata stats in the sce

    # Filter the intact cells in the object

    print("making heatmaps")
    # Make the barcoding heatmaps

    # Make merged stats

    # split_stats_write_function

    # Generate initial dynamic UMAP data - TDL

    # create html report


  } else {

    warning("Mode not recognised please choose from - sinlge or merge")
    warning("Single - process sublibraries independantly")
    warning("merge - merge sublibrary output into one sce object")
    warning("Exiting")

    return()
  }
}
