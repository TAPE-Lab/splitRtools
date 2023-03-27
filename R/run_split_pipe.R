#' @rdname run_split_pipe
#'
#' @title Run the splitRtools split-seq pipeline
#'
#' @param mode String of either 'single' or 'merge', this will determine whether to process each sublibrary separately or in 'merge' mode do this and then process all the sublibraries combined together.
#' @param n_sublibs Integer providing the number of sublibrary folders in the `data_folder`, this is a control number that will cause a exit if it does not match.
#' @param data_folder File path specifying the absolute path of the folder containing the zUMIs pipeline outputs for the sublibraries. The zUMIs output folder names must be named the same as the zUMIs experiment name specified in the .yaml file.
#' @param output_folder File path specifying where to store the outputs of the pipeline. Folder will be created if it does not exist.
#' @param filtering_mode String of either 'knee' or 'manual' to filter data based in UMI per cell. Knee fits model from `DropletUtils`, 'manual' let user specify the cutoff selected in `filter_value`.
#' @param filter_value Integer specifying the UMI cutoff to filter cells/barcodes by based on UMI per cell if `filtering_mode` is set to 'manual'.
#' @param count_reads Boolean specifying whether to count the reads on R1 of sublibrary folders in the directory `fastq_path`.
#' @param total_reads This currently only works for 'single' mode.
#' @param fastq_path File path specifying a master folder with FASTQ files within folder labelled by sublibrary zUMI experiment name, format like `data_folder`. Only used if `count_reads = TRUE`. The pipeline will then count the number of FASTQ reads in R1 of each folder. This can be very slow.
#' @param rt_bc File path specifying the barcode layout for the RT1 plate, see file included in package for example.
#' @param lig_bc File path specifying the barcode layout for the L2 and L3 plates, see file included in package for example.
#' @param sample_map File path specifying the layout of sample loading in the RT plate, indicating which samples were placed in which wells.
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return SCE outputs, summary figures and reports from zUMIs processing of split-seq FASTQ data in `output_folder` filepath.
#'
#' @import SingleCellExperiment
#' @export

# Bugs
# There is a bug where two lines are added to each SCE lib info metadata slot
# The count_reads needs to be spun into a csv submission, doesn't work with merge.
# Need to the cell_barcodes to include sublib index or there will be conflicts

# TODO

# Test SCE2anndata with py_pickle as may preserve some additional functionality

run_split_pipe <- function(
  mode = 'single',
  n_sublibs = 1,
  data_folder,
  output_folder,
  filtering_mode = "knee",
  filter_value = 1000,
  count_reads = FALSE,
  total_reads,
  fastq_path,
  rt_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  lig_bc = "../test_data_sp_5_miseq/barcodes_v1.csv",
  sample_map = "../test_data_sp_5_miseq/cell_metadata.xlsx"
){
    start_time <- Sys.time()
    message(paste0("Running pipeline config in mode - ",mode))
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

        message("Running in multi sublib mode - sublibrary input analysed separately and then merged!")

      }

      # initialize the exp_name lists for merging later
      exp_name_list <- list()

      for(i in 1:length(dirs)){

      message(paste0("extracting zUMIs data from sublibrary -- ", dirs[i]))

      # Extract the experiment name from the sublib_folder
      exp_name <- rev(stringr::str_split(dirs[i], pattern = "/")[[1]])[1]

      # Collect the exp name for later merging
      exp_name_list_tmp <- list(exp_name)
      exp_name_list <- append(exp_name_list, exp_name_list_tmp)

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

      # If reads to be counted
      if(count_reads == TRUE){

        # Exctract the FASTQ path
        fastq_path_dir <- normalizePath(fastq_path)

        # Append the identical name as the sublib folder
        fastq_path_abs <- paste0(fastq_path_dir, "/", exp_name)

        # Extract the total raw read info
        # Found in general_utils.R
        total_reads <- get_read_count(fastq_path = fastq_path_abs)

      }

      # Plot the raw reads -- this is edited out for time
      # This takes a really long time and need to be modified
      # Not sustainable in current format
      #raw_bc_hist(sub_lib_fp = dirs[1],
      #           output_folder = output_folder_abs,
      #            exp_name = exp_name)

      # Binning stats visualisations
      # TODO

      # Get the overall run_stats
      # Found in general_utils.R
      library_stats_df <- get_seq_run_info(total_reads = total_reads,
                                      sub_lib_fp = dirs[i],
                                      output_folder = output_folder_abs,
                                      exp_name = exp_name)

      # Get the dge_mtx
      dge_mtx_fp <- paste0(exp_name, ".dgecounts.rds")

      # Get the gene_names
      gene_names_fp <- paste0(exp_name, ".gene_names.txt")


      # create the unfiltered SCE object
      # Found in sce_utils.R
      sce_split = gen_split_sce(sub_lib_fp = dirs[i],
                                output_folder = output_folder_abs,
                                dge_mtx = dge_mtx_fp,
                                gene_names = gene_names_fp,
                                exp_name = exp_name,
                                library_stats_df = library_stats_df)

      # Label the sce with sample and sublibrary metadata
      # Found in sce_utils.R
      sce_split_lab <- label_sce_data(sce_split = sce_split,
                                      rt_bc_map = rt_bc,
                                      lig_bc = lig_bc,
                                      sample_map = sample_map,
                                      output_folder = output_folder_abs,
                                      exp_name = exp_name)


      # Filter the intact cells in the object with dropletUTILS
      # Found in sce_utils.R
      sce_split_lab_filt <- filter_split_sce(sce_split = sce_split_lab,
                                             output_folder = output_folder_abs,
                                             exp_name = exp_name,
                                             filtering_mode = filtering_mode,
                                             filter_value = filter_value)

      print("making heatmaps")
      # Make the heatmaps
      # Found in plotting_util.R
      generate_barcoding_heatmaps(sce_split = sce_split_lab_filt,
                                  output_folder = output_folder_abs,
                                  exp_name = exp_name)

      print("making stats")
      # Generate read level and cell statistics summary
      # Taking the filtered sce as input
      # Found in general_utils.R
      sce_split_lab_filt_stats <- split_stats_output(sce_split = sce_split_lab_filt,
                                                     output_folder = output_folder_abs,
                                                     exp_name = exp_name)


      # Make downsampling data and embed into the sce_bject
      # Change this to a UMI/gene vs reads plot
      # sce_split_lab_filt_stats_ds <- split_stats_downsample_rds(sce_split = sce_split_lab_filt,
                                     #   output_folder = output_folder_abs,
                                     #   exp_name = exp_name)

      # split_stats_write_function to write complete objects to file
      # Found in sce_utils.R
      write_sce_split_lab_filt_stats(sce_split = sce_split_lab_filt_stats,
                                     output_folder = output_folder_abs,
                                     exp_name = exp_name)


      # create html report
      # Using something like this
      #rrmarkdown::draft("dashboard.Rmd",
      #                   outfile = "split_tools",
      #                   package = "flexdashboard")

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

    # Take the SCEs that were generate earlier in single mode
    merge_sce <- merge_sce_sublibs(exp_name_list = exp_name_list,
                                   output_folder = output_folder_abs)


    # Filter the intact cells in the object with dropletUTILS
    sce_split_lab_filt <- filter_split_sce(sce_split = merge_sce,
                                           output_folder = output_folder_abs,
                                           exp_name = merge_out_dir,
                                           filtering_mode = filtering_mode,
                                           filter_value = filter_value)

    print("making heatmaps")
    # Make the heatmaps
    # This doesn't work anymore, need a count based approach
    generate_barcoding_heatmaps(sce_split = sce_split_lab_filt,
                                output_folder = output_folder_abs,
                                exp_name = merge_out_dir)


    # split_stats_gen_function
    sce_split_lab_filt_stats <- sce_merge_stats(sce_split = sce_split_lab_filt,
                                                output_folder = output_folder_abs,
                                                exp_name = merge_out_dir)

    # write the complete merged sce to file
    write_sce_split_lab_filt_stats(sce_split = sce_split_lab_filt_stats,
                                   output_folder = output_folder_abs,
                                   exp_name = merge_out_dir)

    # Produce basic Seurat analysis for dynamic investigation
    # create html report here

  }

    end_time <- Sys.time()

    difference <- difftime(end_time, start_time, units='mins')
    print(paste0("Pipeline completed in: ",round(difference, digits = 2), ' minutes'))

    return()

}
