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
#' @param fastq_path A string representing the raw FastQ path
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
                                                     "library_info.csv"))


  return(seq_run_metrics_df)

}


