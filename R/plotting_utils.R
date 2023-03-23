#' @rdname raw_bc_hist
#'
#' @title generate hist from raw bc_data
#'
#' @param sub_lib_fp a filepath to the sublibrary data folder
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return complete exit value
#'
#' @import ggplot2
#' @import readr
#'
#'
raw_bc_hist <- function(sub_lib_fp,
                        output_folder,
                        exp_name
){

  message("loading raw barcodes read counts")

  # Load the raw counts
  bc_stats_df <- readr::read_delim(paste0(sub_lib_fp, paste0("/",exp_name,".BCstats.txt")),
                                  col_names = FALSE,
                                  delim = "\t")
  # Rename the cols
  colnames(bc_stats_df) <- c("BC", "n_reads")

  # Cutoff at 100 reads like in the pipeline
  bc_stats_df_filt <- bc_stats_df %>% dplyr::filter(n_reads > 100)

  # Create histogram of raw barcode read_depth
  gg_hist <- ggplot2::ggplot(bc_stats_df_filt, aes(x=n_reads)) +
              ggplot2::geom_histogram(color="darkblue", fill="lightblue", binwidth = 500) +
              ggplot2::geom_vline(aes(xintercept=median(n_reads)),
                            color="blue", linetype="dashed", size=1) +
              ggplot2::scale_y_continuous(trans='log10') +
              ggplot2::theme_classic()

  ggsave(filename = paste0(output_folder, '/',exp_name ,"/gplots/1_raw_bc_reads.png"),
        plot = gg_hist,
        device = "png")

  return()

}

#' @rdname generate_barcoding_heatmaps
#'
#' @title Generate plate heatmaps of cell barcoding
#'
#' @param sce_split the sce representation of the object
#'
#' @author James Opzoomer \email{james.opzoomer@gmail.com}
#'
#' @return Null - generates the barcoding plate distribution heatmaps
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import stringr
#'
#'
generate_barcoding_heatmaps <- function(sce_split,
                                        output_folder,
                                        exp_name){

  # This is currently cusing errors in with the v2 libraries

  message("Generating QC heatmaps on filtered objects.")

  # Extract the barcoding:
  # 1) Count information
  # 2) UMI information

  # Extract count information by splitting well info by _
  well_locs <- sce_split$well_indexes
  well_locs_list <- stringr::str_split(well_locs, pattern = "_")
  rt_locs <- as.numeric(sapply(well_locs_list, "[[", 1))
  r1_locs <- sapply(well_locs_list, "[[", 2)
  r2_locs <- sapply(well_locs_list, "[[", 3)

  #################### count information #####################################

  # Create a dataframe from the count object
  rt_counts_layout <- data.frame(well_position = 1:48, counts = 0)

  counts <- table(rt_locs)

  # Count the occurrence of the numbers
  rt_counts_layout$counts[match(names(counts), df$well_position)] <- as.numeric(counts)

  # Convert into a x by x matrix for plotting
  # RT is 4x12
  rt_loc_mat <- matrix(rt_counts_layout$counts, nrow = 4, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(rt_counts_layout$counts)), c("gray", "floralwhite", "forestgreen"))


  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/rt_barcoding_layout.png"),
                  res = 400, width = 15, height = 5, units = "cm")

  draw(Heatmap(rt_loc_mat,
                cluster_columns = F,
                cluster_rows = F,
                col = pal_1,
                rect_gp = grid::gpar(col = "black", lwd = 1),
                column_title = "RT cell counts",
                heatmap_legend_param = list(
                title = "cell number"))
  )
  dev.off()

  write_csv(rt_counts_layout, file=paste0(output_folder, '/',exp_name ,"/reports/rt_barcoding_layout.csv"))


  # Ligation 1 is 8x12
  r1_counts_layout <- data.frame("well_position" = NULL, counts = NULL)

  for(i in 1:96){
    count = table(r1_locs)[names(table(r1_locs)) == i] # Count number of elements equal to certain value
    if(length(count)==0){count = 0}
    count_row = data.frame("well_position" = i, counts = count)
    r1_counts_layout = rbind(r1_counts_layout, count_row)
  }

  r1_loc_mat <- matrix(r1_counts_layout$counts, nrow = 8, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(r1_counts_layout$counts)), c("gray", "floralwhite", "forestgreen"))

  message("writing RT heatmap!")
  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/ligation_1_barcoding_layout.png"),
      res = 400, width = 15, height = 10, units = "cm")

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(r1_loc_mat,
               cluster_columns = F,
               cluster_rows = F,
               col = pal_1,
               rect_gp = grid::gpar(col = "black", lwd = 1),
               column_title = "Ligation_1 cell counts",
               heatmap_legend_param = list(
                 title = "cell number"))
  )
  dev.off()

  # Ligation 2 is 8x12
  r2_counts_layout <- data.frame("well_position" = NULL, counts = NULL)

  for(i in 1:96){
    count = table(r2_locs)[names(table(r2_locs)) == i] # Count number of elements equal to certain value
    if(length(count)==0){count = 0}
    count_row = data.frame("well_position" = i, counts = count)
    r2_counts_layout = rbind(r2_counts_layout, count_row)
  }

  r2_loc_mat <- matrix(r2_counts_layout$counts, nrow = 8, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(r2_counts_layout$counts)), c("gray", "floralwhite", "forestgreen"))

  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/ligation_2_barcoding_layout.png"),
      res = 400, width = 15, height = 10, units = "cm")

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(r2_loc_mat,
               cluster_columns = F,
               cluster_rows = F,
               col = pal_1,
               rect_gp = grid::gpar(col = "black", lwd = 1),
               column_title = "Ligation_2 cell counts",
               heatmap_legend_param = list(
                 title = "cell number"))
  )
  dev.off()

  ligation_map_counts <- cbind(r1_counts_layout,
                               r2_counts_layout$counts)

  colnames(ligation_map_counts) <- c("well_id", "ligation_1_cell_count", "ligation_2_cell_count")

  write_csv(ligation_map_counts, file=paste0(output_folder, '/',exp_name ,"/reports/ligation_barcoding_layout.csv"))

  #################### transcript information #################################

  # Combine the UMI count with the locs
  umi_rt <- data.frame(rt_id = rt_locs, umi = sce_split$total)

  # Summarize median by id
  med_umi_rt <- umi_rt %>% group_by(rt_id) %>% summarise(median(umi))

  # extract the counts, reorder and identify and fill in 0s
  rt_umi_layout <- data.frame("well_position" = NULL, counts = NULL)

  for(i in 1:48){
    umi_count = med_umi_rt$`median(umi)`[which(med_umi_rt$rt_id == i)] # Count number of elements equal to certain value
    if(length(umi_count)==0){umi_count = 0}
    count_row = data.frame("well_position" = i, counts = umi_count)
    rt_umi_layout = rbind(rt_umi_layout, count_row)
  }

  rt_umi_mat <- matrix(rt_umi_layout$counts, nrow = 4, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(rt_umi_layout$counts)), c("gray", "floralwhite", "darkorange"))

  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/rt_umi_layout.png"),
      res = 400, width = 15, height = 5, units = "cm")

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(rt_umi_mat,
               cluster_columns = F,
               cluster_rows = F,
               col = pal_1,
               rect_gp = grid::gpar(col = "black", lwd = 1),
               column_title = "RT median UMI per well",
               heatmap_legend_param = list(
                 title = "Median UMI per well"))
  )
  dev.off()

  ### This isnt working properly
  ## LIGATION UMI HEATMAPS
  ## Ligation 1
  # Combine the UMI count with the locs
  umi_r1 <- data.frame(r1_id = r1_locs, umi = sce_split$total)

  # Summarize median by id
  med_umi_r1 <- umi_r1 %>% group_by(r1_id) %>% summarise(median(umi))

  # extract the counts, reorder and identify and fill in 0s
  r1_umi_layout <- data.frame("well_position" = NULL, counts = NULL)

  for(i in 1:96){
    umi_count = med_umi_r1$`median(umi)`[which(med_umi_r1$r1_id == i)] # Count number of elements equal to certain value
    if(length(umi_count)==0){umi_count = 0}
    count_row = data.frame("well_position" = i, counts = umi_count)
    r1_umi_layout = rbind(r1_umi_layout, count_row)
  }

  r1_umi_mat <- matrix(r1_umi_layout$counts, nrow = 8, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(r1_umi_layout$counts)), c("gray", "floralwhite", "darkorange"))

  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/ligation_1_umi_layout.png"),
      res = 400, width = 15, height = 10, units = "cm")

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(r1_umi_mat,
               cluster_columns = F,
               cluster_rows = F,
               col = pal_1,
               rect_gp = grid::gpar(col = "black", lwd = 1),
               column_title = "Ligation_1 median UMI per well",
               heatmap_legend_param = list(
                 title = "Median UMI per well"))
  )
  dev.off()

  ## Ligation 2
  # Combine the UMI count with the locs
  umi_r2 <- data.frame(r2_id = r2_locs, umi = sce_split$total)

  # Summarize median by id
  med_umi_r2 <- umi_r2 %>% group_by(r2_id) %>% summarise(median(umi))

  # extract the counts, reorder and identify and fill in 0s
  r2_umi_layout <- data.frame("well_position" = NULL, counts = NULL)

  for(i in 1:96){
    umi_count = med_umi_r2$`median(umi)`[which(med_umi_r2$r2_id == i)] # Count number of elements equal to certain value
    if(length(umi_count)==0){umi_count = 0}
    count_row = data.frame("well_position" = i, counts = umi_count)
    r2_umi_layout = rbind(r2_umi_layout, count_row)
  }

  r2_umi_mat <- matrix(r2_umi_layout$counts, nrow = 8, ncol = 12, byrow = TRUE)

  # Counts pal
  pal_1 <- colorRamp2(c(0,1,max(r2_umi_layout$counts)), c("gray", "floralwhite", "darkorange"))

  # Heatmap
  png(file=paste0(output_folder, '/',exp_name ,"/gplots/ligation_2_umi_layout.png"),
      res = 400, width = 15, height = 10, units = "cm")

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(r2_umi_mat,
               cluster_columns = F,
               cluster_rows = F,
               col = pal_1,
               rect_gp = grid::gpar(col = "black", lwd = 1),
               column_title = "Ligation_2 median UMI per well",
               heatmap_legend_param = list(
                 title = "Median UMI per well"))
  )
  dev.off()

  return()
}

