if (getRversion() >= "2.15.1") utils::globalVariables(c(".")) ## Disables notes about '.' due to magrittr

#' Read the data in an eager input TSV file
#'
#' @param tsv_fn character. The path to the input TSV file used in your eager run
#'
#' @return A tibble containing the eager TSV information. This tibble includes three additional columns: initial_merge additional_merge strandedness_clash.
#'    The initial_merge column is logical specifying if merging of libraries of the same UDG/strandedness combination was done after deduplication.
#'    The additional_merge column is logical specifying if merging of libraries of the same strandedness but different UDG was done after bam trimming.
#'    The strandedness_clash column is a logical specifying if libraries of differing strandedness exist that share the same Sample_Name.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_input_tsv_data <- function(tsv_fn) {
  ## End execution if TSV file is not found
  check_file_exists(tsv_fn, "read_input_tsv_data")

  tsv_data <- readr::read_tsv(tsv_fn, col_types = "cciiccccccc") %>%
    ## Remove unexpected columns in case any are there (eager3+)
    dplyr::select(.data$Sample_Name, .data$Library_ID, .data$Strandedness, .data$UDG_Treatment, .data$R1, .data$R2, .data$BAM)

  merged_after_dedup <- tsv_data %>%
    ## Keep only number of libraries, their construction and udg_treatment
    dplyr::select(.data$Sample_Name, .data$Library_ID, .data$Strandedness, .data$UDG_Treatment) %>%
    ## Remove non-unique rows to exclude same library over multiple lanes from the count.
    dplyr::distinct() %>%
    ## Get number of libraries per sample/strandedness/UDG combination (merged_bams/initial)
    dplyr::group_by(.data$Sample_Name, .data$Strandedness, .data$UDG_Treatment) %>%
    dplyr::summarise(
      .groups='keep',
      n=dplyr::n(),
      initial_merge=ifelse(.data$n>1, T,F)
    ) %>%
    dplyr::select(-.data$n)

  merged_after_trimming <- merged_after_dedup %>%
    ## Get number of libraries per sample/strandedness/UDG combination (merged_bams/initial)
    dplyr::group_by(.data$Sample_Name, .data$Strandedness) %>%
    dplyr::summarise(
      .groups='keep',
      n=dplyr::n(),
      additional_merge=ifelse(.data$n>1, T,F)
    ) %>%
    dplyr::select(-.data$n)

  ## Get number of those combinations per sample/strandedness (merged_bams/additional)
  multiple_strandedness_per_ind <- merged_after_trimming %>%
    dplyr::group_by(.data$Sample_Name) %>%
    dplyr::summarise(
      .groups='keep',
      n=dplyr::n(),
      strandedness_clash=ifelse(.data$n>1, T,F)
    ) %>%
    dplyr::select(-.data$n)

  ## Get number of those combinations per sample to check if multiple strandedness libraries share a Sample_Name
  decision_tree <- dplyr::full_join(multiple_strandedness_per_ind, merged_after_trimming, by="Sample_Name") %>%
    dplyr::full_join(., merged_after_dedup, by=c("Sample_Name", "Strandedness")) %>%
    dplyr::full_join(., tsv_data, by=c("Sample_Name", "Strandedness", "UDG_Treatment")) %>%
    ## Reorder columns for convenience
    dplyr::select(.data$Sample_Name, .data$Library_ID, .data$Strandedness, .data$UDG_Treatment, .data$R1, .data$R2, .data$BAM, .data$initial_merge, .data$additional_merge, .data$strandedness_clash)

  ## Return tibble with TSV data including three columns with which merging steps are taken, if any.
  decision_tree
}

#' Infer the bam names for different merged states from a TSV parsed with \link[eagerR]{read_input_tsv_data}
#'
#' @param parsed_tsv Tibble. A tibble generated using \link[eagerR]{read_input_tsv_data}.
#' @param run_trim_bam logical. Was `--run_trim_bam` used in the eager run?
#'
#' @return A Tibble containing two additional columns with the names of the merged bams for a sample after each merging step.
#' @export
infer_merged_bam_names <- function(parsed_tsv, run_trim_bam = NA) {
  ## Validate run_trim_bam input
  validate_user_input(run_trim_bam, "infer_merged_bam_names", "run_trim_bam", c(T,F))

  ## BUILT WITH 2.4.3+ in mind. names might need changing for previous/future versions

  ## Technically, a user might have skipped deduplication and bams might have yet another name, but I will look into that another day.
  output <- parsed_tsv %>%
    dplyr::mutate(
      ## Bam name of initial merge
      initial_merge_bam = dplyr::case_when(
        ## No merging, no trimming BOB001.A0101.SG1_rmdup.bam
        initial_merge == F & run_trim_bam == F ~ paste0(Library_ID, "_rmdup.bam"),
        ## No merging, only trimming  BOB001.A0101.SG1.trimmed_doublestrand.bam
        initial_merge == F & run_trim_bam == T ~ paste0(Library_ID, ".trimmed.bam"),
        ## Library merging, but no trimming LTN001_udghalf_libmerged_rmdup.bam
        initial_merge == T & run_trim_bam == F ~ paste0(Sample_Name, "_udg", UDG_Treatment,"_libmerged_rmdup.bam"),
        ## Library merging, and also trimming VIX001_libmerged.trimmed_doublestrand.bam
        initial_merge == T & run_trim_bam == T ~ paste0(Sample_Name, "_libmerged.trimmed.bam"),
        TRUE ~ NA_character_ ## Catch-all clause
      ),

      ## Bam from additional merge
      additional_merge_bam = dplyr::case_when(
        ## If no additional merging, then this is same as initial_merge bam
        additional_merge == F ~ initial_merge_bam,
        ## If merge of UDG treatments needed, check for ssDNA dsDNA clash. If so, set to NA for safety.
        additional_merge == T & strandedness_clash == T ~ NA_character_,
        ## If merge of UDG treatments needed, trimmed/untrimmed status does not matter for naming
        additional_merge == T & strandedness_clash == F ~ paste0(Sample_Name, "_libmerged_add.bam"),
        TRUE ~ NA_character_ ## Catch-all clause
      ),

      ## Sexdet bam is renamed in sexdeterrmine_prep, so gotta account for that
      sexdet_bam = stringr::str_replace(.data$additional_merge_bam, ".bam$", paste0("_", .data$Strandedness,"strand.bam"))
    )

  output
}
