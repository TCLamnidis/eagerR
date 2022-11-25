#' Read nf-core/eager Sex Determination results JSON
#'
#' @param sexdet_json character. The path to the input sex determination result json
#'
#' @return A Tibble containing the sex determination results.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_sexdet_json <- function(sexdet_json) {
  ## End execution if input file is not found
  check_file_exists(sexdet_json, "read_sexdet_json")

  sexdet_json_data <- jsonlite::read_json(sexdet_json, simplifyVector=T) ## Raw JSON data
  sexdet_results <- do.call(rbind, sexdet_json_data[-1]) %>%
    tibble::as_tibble(rownames="input_bam") %>%
    ## Standardise column names cause sexdeterrmine header names are a mess >.< Looking at you, past-me...
    dplyr::rename_with(., ~standardise_column_names(..1, prefix="sexdet")) %>%
    ## Flatten columns from single-value lists to value columns
    ## Unnest after renaming so column names are more likely to stay consistent through time :crossed_fingers:
    tidyr::unnest(
      cols = c(
        `.data$sexdet_snps_autosomal`,
        `.data$sexdet_xsnps`,
        `.data$sexdet_ysnps`,
        `.data$sexdet_nr_aut`,
        `.data$sexdet_nrx`,
        `.data$sexdet_nry`,
        `.data$sexdet_ratex`,
        `.data$sexdet_ratey`,
        `.data$sexdet_rateerrx`,
        `.data$sexdet_rateerry`
      )
    )
  sexdet_results
}

#' Read nf-core/eager nuclear contamination results JSON
#'
#' @param angsd_cont_json character. The path to the input nuclear contamination result json
#'
#' @return A Tibble containing the per-library nuclear contamination results.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_angsd_cont_json <- function(angsd_cont_json) {
  ## End execution if input file is not found
  check_file_exists(angsd_cont_json, "read_angsd_cont_json")

  angsd_cont_json_data <- jsonlite::read_json(angsd_cont_json, simplifyVector=T) ## Raw JSON data
  angsd_results <- do.call(cbind, angsd_cont_json_data['data']) %>%
    tibble::as_tibble(rownames="library_id") %>%
    ## All values are in a list column called 'data'. To get multiple columns I need to transpose this list, turn it to a dataframe and unnest it to multiple list columns.
    dplyr::mutate(r = purrr::map(.data$data, ~ data.frame(t(.)))) %>%
    tidyr::unnest(.data$r) %>%
    ## Drop the fully nested column
    dplyr::select(-.data$data) %>%
    ## Standardise column names
    dplyr::rename_with(., ~standardise_column_names(..1, prefix="angsd")) %>%
    ## Finally unnest the final nested columns to get a standard non-listed tibble
    ## Unnest after renaming so column names are more likely to stay consistent through time :crossed_fingers:
    tidyr::unnest(
      cols = c(
        `.data$angsd_num_snps`,
        `.data$angsd_method1_mom_estimate`,
        `.data$angsd_method1_mom_se`,
        `.data$angsd_method1_ml_estimate`,
        `.data$angsd_method1_ml_se`,
        `.data$angsd_method2_mom_estimate`,
        `.data$angsd_method2_mom_se`,
        `.data$angsd_method2_ml_estimate`,
        `.data$angsd_method2_ml_se`
      )
    ) %>%
    ## Make N/As into native R NA
    dplyr::na_if("N/A")
  angsd_results
}

#' Read nf-core/eager EIGENSTRAT SNP coverage results JSON
#'
#' @param snp_cov_json character. The path to the input eigenstrat snp coverage result json
#'
#' @return A Tibble containing the per-library eigenstrat snp coverage results.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_snp_coverage_json <- function(snp_cov_json) {
  ## End execution if input file is not found
  check_file_exists(snp_cov_json, "read_snp_coverage_json")

  snp_coverage_json_data <- jsonlite::read_json(snp_cov_json, simplifyVector=T) ## Raw JSON data
  snp_coverage_results <- do.call(cbind, snp_coverage_json_data['data']) %>%
    tibble::as_tibble(rownames="sample_name") %>%
    ## All values are in a list column called 'data'. To get multiple columns I need to transpose this list, turn it to a dataframe and unnest it to multiple list columns.
    dplyr::mutate(r = purrr::map(.data$data, ~ data.frame(t(.)))) %>%
    tidyr::unnest(.data$r) %>%
    ## Drop the fully nested column
    dplyr::select(-.data$data) %>%
    ## Standardise column names
    dplyr::rename_with(., ~standardise_column_names(..1, prefix="snpcov")) %>%
    ## Finally unnest the final nested columns to get a standard non-listed tibble
    ## Unnest after renaming so column names are more likely to stay consistent through time :crossed_fingers:
    tidyr::unnest(
      cols = c(
        `.data$snpcov_covered_snps`,
        `.data$snpcov_total_snps`
        )
      )

  snp_coverage_results
}

#' Read nf-core/eager DamageProfiler results JSON
#'
#' @param dmgprof_json character. The path to the input DamageProfiler result json for the specified library.
#' @param library_id character. The Library_ID of the library the results are based on.
#'
#' @return A Tibble containing the DamageProfiler results for the specified library.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_damageprofiler_json <- function(dmgprof_json, library_id = '') {
  ## Damageprofiler results are in files for each library, so parsing the JSON is done individually,
  ##   and can be made into a tibble with purrr::map_dfr/rbind

  ## End execution if input file is not found
  check_file_exists(dmgprof_json, "read_damageprofiler_results")

  ## lib_id MUST be provided
  if (library_id == '' ) {
    stop(paste0("[read_damageprofiler_results()]: No library id provided. Please provide a valid library id."))
  }

  dmgprof_json_data <- jsonlite::read_json(dmgprof_json, simplifyVector=F) ## Raw JSON data

  ## Calculate number of reads used in damage calculation for weighted summing across libraries later.
  ##  Calculated as the sum of all fw and rv observations in the length distribution
  n_reads <- dmgprof_json_data$lendist_fw %>% unlist() %>% sum()
  n_reads <- n_reads + dmgprof_json_data$lendist_rv %>% unlist() %>% sum()

  ## Create a tibble with all the listed values from the JSON and the library ID and number of reads
  dmgprof_results <- tibble::tibble(
      ## Add library_id and number of reads to table.
      library_id = library_id,
      num_reads = n_reads,
      lendist_fw = dmgprof_json_data['lendist_fw'],
      lendist_rv = dmgprof_json_data['lendist_rv'],
      summary_stats = dmgprof_json_data['summary_stats'],
      dmg_5p = dmgprof_json_data['dmg_5p'],
      dmg_3p = dmgprof_json_data['dmg_3p']
    ) %>%
    ## Unnest dmg columns
    tidyr::unnest(cols=c(.data$dmg_5p,.data$dmg_3p)) %>%
    dplyr::mutate(x=dplyr::row_number()) %>%
    ## Use row number to infer the bp
    tidyr::pivot_wider(values_from=c(.data$dmg_5p, .data$dmg_3p), names_from=.data$x) %>%
    ## Finally, standardise the column names
    dplyr::rename_with(., ~standardise_column_names(..1, prefix="damage"))

  dmgprof_results
}

#' Read all Damageprofiles JSONs in a DamageProfiler results directory into a tibble.
#'
#' @param result_dir character. The path to the root DamageProfiler results directory for an eager run.
#'
#' @return A Tibble containing the JSON results (compiled using \link[eagerR]{read_damageprofiler_json}) for all
#'     libraries in the eager run.
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_damageprofiler_jsons_from_dir <- function(result_dir) {
  ## Check the dir exists
  check_file_exists(result_dir, "read_damageprofiler_json_from_dir")

  ## Find all files in result dir that end in .json and read into a tibble column
  dmg_json_list <- list.files(result_dir, pattern="\\.json$", full.names = T, recursive=T) %>%
  tibble::as_tibble() %>%
  ## Infer the library ID from the result directory name. First split at "_rmdup", then take the basename of the first element.
    ## This will not work if the user's eager-results directory path includes '_rmdup'. This should be a very niche case... :crossed_fingers:
  dplyr::mutate(lib_id=basename(stringr::str_split(.data$value, "_rmdup", simplify=TRUE)[,1]))

  damage_results <- dmg_json_list %>%
    purrr::pmap_dfr(., ~read_damageprofiler_json(..1, ..2))
  damage_results
}

#' Read nf-core/eager endorspy results JSON
#'
#' @param endorspy_json character. The path to the input endorspy result json for the specified library.
#' @param library_id character. The Library_ID of the library the results are based on.
#'
#' @return A Tibble containing the endorspy results for the specified library.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_endorspy_json <- function(endorspy_json, library_id = '') {
  ## Damageprofiler results are in files for each library, so parsing the JSON is done individually,
  ##   and can be made into a tibble with purrr::map_dfr/rbind

  ## End execution if input file is not found
  check_file_exists(endorspy_json, "read_endorspy_json")

  ## lib_id MUST be provided
  if (library_id == '' ) {
    stop(paste0("[read_endorspy_json()]: No library id provided. Please provide a valid library id."))
  }

  endorspy_json_data <- jsonlite::read_json(endorspy_json, simplifyVector=F) ## Raw JSON data

  endorspy_results <- do.call(cbind, endorspy_json_data['data']) %>%
    tibble::as_tibble() %>%
    ## All values are in a list column called 'data'. To get multiple columns I need to transpose this list, turn it to a dataframe and unnest it to multiple list columns.
    dplyr::mutate(r = purrr::map(.data$data, ~ data.frame(t(.)))) %>%
    tidyr::unnest(.data$r) %>%
    ## Drop the fully nested column
    dplyr::select(-.data$data) %>%
    ## Add library ID
    dplyr::mutate(
      library_id = library_id
    ) %>%
    ## Standardise column names
    dplyr::rename_with(., ~standardise_column_names(..1, prefix="endorspy")) %>%
    ## Finally unnest the final nested columns to get a standard non-listed tibble
    ## Unnest after renaming so column names are more likely to stay consistent through time :crossed_fingers:
    tidyr::unnest(
      cols = dplyr::everything() ## Should work even if no '_post' column exists
    )
}

#' Read all endorspy JSONs in an endorspy results directory into a tibble.
#'
#' @param result_dir character. The path to the endorspy results directory for an eager run.
#'
#' @return A Tibble containing the JSON results (compiled using \link[eagerR]{read_endorspy_json}) for all
#'     libraries in the eager run.
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
read_endorspy_jsons_from_dir <- function(result_dir) {
  ## Check the dir exists
  check_file_exists(result_dir, "read_endorspy_jsons_from_dir")

  ## Find all files in result dir that end in .json and read into a tibble column
  endorspy_json_list <- list.files(result_dir, pattern="_endogenous_dna_mqc\\.json$", full.names = T, recursive=T) %>%
    tibble::as_tibble() %>%
    ## Infer the library ID from the result file name. Take the basename, and then remove the suffix.
    dplyr::mutate(lib_id=basename(.data$value) %>% sub("_endogenous_dna_mqc.json$", '', .))

  endorspy_results <- endorspy_json_list %>%
    purrr::pmap_dfr(., ~read_endorspy_json(..1, ..2))

  endorspy_results
}
