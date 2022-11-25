#' Check that a file exists, else print an error and exit.
#'
#' @param tested_file character. The path to the file to be tested.
#' @param function_name character. The name of the function calling this one. To be printed in the error for debugging.
#'
#' @noRd
check_file_exists <- function(tested_file, function_name) {
  if (!file.exists(tested_file)) {
    stop(paste0("[",function_name,"()]: File '", tested_file, "' not found."))
  }
}

#' Validate that a user's input is among the allowed choices, else throw an error.
#'
#' @param tested_input The input given by the user.
#' @param function_name character. The name of the function calling this one. To be printed in the error for debugging.
#' @param param_name characer. The name of the paramater being validated. To be printed in the error for debugging.
#' @param valid_entries A vector of the valid inputs for the parameter.
#'
#' @noRd
validate_user_input <- function(tested_input, function_name, param_name, valid_entries) {
  if (! tested_input %in% valid_entries) {
    stop(paste0("[", function_name, "()]: Parameter ", param_name, " must be one of: ", paste0(valid_entries, collapse=","),". Invalid value provided:'", tested_input, "'"))
  }
}

#' Standardise the naming of columns
#'
#'  Will turn the given column name to lowercase, replace any spaces with underscores
#'  and add the desired prefix, if any.
#'
#' @param col_name character. The name of the column to be standardise.
#' @param prefix character. A desired prefix to be added to the column name.
#'
#' @return character. The standardised column name.
#' @noRd
standardise_column_names <- function(col_name, prefix = '') {
  ## Replace spaces with underscores and make all letters lowercase
  if (is.na(prefix)) {
    stop(paste0("[standardise_column_names()]: Column name prefix cannot be NA."))
  } else if (prefix == '') {
    new_col_name <- gsub(" ", "_", col_name) %>% tolower()
  } else {
    new_col_name <- paste(prefix, gsub(" ", "_", col_name) %>% tolower(), sep="_")
  }
  new_col_name
}
