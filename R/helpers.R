check_file_exists <- function(tested_file, function_name) {
  if (!file.exists(tested_file)) {
    stop(paste0("[",function_name,"()]: File '", tested_file, "' not found."))
  }
}

validate_user_input <- function(tested_input, function_name, param_name, valid_entries) {
  if (! tested_input %in% valid_entries) {
    stop(paste0("[", function_name, "()]: Parameter ", param_name, " must be one of: ", paste0(valid_entries, collapse=","),". Invalid value provided:'", tested_input, "'"))
  }
}

standardise_column_names <- function(col_name, prefix = '') {
  ## Replace spaces with underscores and make all letters lowercase
  new_col_name <- paste(prefix, gsub(" ", "_", col_name) %>% tolower(), sep="_")
  new_col_name
}
