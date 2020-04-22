library(tidyverse)
library(phyloseq)
library(recipes)

#' @title Generating the a matrix 
#TODO: Add phyloseq object support  
generate_a_matrix <- function(tab, taxlevel, drop_unknown = FALSE){
  # processing the tab data frame 
  if (class(tab) == "taxonomyTable"){
    
  }
  tab_tib <- tab %>% as("matrix") %>% as.data.frame() %>% rownames_to_column(var = "tax_id") %>% 
    as_tibble() %>% select(c(tax_id, !!quo(taxlevel))) %>% 
    mutate_at(vars(-tax_id), function(.x){replace_na(as.character(.x), "Unknown")})
  # make dummy variables 
  dummy <- tab_tib %>% recipe(tax_id ~ .) %>% step_dummy(-tax_id) %>% prep(training = tab_tib) %>% 
    bake(new_data = tab_tib)
  if (drop_unknown == T){
    dummy <- dummy %>% select(-ends_with("Unknown"))
  }
  return(dummy)
}



