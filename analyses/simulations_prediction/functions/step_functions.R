library(recipes)

step_aggregate <- function(
    recipe, 
    ...,
    role = NA,
    trained = FALSE,
    agg_level = "Genus",
    tax_table = NULL,
    index = NULL,
    skip = FALSE,
    id = rand_id("aggregate")
){
    terms <- ellipse_check(...)
    add_step(
        recipe,
        step_aggregate_new(
            terms = terms, 
            role = role, 
            trained = trained, 
            agg_level = agg_level,
            tax_table = tax_table,
            index = index,
            skip = skip,
            id = id
        )
    )
}

step_aggregate_new <- function(terms, role, trained, agg_level, tax_table, index, skip, id){
    step(
        subclass = "aggregate",
        terms = terms, 
        role = role, 
        trained = trained,
        agg_level = agg_level,
        tax_table = tax_table,
        index = index,
        skip = skip,
        id = id
    )
}

# prep method
# This method would then 
prep.step_aggregate <- function(x, training, info = NULL, ...){
    col_names <- terms_select(terms = x$terms, info = info)
    
    if (length(col_names) != ncol(training)){
        rlang::abort("Transformation must be done across all variables")
    }
    # generate the index matrix from taxtable 
    index_mat <- taxtab2A(tax = x$tax_table, agg_level = x$agg_level, full = TRUE)
    index <- vector(mode = "list", length = ncol(index_mat))
    for (i in seq_along(colnames(index_mat))){
        index[[i]] <- which(index_mat[,i] == 1) %>% unname(which(testing[,1] == 1))
    }
    names(index) <- colnames(index_mat)
    step_aggregate_new(
        terms = terms, 
        trained = TRUE, 
        role = x$role, 
        agg_level = x$agg_level,
        tax_table = x$tax_table,
        index = index, 
        skip = x$skip,
        id = x$id
    )
}

# bake method
bake.step_aggregate <- function(object, new_data, ...){
    new_data <- map_dfc(object$index, ~{
        rowSums(new_data[,.x])
    })
    new_data
}
