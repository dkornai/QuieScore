##### USER FACING FUNCTIONS #####

#' @title processInput
#' @description  This function takes user supplied data frames of gene expression data, and formats them to be suitable for the 'QuiescenceScore' function
#'
#' @param input_object dataframe containing gene expression data, columns are samples, rows are genes
#' @param cancer_type the 3 or 4 letter code of the cancer type in the submitted samples. Can only be 1 type.
#' @param gene_naming which type of naming scheme are used to specify the genes? set to 'name' for names, and 'ensg' for ENSG
#' @param log_transformed Optional Parameter, the default value is TRUE. If the expression data needs to get log2 transformed, manually specify FALSE
#'
#' @return A list formatted object that can be passed to the "QuiescenceScore" function
#' @export
#'
#' @examples
#' patient1_data_processed <- processInput(patient1_data_raw, cancer_type = "LUAD", gene_naming = "name")
#' patient2_data_processed <- processInput(patient2_data_raw, cancer_type = "BRCA", gene_naming = "ensg")
#' patient3_data_processed <- processInput(patient3_data_raw, cancer_type = "PAAD", gene_naming = "ensg", log_transformed = FALSE)

processInput <- function(input_object, cancer_type, gene_naming, log_transformed = TRUE){
  # check if the input arguments are of the correct types, and have valid values
  error_msg <- checkInput_processInput(input_object, cancer_type, gene_naming, log_transformed)
  if (length(error_msg) != 1){
    cat(error_msg)
    stop()
  }

  # implements log transformation if required
  if (log_transformed == FALSE){
    input_object <- log2(input_object + 1)
  }

  # removes rows with only NA values, and warns the user that removals have occurred
  input_gene_number <- nrow(input_object)
  complete_input_gene_names <- rownames(input_object)
    # collect rows which are not complete
  is_not_complete <- !complete.cases(input_object)
  input_object <- na.omit(input_object)

  if(nrow(input_object) < input_gene_number){
    cat("WARNING: data contains NA values for genes:\n")
      # print the names of the removed genes
    genes_with_NA <- complete_input_gene_names[is_not_complete]
    for (i in 1:length(genes_with_NA)){
      cat(genes_with_NA[i])
      cat(" \n")
    }
    cat("these genes were removed\n\n");
  }

  # check if gene names contain ENSG tags, and warn the user if the input data has been mislabeled
  if (gene_naming == "name"){
    if (sum(as.numeric(grepl("ENSG0", rownames(input_object)))) > 0){
      cat("WARNING: gene list contains ENSG gene names, but the program is expecting symbols")
    }
  }

  return(list(input_object, gene_naming, cancer_type))
}

#' @title QuiescenceScore
#' @description This function calculates general and sub-type specific quiescence scores using a variety of methods, and also generates a percentile figure which compares the score to a TCGA dataset from the same cancer type.
#'
#' @param input_data a file generated using the 'processInput' function
#' @param genelist_eval specify whether the quiescence phenotype should be assessed using the full 139 genes or a core set of 18 genes. Default is "full", manually specify if the "core" set should be used
#' @param method_eval specify the methodology used to calculate quiescence scores from the expression data. Default mode is "zscore", manually set to "gsva" or "ssgsea" if other methods should be used.
#' @param subtype_score specify whether the samples should be scored according to their quiescence subtypes. Default is TRUE, manually set to FALSE if only the main quiescence score should be calculated
#' @param vis specify whether the results of the analysis should be visualized in a PCA plot. Default is TRUE, manually set to FALSE if PCA visualizations are not needed.
#' @param parallel specify the number of CPU threads that will be used in the calculation. Also impacts memory usage, with larger values requiring larger memory. If an analysis fails due to a lack of memory, lower this value
#'
#' @return a data frame that details the general and optionally sub-type specific quiescence score of samples, and the percentile value of the score compared to a TCGA data set from the same cancer.
#' @export
#'
#' @import GSVA
#' @import sva
#' @import ggnewscale
#' @import ggplot2
#' @import ggpubr
#' @examples
#' patient1_results <- QuiescenceScore(patient1_data_processed)
#' patient2_results <- QuiescenceScore(patient2_data_processed, genelist_eval = "core")
#' patient3_results <- QuiescenceScore(patient3_data_processed, subtype_score = TRUE, vis = FALSE)
#' patient4_results <- QuiescenceScore(patient4_data_processed, subtype_score = TRUE, parallel = 8)

QuiescenceScore <- function(input_data, genelist_eval = "full", method_eval = "zscore", subtype_score = FALSE, vis = TRUE, parallel = 4){
  # check if the input arguments are of the correct types, and have valid values
  error_msg <- checkInput_QuiescenceScore(input_data, genelist_eval, method_eval, subtype_score, vis, parallel)
  if (length(error_msg) != 1){
    cat(error_msg)
    stop()
  }

  # set internal cancer type variable
  cancer_type <- unlist(input_data[3])

  # load in the required gene lists
  gene_naming <- input_data[2]
  if (gene_naming == "name"){
      # load TCGA and filter to be cancer type specific
    TCGA_expr <- TCGA_manygenes_name[TCGA_manygenes_name$CancerType %in% cancer_type,]
      # load the genes used for identifying subtypes of quiescence
    q_subtype_genes <- qlist_subt_name
      # load the full or core list for identifying general quiescence
    if (genelist_eval == "full"){ gene_list <- genelist_full_name }
    else if (genelist_eval == "core"){ gene_list <- genelist_core_name }
  }
  else if (gene_naming == "ensg"){
      # load TCGA and filter to be cancer type specific
    TCGA_expr <- TCGA_manygenes_ensg[TCGA_manygenes_ensg$CancerType %in% cancer_type,]
      # load the genes used for identifying subtypes of quiescence
    q_subtype_genes <- qlist_subt_ENSG
      # load the full or core list for identifying general quiescence
    if (genelist_eval == "full"){ gene_list <- genelist_full_ENSG }
    else if (genelist_eval == "core"){ gene_list <- genelist_core_ENSG }
  }

  # load the user's expression data into a dedicated data frame
  user_expr <- as.data.frame(input_data[1])
  # create the batch annotation vector
  user_batch <- c(rep(FALSE,nrow(TCGA_expr)),rep(TRUE,ncol(user_expr)))

  # calculation of quiescence scores, and percentiles
  cat("General Quiescence score calculation started\n")
  qscore_user_raw <- qScore(user_expr, method_eval, gene_list, parallel)

  cat("Batch correcting TCGA samples and user data...\n")
  bcorr_expr <- batNorm(user_expr, TCGA_expr)
  qscore_user_bcorr <- qScore_percentile(bcorr_expr, method_eval, gene_list, user_batch, "user", parallel)

  results <- cbind(qscore_user_raw, qscore_user_bcorr)
  rownames(results) <- colnames(user_expr)
  colnames(results) <- c("q_score_raw", "q_score_normalized", "q_percentile")
  cat("  Calculation finished\n")

  # optional scoring of quiescence sub types
  if (subtype_score == TRUE){
    qscore_subtypes <- qScore_subtype(bcorr_expr, method_eval, q_subtype_genes, user_batch, parallel)
    results <- cbind(results, qscore_subtypes)
  }

  # PCA visualization
  name <- deparse(substitute(input_data))
  if (vis == TRUE){
    plot <- plot_PCA(name, bcorr_expr, method_eval, gene_list, user_batch, parallel)
    print(plot)
  }

  cat("analysis finished\n")
  return(results)
}
