##### BACKGROUND FUNCTIONS #####

# performs batch normalization to the TCGA data to make results comparable
library(sva)

batNorm <- function(expr_user, expr_TCGA){
  user_data <- t(expr_user)
  # select genes which are reported in the TCGA and in the expression data of interest
  TCGA_genes <- colnames(expr_TCGA)
  Userdata_genes <- colnames(user_data)
  Common_genes <- TCGA_genes[TCGA_genes %in% Userdata_genes]
  expr_TCGA <- expr_TCGA[,colnames(expr_TCGA) %in% Common_genes]
  user_data <- user_data[,colnames(user_data) %in% Common_genes]
  combined_data <- rbind(expr_TCGA, user_data)
  # save batch annotation
  batch <- c(rep("TCGA", nrow(expr_TCGA)), rep("External_dataset", nrow(user_data)))
  batch <- data.frame(batch)
  rownames(batch) <- c(rownames(expr_TCGA),rownames(user_data))
  # setup correction
  combined_data <- data.matrix(combined_data)
  combined_data <- t(combined_data)
  modcombat = model.matrix(~1, data=batch)
  batch = batch$batch
  # apply the Combat function
  scaled_data = sva::ComBat(dat=combined_data, batch=batch, mod=modcombat, par.prior=TRUE)
  scaled_data <- data.frame(scaled_data)

  return(scaled_data)
}

# percentile calculation
percentile <- function(compare_to, compare) {
  percentile_vec <- vector(length = length(compare))
  pre_sorted <- sort(compare_to)
  for (i in 1:length(compare)){
    numerator <- length(pre_sorted[compare_to < compare[i]])
    denominator <- length(compare_to)
    percentile_vec[i] <- round(numerator/denominator,3)*100
  }

  return(percentile_vec)
}

# quiescence score calculation for a given data set
library(GSVA)

qScore <- function(input_expr, method_eval, genelist, parallel){
  # filter expression data so that only relevant genes are included
  user_expr <- as.matrix(input_expr[which(rownames(input_expr) %in% unlist(genelist, recursive = TRUE)), ])

  # warn the user that genes used for inferring the given quiescence type are not all present
  if(length(which(rownames(user_expr) %in% unlist(genelist[1]))) < 2){
    cat("WARNING: Less than 2 upregulated genes required present in the user data!\n") }
  if(length(which(rownames(user_expr) %in% unlist(genelist[2]))) < 2){
    cat("WARNING: Less than 2 downregulated genes required present in the user data!\n") }

  # calculate quiescence scores with one of the methods
  if(method_eval == "zscore"){
    q_score <- GSVA::gsva(user_expr, genelist, method = "zscore", parallel.sz=parallel) }
  else if (method_eval == "ssgsea"){
    q_score <- GSVA::gsva(user_expr, genelist, method = "ssgsea",  parallel.sz=parallel) }
  else if (method_eval == "gsva"){
    q_score <- GSVA::gsva(user_expr, genelist, mx.diff= TRUE,  parallel.sz=parallel) }

  # formatting and cleanup
  q_score <- data.frame(t(q_score))
  q_score$score <- q_score$X1 - q_score$X2
  q_score$X1 <- NULL; q_score$X2 <- NULL

  return(q_score)
}

# calculation of batch corrected quiescence score, and the percentile figure of the score. Function can be configured to return values for the user data or TCGA data
qScore_percentile <- function(input_expr, method_eval, genelist, batch_vector, user_or_TCGA, parallel){
  # calculate q_scores for all samples
  qscore <- qScore(input_expr, method_eval, genelist, parallel)
  # split results
  user_qscore <- qscore[batch_vector == TRUE,]
  TCGA_qscore <- qscore[batch_vector == FALSE,]
  # calculate percentiles
  percentile_user_vec <- percentile(TCGA_qscore, user_qscore)
  percentile_TCGA_vec <- percentile(TCGA_qscore, TCGA_qscore)
  # make output the user or TCGA results, depending on "user_or_TCGA"
  if (user_or_TCGA == "user"){
    output <- matrix(nrow = length(user_qscore), ncol = 2)
    output[,1] <- user_qscore
    output[,2] <- percentile_user_vec
  }
  else if (user_or_TCGA == "TCGA"){
    output <- matrix(nrow = length(TCGA_qscore), ncol = 2)
    output[,1] <- TCGA_qscore
    output[,2] <- percentile_TCGA_vec
  }

  return(data.frame(output))
}

# calculation of all the subtype specific quiescence scores and percentiles
qScore_subtype <- function(input_bcorr_expr, method_eval, q_subtype_genes, batch_vector, parallel){
  # set up output table
  qscore_table <- matrix(nrow = sum(batch_vector), ncol = 10)
  colnames(qscore_table) <- c("q_cdk", "q_cdk_percentile", "q_cont", "q_cont_percentile","q_mek",
                              "q_mek_percentile","q_ser", "q_ser_percentile","q_spn", "q_spn_percentile")

  qtypes <- c("CDK4/6 inhibition", "Contact inhibition", "MEK inhibition", "Serum starvation", "Spontaneous")
  qscore_location <- seq(1, 9, by = 2)

  # for each subtype
  for (i in 1:5){
    cat(paste0(" Subtype calculations initiated for: ", qtypes[i], "\n"))
    # get gene list specific to subtype
    gene_list_subtype <- unlist(q_subtype_genes[i], recursive = FALSE)
    # calculate score and add to output
    q_score_subtype <- qScore_percentile(input_bcorr_expr, method_eval, gene_list_subtype, batch_vector, "user", parallel)
    qscore_table[,qscore_location[i]]   <- q_score_subtype[,1]
    qscore_table[,qscore_location[i]+1] <- q_score_subtype[,2]
  }

  return(qscore_table)
}

# checks the input to the processInput function, and returns an error message if the parameters are of the incorrect type or value
checkInput_processInput <- function(input_object, cancer_type, gene_naming, log_transformed){
  ## check that arguments are in the expected format, and stop the program if required
  error_msg <- ""
  # check if input_object is dataframe
  if (is.data.frame(input_object) == FALSE){
    error_msg <- c(error_msg, "ERROR: the samples with gene expression data must be supplied as a dataframe\n")
  }
  # check if cancer_type is valid
  if (sum(cancertypes %in% cancer_type) != 1){
    error_msg <- c(error_msg, "ERROR: please use the 3 or 4 letter code of one of the 31 supported cancer types\n")
  }
  # check if gene_naming is either "name" or "ensg"
  if (gene_naming != "name" & gene_naming != "ensg"){
    error_msg <- c(error_msg,"ERROR: gene_naming must be 'name' or 'ensg' \n depending on wether the genes are annotated using their names or ensemble codes\n")
  }
  # check if log_transformed is either true or false
  if (is.logical(log_transformed) == FALSE){
    error_msg <- c(error_msg, "ERROR: log_transformed must be TRUE or FALSE \n depending on wether the expression data has been log2 transformed already\n")
  }

  return(error_msg)
}

# checks the input to the main QuiescenceScore function, and returns an error message if the parameters are of the incorrect type or value
checkInput_QuiescenceScore <- function(input_data, genelist_eval, method_eval, subtype_score, vis, parallel){
  error_msg <- ""
  # check that the input data is in the expected format
  if (is.list(input_data) == FALSE | length(input_data) != 3){
    error_msg <- c(error_msg,"ERROR: data supplied to function is in incorrect format \n please only supply data generated using the 'processInput' function \n")
  }
  # check that the gene list used in the downstream analysis is in a valid format
  if (genelist_eval != "full" & genelist_eval != "core"){
    error_msg <- c(error_msg,"ERROR: genelist_eval must be 'full' or 'core' \n depending on whether the queiescence scores should be assesed using the full or core gene set \n")
  }
  # check that method_eval is valid
  if (method_eval != "zscore" & method_eval != "ssgsea" & method_eval != "gsva"){
    error_msg <- c(error_msg,"ERROR: method_eval must be 'zscore' 'gsva' or 'ssgsea' \n depending on the intended mode of calculating quiescence scores\n")
  }
  # check that subtype_score is in the expected format
  if (is.logical(subtype_score) == FALSE){
    error_msg <- c(error_msg,"ERROR: log_transformed must be TRUE or FALSE \n depending on wether the quiescence subtype of the data should be analysed\n")
  }
  # check that the parallel value is a valid integer
  if (parallel%%1 != 0){
    error_msg <- c(error_msg,"ERROR: parallel must be set to a valid integer less than the number of CPU cores in the system\n")
  }


  return(error_msg)
}

# plots a PCA of the gene expression results, and colors the points by their batch (user or TCGA), and their percentile value
library(ggplot2)
library(ggnewscale)

plot_PCA <- function(dataset_name, input_bcorr_expr, method_eval, gene_list, batch_vector, parallel){
  cat("Generating plots... \n")
  # get the quiescence scores and percentiles of the user and TCGA data
  qscore_TCGA <- invisible(qScore_percentile(input_bcorr_expr, method_eval, gene_list, batch_vector, "TCGA", parallel))
  qscore_user <- invisible(qScore_percentile(input_bcorr_expr, method_eval, gene_list, batch_vector, "user", parallel))

  # this vector contains the combined quiescence scores
  combined_q <- c(as.vector(qscore_TCGA[,1]), as.vector(qscore_user[,1]))
  # this vector contains the combined percentile scores
  combined_p <- c(as.vector(qscore_TCGA[,2]), as.vector(qscore_user[,2]))
  annotation_table <- as.data.frame(cbind(batch_vector, combined_q, combined_p))

  # filter table to only include quiescence relevant genes
  combined_results_table <- t(input_bcorr_expr)
  combined_results_table <- combined_results_table[,colnames(combined_results_table) %in% unlist(gene_list, recursive = TRUE)]

  # perform PCA
  pca_vis <- princomp(combined_results_table)

  # extract the amount of variance explained by the first 2 principal components
  pca_eig <- (pca_vis$sdev)^2
  pca_component_1 <- round(((pca_eig[1]/sum(pca_eig))*100), digits = 2)
  pca_component_2 <- round(((pca_eig[2]/sum(pca_eig))*100), digits = 2)

  # filter PCA results to only include coordinates
  pca_vis <- as.data.frame(pca_vis[["scores"]])[,1:2]
  colnames(pca_vis) <- c("c1", "c2")

  # split data frame in preparation for visualization
  pca_vis_tcga <- pca_vis[batch_vector == FALSE,]
  pca_vis_user <- pca_vis[batch_vector == TRUE,]
  data_tcga <- annotation_table[batch_vector == FALSE,]
  data_user <- annotation_table[batch_vector == TRUE,]

  # plot
  plot_vis_cord <- ggplot2::ggplot(pca_vis) +
    # plotting of TCGA samples
    ggplot2::geom_point(data = pca_vis_tcga, ggplot2::aes(x = c1, y = c2, colour = data_tcga$combined_p), alpha = 0.65, shape = 21, size = 2.5) +
    ggplot2::scale_color_gradient("quiescence percentile\nof TCGA samples", low = "blue2", high = "coral3") +
    # plotting of user samples
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data = pca_vis_user, ggplot2::aes(x = c1, y = c2, colour = data_user$combined_p), shape = 18, size = 3.5) +
    ggplot2::scale_color_gradient("quiescence percentile\nof user samples",low = "green", high = "purple") +
    # labels showing percentage explained
    ggplot2::xlab(paste0("PC1 (", pca_component_1, "%)")) +
    ggplot2::ylab(paste0("PC2 (", pca_component_2, "%)")) +
    ggplot2::ggtitle(paste0(dataset_name, ": PCA plot of gene expression data"))

  return(plot_vis_cord)
}
