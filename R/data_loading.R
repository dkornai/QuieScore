##### LOADING DATA #####

# GENE LISTS
# full list
data(downregulated_common)
data(upregulated_common)
genelist_full_name <- list(upregulated_common, downregulated_common)

# full ENSG list
data(downregulated_common_ENSG)
data(upregulated_common_ENSG)
genelist_full_ENSG <- list(upregulated_common_ENSG, downregulated_common_ENSG)

# compact list
data(downregulated_core)
data(upregulated_core)
genelist_core_name <- list(downregulated_core, upregulated_core)

# compact ENSG list
data(downregulated_core_ENSG)
data(upregulated_core_ENSG)
genelist_core_ENSG <- list(downregulated_core_ENSG, upregulated_core_ENSG)

# gene lists for quiescence sub types
data(qlist_subt_name)
data(qlist_subt_ENSG)

# EXPRESSION DATA FROM TCGA
data(TCGA_manygenes_ensg)
data(TCGA_manygenes_name)

# cancer type list
data(cancertypes)
