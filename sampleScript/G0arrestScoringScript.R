## Install library:
library(devtools)
install_github("https://github.com/dkornai/QuieScore") 

# Load library:
library(QuieScore)

# Load sample expression matrix (rows are genes, columns are samples):
load("mat.expr.RData")
head(mat.expr)

## Process the data (need to specify cancer type, gene identifier and whether the data is log2 transformed):
processedData <- processInput(mat.expr, cancer_type = "ESCA", 
                           gene_naming = "name", log_transformed=TRUE)

# Calculate G0 arrest scores:
G0scores <- QuiescenceScore(processedData)
# Display scores:
G0scores
