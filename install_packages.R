if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "edgeR"))
install.packages(c("pROC", "ggplot2", "dplyr", "tidyr")), repos = "http://cran.us.r-project.org")
