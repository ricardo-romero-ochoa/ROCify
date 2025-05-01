#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess

def run_roc_analysis(genes):
    genes_str = ','.join(f'"{gene}"' for gene in genes)
    r_script = f'''
    library(pROC); library(ggplot2); library(dplyr); library(tidyr)

    data <- read.csv("expr_data.csv")

    roc_list <- list()
    auc_results <- data.frame()

    for (gene in c({genes_str})) {{
        roc_obj <- roc(data$Diagnosis, data[[gene]], levels=c(0,1), direction="<", na.action=na.omit)
        auc_ci <- ci.auc(roc_obj)
        auc_results <- rbind(auc_results, data.frame(
            Gene=gene, AUC=round(auc(roc_obj),3),
            CI_95=paste(round(auc_ci[1],3),round(auc_ci[3],3),sep="-")
        ))
        roc_list[[gene]] <- roc_obj
    }}

    png("outputs/roc_curves.png", width=800, height=600)
    plot(roc_list[[1]], col=rainbow(length(roc_list))[1], legacy.axes=TRUE)
    if(length(roc_list) > 1) {{
        for(i in 2:length(roc_list)) plot(roc_list[[i]], col=rainbow(length(roc_list))[i], add=TRUE)
    }}
    dev.off()

    write.csv(auc_results, "outputs/auc_results.csv", row.names=FALSE)

    cutoff_results <- do.call(rbind, lapply(names(roc_list), function(gene) {{
        optimal <- coords(roc_list[[gene]], "best", ret=c("threshold","specificity","sensitivity","accuracy"))
        data.frame(Gene=gene, Cutoff=round(optimal["threshold"],2), Sensitivity=round(optimal["sensitivity"],2),
                   Specificity=round(optimal["specificity"],2), Accuracy=round(optimal["accuracy"],2))
    }}))

    write.csv(cutoff_results, "outputs/optimal_cutoffs.csv", row.names=FALSE)
    '''
    with open("roc_analysis.R", "w") as file:
        file.write(r_script)

    subprocess.run(["Rscript", "roc_analysis.R"])

