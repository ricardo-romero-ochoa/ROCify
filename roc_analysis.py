import subprocess
import shutil
import os

def run_roc_analysis(genes):
    genes_str = ','.join(f'"{gene}"' for gene in genes)

    r_script = f'''
    library(pROC)
    library(ggplot2)
    library(dplyr)
    library(tidyr)

    data <- read.csv("expr_data.csv")

    if (!dir.exists("outputs")) dir.create("outputs")

    roc_list <- list()
    auc_results <- data.frame()

    for (gene in c({genes_str})) {{
        tryCatch({{
            if (!(gene %in% colnames(data))) {{
                message(paste("Gene not found in expr_data.csv:", gene))
                next
            }}
            roc_obj <- roc(data$Diagnosis, data[[gene]], levels=c(0,1), direction="<", na.action=na.omit)
            auc_ci <- ci.auc(roc_obj)
            auc_results <- rbind(auc_results, data.frame(
                Gene=gene,
                AUC=round(auc(roc_obj), 3),
                CI_95=paste(round(auc_ci[1], 3), round(auc_ci[3], 3), sep="-")
            ))
            roc_list[[gene]] <- roc_obj
        }}, error = function(e) {{
            message(paste("Error with gene", gene, ":", e$message))
        }})
    }}

    if (length(roc_list) > 0) {{
        png("outputs/roc_curves.png", width=800, height=600)
        plot(roc_list[[1]], col=rainbow(length(roc_list))[1], legacy.axes=TRUE, main="ROC Curves")
        if (length(roc_list) > 1) {{
            for (i in 2:length(roc_list)) {{
                plot(roc_list[[i]], col=rainbow(length(roc_list))[i], add=TRUE)
            }}
        }}
        abline(a=0, b=1, lty=2, col="gray")
        legend("bottomright", legend=names(roc_list), col=rainbow(length(roc_list)), lwd=2)
        dev.off()
    }} else {{
        message("No valid ROC curves generated.")
    }}

    write.csv(auc_results, "outputs/auc_results.csv", row.names=FALSE)

    cutoff_results <- do.call(rbind, lapply(names(roc_list), function(gene) {{
        coords_out <- coords(roc_list[[gene]], "best", ret=c("threshold","specificity","sensitivity","accuracy"))
        data.frame(
            Gene=gene,
            Cutoff=round(coords_out["threshold"], 2),
            Sensitivity=round(coords_out["sensitivity"], 2),
            Specificity=round(coords_out["specificity"], 2),
            Accuracy=round(coords_out["accuracy"], 2)
        )
    }}))

    write.csv(cutoff_results, "outputs/optimal_cutoffs.csv", row.names=FALSE)
    '''

    with open("roc_analysis.R", "w") as f:
        f.write(r_script)

    rscript_path = shutil.which("Rscript") or "C:\\Program Files\\R\\R-4.5.1\\bin\\Rscript.exe"
    if not os.path.exists(rscript_path):
        raise FileNotFoundError("Rscript not found. Please ensure R is installed and in your PATH.")

    print("Running ROC analysis R script...")
    subprocess.run([rscript_path, "roc_analysis.R"], check=True)
    # Optional cleanup:
    os.remove("roc_analysis.R")



