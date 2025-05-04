import subprocess
import shutil
import os

def run_data_pipeline(tcga_project, genes):
    genes_str = ','.join(f'"{gene}"' for gene in genes)

    # Create install_packages.R call to install required R libraries
    rscript_path = shutil.which("Rscript") or "C:\\Program Files\\R\\R-4.5.0\\bin\\Rscript.exe"
    if not os.path.exists(rscript_path):
        raise FileNotFoundError("Rscript not found. Please ensure R is installed and in your PATH.")

    subprocess.run([rscript_path, "install_packages.R"], check=True)

    r_script = f'''
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(edgeR)
    library(dplyr)
    library(tibble)

    # Step 1: Download and prepare TCGA data
    query <- GDCquery(
      project = "{tcga_project}",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query)
    data <- GDCprepare(query)

    # Step 2: Normalize using edgeR
    count_matrix <- assays(data)$unstranded
    rownames(count_matrix) <- rowData(data)$gene_name
    count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

    dge <- DGEList(counts = count_matrix)
    dge <- calcNormFactors(dge)
    logcpm <- cpm(dge, log=TRUE, prior.count=3)

    # Step 3: Subset safely and transpose (handles 1 or more genes)
    genes_to_use <- c({genes_str})
    expr_matrix <- logcpm[genes_to_use, , drop = FALSE]
    expr_subset <- as.data.frame(t(expr_matrix))
    expr_subset <- expr_subset %>% tibble::rownames_to_column("Sample_ID")

    # Step 4: Add diagnosis from sample ID
    expr_subset <- expr_subset %>%
      mutate(Diagnosis = ifelse(substr(Sample_ID, 14, 15) == "01", 1,
                         ifelse(substr(Sample_ID, 14, 15) == "11", 0, NA))) %>%
      filter(!is.na(Diagnosis))

    # Step 5: Export to CSV
    write.csv(expr_subset, "expr_data.csv", row.names = FALSE)
    '''

    with open("process_data.R", "w") as f:
        f.write(r_script)

    print("Running TCGA data preparation R script...")
    subprocess.run([rscript_path, "process_data.R"], check=True)
    os.remove("process_data.R")





