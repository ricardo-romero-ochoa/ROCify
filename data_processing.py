#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess

def run_data_pipeline(tcga_project, genes):
    genes_str = ','.join(f'"{gene}"' for gene in genes)
    r_script = f'''
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(edgeR)
    library(dplyr)

    query <- GDCquery(
      project = "{tcga_project}",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query)
    data <- GDCprepare(query)

    count_matrix <- assays(data)$unstranded
    rownames(count_matrix) <- rowData(data)$gene_name
    count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

    dge <- DGEList(counts = count_matrix)
    dge <- calcNormFactors(dge)
    logcpm <- cpm(dge, log=TRUE, prior.count=3)

    expr_subset <- logcpm[c({genes_str}), ] %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Sample_ID")
    expr_subset <- expr_subset %>% mutate(Diagnosis = ifelse(substr(Sample_ID, 14, 15) == "01", 1, ifelse(substr(Sample_ID, 14, 15) == "11", 0, NA))) %>% filter(!is.na(Diagnosis))
    write.csv(expr_subset, "expr_data.csv", row.names = FALSE)
    '''
    with open("process_data.R", "w") as file:
        file.write(r_script)

    subprocess.run(["Rscript", "process_data.R"])

